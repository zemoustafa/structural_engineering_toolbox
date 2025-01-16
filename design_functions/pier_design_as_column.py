import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sectionproperties.pre.library import concrete_rectangular_section

from concreteproperties import Concrete, ConcreteLinear, ConcreteSection, RectangularStressBlock, SteelBar, SteelElasticPlastic

from concreteproperties.design_codes import AS3600
from concreteproperties.results import MomentInteractionResults
from shapely.geometry import LineString

def __round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def __find_intersection(f_m_x, f_n_x, m_star, n_star):
    """
    Find the intersection of a line with a curve.

    Args:
        f_m_x (list of floats): x-coordinates of the curve.
        f_n_x (list of floats): y-coordinates of the curve.
        m_star (float): x-coordinate of the end point of the line.
        n_star (float): y-coordinate of the end point of the line.

    Returns:
        Point or MultiPoint: Intersection point(s) between the extended line and the curve.
    """
    # Combine x and y coordinates into a list of tuples for the curve
    curve_coords = list(zip(f_m_x, f_n_x))
    
    # Define the curve as a LineString
    curve = LineString(curve_coords)
    
    # Define the line extending from (0, 0) through (m_star, n_star)
    direction = (m_star, n_star)
    extended_line = LineString([(0, 0), (10 * direction[0], 10 * direction[1])])  # Extend by a factor of 10 or more

    # Check if the line is within the curve
    line_segment = LineString([(0, 0), direction])
    if curve.contains(line_segment):
        # If the line is within the curve, extend it and find the intersection
        intersection = curve.intersection(extended_line)
    else:
        # If the line is outside, find the intersection directly
        intersection = curve.intersection(extended_line)

    return intersection
    
def __calculate_effective_shear_depth(d, cover, v_bar_area, v_bar_cts, h_bar_dia):
    '''
    Determine effective shear depth of deep beam section

    Parameters:
    d (float): total depth of section
    cover (float): cover to reinforcement
    v_bar_area (float): area of single vertical bar
    v_bar_cts (float): spacing of vertical reinforcement
    h_bar_dia (float): diameter of horziontal bar or ligs

    Returns:
    v_uc (float): concrete contribution to shear strength
    v_us (float): steel contribution to shear strength

    '''
    flexural_tension_zone = d / 2 # assume half of the section goes into flexural tension
    d_s = d - 2*cover - 2*h_bar_dia - 2*(0.5 * v_bar_dia) # distance from first and last layers of vert reo
    num_v_bar_layers = math.ceil(d_s / v_bar_cts) # number of layers of vert reo
    actual_spacing = d_s / num_v_bar_layers # actual spacing of bars
    first_layer = cover + h_bar_dia + 0.5*v_bar_dia
    n = 0 # number of reo layers in flexural tension zone
    Ast = [] # list of Ast in each layer i
    d_si = [] # list of d_si from extreme compression fibre
    current_layer = first_layer
    for i in range(num_v_bar_layers):
        if current_layer >= flexural_tension_zone:
            Ast.append(2 * v_bar_area)
            d_si.append(current_layer)
        current_layer = current_layer + actual_spacing

    dsy = sum(Ast[i] * d_si[i] for i in range(len(Ast))) / sum(Ast)
    dvy = max(0.72*d, 0.9*dsy)
    ast_total = round( sum(Ast[i] for i in range(len(Ast))), 0 ) # ast within tensile region

    return dvy, ast_total

def column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts):
    '''
    Determines concrete contribution to shear capacity of column based on beam shear to Section 8

    Parameters:
    fc (float): concrete strength
    cover (float): cover to reinforcement
    d (float): total depth of section
    b (float): total width of section
    v_bar_dia (float): diameter of vertical bars
    v_bar_cts (float): spacing of vertical bars
    h_bar_dia (float): diameter of horizontal bars
    h_bar_cts (float): spacing of horizontal bars

    '''
    # calculate reo area and reo rate in section
    fsy = 500 # MPa
    h_bar_area = math.pi * h_bar_dia**2 / 4 # mm2
    v_bar_area = math.pi * v_bar_dia**2 / 4 # mm2
    a_ct = d * b / 2 # mm2 - area of concrete in tensile region
    n_bars_v = __round_up_to_even(2 * ( d - 2 * cover - 2 * h_bar_dia ) / v_bar_cts) # no. vertical bars in tensile region
    n_bars_h = round(2 * ( h - 2 * cover ) / h_bar_cts, 0) # no. horizontal bars in section

    rho = ( n_bars_v * v_bar_area / (d * b) ) * 100 # reinforcement rate

    # Concrete contribution to shear cl 8.2.4
    dvy, ast_total = __calculate_effective_shear_depth(d, cover, v_bar_area, v_bar_cts, h_bar_dia)
    dvs = d / 2 - dvy
    strain_x = ( (m_star*1000000)/dvy + (v_star*1000) ) / ( 2*( dvs/ dvy * 200000 * ast_total))

    # Simple method for kv and theta v
    theta_v = 36 # angle of inclination of concrete strut
    theta_v_rads = math.radians(theta_v)
    A_sv_min = h_bar_cts * 0.08 * math.sqrt(fc) * b / fsy  # Minimum transverse shear reinforcement - Cl 8.2.1.7
    A_sv = h_bar_area * n_bars_h
    kv = 200/(100 + 1.3*d) if A_sv/h_bar_cts < A_sv_min/h_bar_cts else 0.15

    # concrete contribution to shear strength
    v_uc = round(kv * b * dvy * min(math.sqrt(fc), 8) / 1000, 0)

    # Transverse shear and torsion reinforcement contribution - Cl 8.2.5
    v_us = ((2 * h_bar_area) * fsy * dvy / h_bar_cts) * (1/math.tan(theta_v_rads)) / 1000

    # return values
    return v_uc, v_us

def moment_magnification_factor():
    pass

def column_buckling_load(d, v_bar_dia, h_bar_dia, cover, p_g, p_q):

    # Determine if column is short or slender
    # Slenderness - Cl 10.5
    
    # Radius of gyration - Cl 10.5.2
    rx = 0.3 * d

    # Column buckling load - Cl 10.4.4
    k_u = 0.545
    phi = 0.65
    m_c = 0
    d_o = d - cover - h_bar_dia - v_bar_dia/2
    alpha_s = 0.6
    E_c = 40000
    beta_d = p_g / (p_g + p_q)
    n_c = (math.pi**2 * h**2) * ((182 * d_o * m_c) / (1 + beta_d))
    
    print('r = ' + str(rx) + ' mm')

def moment_interaction_design(fc, cover, d, b, v_bar_dia, v_bar_cts, h_bar_dia, n_star, m_star):
    '''
    Checks column using moment interaction diagram to Section 10
    
    Parameters:
    fc (float): concrete strength
    cover (float): cover to reinforcement
    d (float): total depth of section
    b (float): total width of section
    v_bar_dia (float): diameter of vertical bars
    v_bar_cts (float): spacing of vertical bars
    h_bar_dia (float): diameter of horizontal bars
    h_bar_cts (float): spacing of horizontal bars

    Returns:
    pass_fail (str): returns 'Pass' if applied load and moment falls within diagram, otherwise 'Fail'
    moment_interaction_plot (Axes): plot of moment interaction diagram
    concrete_section_plot (Axes): plot showing concrete section and reo layout
    phi_n (float): intersection coordinate N with the curve
    phi_m (float): intersection coordinate M with the curve

    '''
    # Determine if column is short or slender
    # and therefore whether moment magnifier is required
    rx = 0.3 * d # Radius of gyration - Cl 10.5.2
    slenderness = 'Short' if h/rx <= 22 else 'Slender'
    magnifier_required = False if slenderness == 'Short' else True
    # Next step - add magnifier if required (for slender columns)

    # Calculate bar layout and geometry
    bar_area = math.pi * v_bar_dia**2 / 4
    n_bars_x = 2
    n_bars_y = int( round((d - 2 * cover - 2 * h_bar_dia) / v_bar_cts, 0) )

    # Calculate factors
    alpha_1 = max(0.72, min(1 - 0.003 * fc, 0.85))
    alpha_2 = max(0.85 - 0.0015*fc, 0.67)
    alpha_options = [alpha_1, alpha_2]
    gamma = 0.97-0.0025*fc

    # Set up moment interaction diagram using concreteproperties library
    design_code = AS3600() # Create AS3600 object
    concrete = design_code.create_concrete_material(compressive_strength=50) # Assign concrete section to design code
    steel = design_code.create_steel_material() # Assign steel material to design code

    # Define concrete section via sectionproperties library
    geom = concrete_rectangular_section(
        b=b, # x
        d=d, # y
        dia_top=v_bar_dia,
        area_top=bar_area,
        n_top=n_bars_x,
        c_top=cover + h_bar_dia,
        dia_bot=v_bar_dia,
        area_bot=bar_area,
        n_bot=n_bars_x,
        c_bot=cover + h_bar_dia,
        dia_side=v_bar_dia,
        area_side=bar_area,
        n_side=n_bars_y - 2,
        c_side=cover + h_bar_dia,
        conc_mat=concrete,
        steel_mat=steel,
    )
    conc_sec = ConcreteSection(geom) # Create concrete section
    design_code.assign_concrete_section(concrete_section=conc_sec) # Assign concrete section above to design code
    concrete_section_plot = conc_sec.plot_section(title='Pier')

    # Create moment interaction diagram. Returns -> tuple[MomentInteractionResults, MomentInteractionResults, list[float]]
    f_mi_res, mi_res, phis = design_code.moment_interaction_diagram(progress_bar=False, n_points=18)

    # Split up unfactored results and factored results
    results_list = mi_res.get_results_lists(moment='m_x')
    f_results_list = f_mi_res.get_results_lists(moment='m_x')
    f_n_x = [x / 1000 for x in f_results_list[0]] # Convert to kN
    f_m_x = [x / 1000000 for x in f_results_list[1]] # Convert to kNm

    # Save results in dataframe for future comparison with other software
    n_x = results_list[0]
    m_x = results_list[1]

    # Create dataframe with results
    df = pd.DataFrame({
        'N': [divmod(x, 1000)[0] for x in n_x], 
        'f_N': f_n_x, 
        'M': [divmod(x, 1000000)[0] for x in m_x],
        'f_M': f_m_x,
        'phi': phis
    })

    # Check applied axial load and moment fall within diagram
    point_in_diagram = f_mi_res.point_in_diagram(n_star, m_star) # Check if point is within diagram
    pass_fail = 'Pass' if point_in_diagram == True else 'Fail' 

    # Initialize phi_n and phi_m
    phi_n, phi_m = None, None

    # Find the intersection point
    intersection_point = __find_intersection(f_m_x, f_n_x, m_star, n_star)
    
    # Extract phi_n and phi_m from the intersection point
    if not intersection_point.is_empty:
        phi_n, phi_m = intersection_point.y, intersection_point.x

    # Plot the point (n_star, m_star) and the line from (0,0) to (n_star, m_star) using matplotlib
    plt.figure()
    plt.plot(f_m_x, f_n_x, label='Interaction Curve') # Plots interaction curve
    plt.plot(m_star, n_star, 'ro', label='Applied Load') # Plots point where N and M are applied
    plt.plot([0, m_star], [0, n_star], 'b-', label='Load Line') # Plots solid line from (0,0) to (n_star, m_star)

    # Plot the extension as a dashed line if inside the curve
    if point_in_diagram:
        plt.plot([m_star, phi_m], [n_star, phi_n], 'b--', label='Extension Line')

    plt.ylabel('Axial Load (kN)')
    plt.xlabel('Moment (kNm)')
    plt.title('Moment Interaction Diagram')
    plt.legend()
    plt.grid(True)
    plt.show()

    return pass_fail, concrete_section_plot, phi_n, phi_m


# input parameters
fc = 50
d = 2500
b = 200
h = 4200
v_bar_dia = 16
h_bar_dia = 12
cover = 30
v_bar_cts = 200
h_bar_cts = 200
n_star = 3000
m_star = 6000
v_star = 500
mu_sp = 2.6

# Check column interaction diagram
pass_fail, concrete_section_plot, phi_n, phi_m = moment_interaction_design(fc, cover, d, b, v_bar_dia, v_bar_cts, h_bar_dia, n_star, m_star)

# Check column shear capacity
v_uc, v_us = column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts)


print('Result: ' + pass_fail)
print('Phi N: ' + str(round(phi_n, 1)))
print('Phi M: ' + str(round(phi_m, 1)))
plt.show()

print('Concrete Shear Capacity: ' + str(round(v_uc, 1)) + ' kN')
print('Steel Shear Capacity: ' + str(round(v_us, 1)) + ' kN')

