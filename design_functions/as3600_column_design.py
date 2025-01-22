import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sectionproperties.pre.library import concrete_rectangular_section

from concreteproperties import Concrete, ConcreteLinear, ConcreteSection, RectangularStressBlock, SteelBar, SteelElasticPlastic

from concreteproperties.design_codes import AS3600
from concreteproperties.results import MomentInteractionResults
from shapely.geometry import LineString, Point, Polygon

def __round_up_to_even(f):
    return math.ceil(f / 2.) * 2


def __find_intersection(f_m, f_n, m_star, n_star):
    """
    Find the intersection of a line with a curve.

    Args:
        f_m_x (list of floats): x-coordinates of the curve.
        f_n_x (list of floats): y-coordinates of the curve.
        m_star (float): x-coordinate of the end point of the line.
        n_star (float): y-coordinate of the end point of the line.

    Returns:
        Point: Intersection point between the extended line and the curve.
    """
    # Combine x and y coordinates into a list of tuples for the curve
    curve_coords = list(zip(f_m, f_n))
    
    # Define the curve as a LineString
    curve = LineString(curve_coords)
    
    # Define the line extending from (0, 0) through (m_star, n_star)
    direction = (m_star, n_star)
    extended_line = LineString([(0, 0), (10 * direction[0], 10 * direction[1])])  # Extend by a factor of 10 or more

    # Find the intersection
    intersection = curve.intersection(extended_line)

    # Extract phi_n and phi_m from the intersection point
    if intersection.geom_type == 'Point':
        phi_n, phi_m = intersection.y, intersection.x
        return phi_n, phi_m
    elif intersection.geom_type == 'MultiPoint':
        point_1 = intersection.geoms[0]
        point_2 = intersection.geoms[1]

        intersection = point_1 if point_1.coords[0][0] > point_2.coords[0][0] else point_2

    phi_n, phi_m = intersection.y, intersection.x

    return phi_n, phi_m


def __is_point_inside_curve(f_m_x, f_n_x, m_star, n_star):
    """
    Determines whether the point (m_star, n_star) lies inside the curve.

    Args:
        f_m_x (list of floats): x-coordinates of the curve.
        f_n_x (list of floats): y-coordinates of the curve.
        m_star (float): x-coordinate of the point to test.
        n_star (float): y-coordinate of the point to test.

    Returns:
        bool: True if the point is inside the curve, False otherwise.
    """
    # Combine coordinates into a list of tuples and ensure closure
    f_m_x.append(f_m_x[0]) # Close the curve
    f_n_x.append(f_n_x[0]) # Close the curve

    curve_coords = list(zip(f_m_x, f_n_x)) # Combine x and y coordinates into a list of tuples

    polygon = Polygon(curve_coords) # Create a polygon from the curve coordinates

    point = Point(n_star, m_star) # Create a point from the input coordinates

    # Check if the line is entirely within the polygon
    is_in_curve = polygon.contains(point)

    return is_in_curve


def __extract_balance_point(mi_res):
    """
    Extracts the balance point from the moment interaction results.

    Args:
        mi_res (MomentInteractionResults): Moment interaction results object

    Returns:
        tuple: Balance point as (moment, axial load) (Nmm, N)
    """
    for result in mi_res.results:
        if round(result.k_u, 3) == 0.545:
            return (result.m_x, result.m_y, result.n)  # Unfactored moment and axial load at balance point
    return None


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
    # strain_x = ( (m_star*1000000)/dvy + (v_star*1000) ) / ( 2*( dvs/ dvy * 200000 * ast_total))

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


def moment_magnification_factor(fc, bracing, d, b, n_c, n_star, m_star_top, m_star_bot, l_e, r, beta_d=0.5):
    '''
    Determines moment amplification factor for slender columns based on Cl 10.4

    Parameters:
    bracing (str): bracing condition of column
    d (float): total depth of section
    b (float): total width of section
    n_c (float): buckling load of column
    n_star (float): applied axial load
    m_star_top (float): applied moment at top of column
    m_star_bot (float): applied moment at bottom of column
    l_e (float): effective length of column
    r (float): radius of gyration of column

    Returns:
    magnification_factor (float): moment magnification factor

    '''
    fc_options = [20, 25, 32, 40, 50, 65, 80, 100] # MPa
    e_c_options = [24000, 26700, 30100, 32800, 34800, 37400, 39600, 42200, 44400] # MPa

    # Check if provided fc is in the options
    if fc not in fc_options:
        raise ValueError(f"Invalid concrete strength (fc). Must be one of {fc_options}")

    if bracing == 'Braced':
        m1 = min(m_star_top, m_star_bot)
        m2 = max(m_star_top, m_star_bot)

        m1_m2_ratio = -1 if m1 / m2 <= 0.05 * d * n_star else m1 / m2 
        km = 0.6 - 0.4 * m1_m2_ratio

        magnification_factor = max(1, km / (1 - n_star / n_c) ) # Cl 10.4.2

    elif bracing == 'Unbraced':
        beta_d = 0 if l_e / r <= 40 and n_star <= m_star_bot / 2 * d and n_star <= m_star_top / 2 * d else beta_d

        i_f = b * d**3 / 12 # moment of inertia of section
    
        # Determine e_c by matching the index of fc in fc_options
        e_c = e_c_options[fc_options.index(fc)] # Modulus of elasticity for concrete
        lambda_uc = 0.8 * e_c * i_f # Ratio of elastic critical buckling load to design load
        magnification_factor = 1 / (1 - (1 + beta_d) / (0.6 * lambda_uc)) # Cl 10.4.3 (b)
    else:
        raise ValueError("Invalid bracing condition. Must be 'Braced' or 'Unbraced'.")
    
    return magnification_factor


def column_buckling_load(d, v_bar_dia, h_bar_dia, cover, m_ub, beta_d):
    '''
    Determine buckling load of column based on Section 10

    Parameters:
    d (float): total depth of section
    v_bar_dia (float): diameter of vertical bars
    h_bar_dia (float): diameter of horizontal bars
    cover (float): cover to reinforcement
    m_ub (float): moment at balance point
    beta_d (float): ratio of G/(G+Q) for column

    Returns:
    n_c (float): buckling load of column (N)

    '''

    # Calculate effective depth of section
    d_0 = d - cover - h_bar_dia - v_bar_dia/2

    # Calculate buckling load - Cl 10.4.4
    n_c = round( (math.pi**2 / h**2) * ((182 * d_0) * 0.65 * m_ub / (1 + beta_d)) / 1000, 1 ) # kN
    return n_c


def moment_interaction_design(fc, cover, d, b, v_bar_dia, v_bar_cts, h_bar_dia, bracing, n_star, m_star_xx_top, m_star_xx_bot, m_star_yy_top, m_star_yy_bot):
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
    bracing (str): bracing condition of column. Value must be 'Braced' or 'Unbraced'
    n_star (float): applied axial load
    m_star_xx_top (float): applied moment at top of column about x axis
    m_star_xx_bot (float): applied moment at bottom of column about x axis
    m_star_yy_top (float): applied moment at top of column about y axis
    m_star_yy_bot (float): applied moment at bottom of column about y axis

    Returns:
    results (tuple): tuple containing the following:
        pass_fail (tuple): tuple containing pass/fail results for x and y axes
        capacity_x (tuple): tuple containing phi_m and phi_n for x axis
        capacity_y (tuple): tuple containing phi_m and phi_n for y axis
        fig (matplotlib.figure.Figure): matplotlib figure object

    '''
    def create_concrete_section(b, d, n_bars_x, n_bars_y):
        bar_area = math.pi * v_bar_dia**2 / 4
        geom = concrete_rectangular_section(
            b=b, d=d,
            dia_top=v_bar_dia, area_top=bar_area, n_top=n_bars_x, c_top=cover + h_bar_dia,
            dia_bot=v_bar_dia, area_bot=bar_area, n_bot=n_bars_x, c_bot=cover + h_bar_dia,
            dia_side=v_bar_dia, area_side=bar_area, n_side=n_bars_y - 2, c_side=cover + h_bar_dia,
            conc_mat=concrete, steel_mat=steel,
        )
        return ConcreteSection(geom)

    def get_moment_interaction_results(conc_sec):
        design_code.assign_concrete_section(concrete_section=conc_sec)
        return design_code.moment_interaction_diagram(progress_bar=False, n_points=18, control_points=[("fy", 1.0)])

    def process_results(f_mi_res):
        f_results_list = f_mi_res.get_results_lists(moment='m_x')
        f_n = [x / 1000 for x in f_results_list[0]]  # Convert to kN
        f_m = [x / 1000000 for x in f_results_list[1]]  # Convert to kNm
        return f_n, f_m

    def calculate_magnified_moment(slenderness, mi_res, m_star, r, axis):
        if slenderness == 'Slender':
            balance_point = __extract_balance_point(mi_res)
            n_c = column_buckling_load(d if axis == 'x' else b, v_bar_dia, h_bar_dia, cover, balance_point[0], beta_d=0.5)
            delta = moment_magnification_factor(fc, bracing, d if axis == 'x' else b, b if axis == 'x' else d, n_c, n_star, m_star_xx_bot if axis == 'x' else m_star_yy_bot, m_star_xx_bot if axis == 'x' else m_star_yy_bot, h, r, beta_d=0.5)
            return m_star * delta
        return m_star

    # Determine if column is short or slender - Cl 10.2
    rx, ry = 0.3 * d, 0.3 * b  # Radius of gyration - Cl 10.5.2
    slenderness_x = 'Short' if h / rx <= 22 else 'Slender'
    slenderness_y = 'Short' if h / ry <= 22 else 'Slender'

    # Calculate bar layout and geometry
    n_bars_x, n_bars_y = 2, int(round((d - 2 * cover - 2 * h_bar_dia) / v_bar_cts, 0))

    # Set up moment interaction diagram using concreteproperties library
    design_code = AS3600()
    concrete = design_code.create_concrete_material(compressive_strength=50)
    steel = design_code.create_steel_material()

    # Create and process moment interaction results for both axes
    conc_sec_x = create_concrete_section(b, d, n_bars_x, n_bars_y)
    f_mi_res_xx, mi_res_xx, _ = get_moment_interaction_results(conc_sec_x)
    f_n_x, f_m_x = process_results(f_mi_res_xx)

    conc_sec_y = create_concrete_section(d, b, n_bars_y, n_bars_x)
    f_mi_res_yy, mi_res_yy, _ = get_moment_interaction_results(conc_sec_y)
    f_n_y, f_m_y = process_results(f_mi_res_yy)

    # Determine max applied moment in each direction
    m_star_x, m_star_y = max(abs(m_star_xx_top), abs(m_star_xx_bot)), max(abs(m_star_yy_top), abs(m_star_yy_bot))

    # Apply moment magnification if column is slender
    m_star_x = calculate_magnified_moment(slenderness_x, mi_res_xx, m_star_x, rx, 'x')
    m_star_y = calculate_magnified_moment(slenderness_y, mi_res_yy, m_star_y, ry, 'y')


    # Check applied axial load and moment fall within diagram
    point_in_diagram_x = __is_point_inside_curve(f_m_x, f_n_x, n_star, m_star_x)
    point_in_diagram_y = __is_point_inside_curve(f_m_y, f_n_y, n_star, m_star_y)
    pass_fail = ('Pass' if point_in_diagram_x else 'Fail', 'Pass' if point_in_diagram_y else 'Fail')

    # Find the intersection points
    phi_n_x, phi_m_x = __find_intersection(f_m_x, f_n_x, m_star_x, n_star)
    capacity_x = (phi_m_x, phi_n_x)
    phi_n_y, phi_m_y = __find_intersection(f_m_y, f_n_y, m_star_y, n_star)
    capacity_y = (phi_m_y, phi_n_y)

    # Plot the interaction curves and applied loads
    fig, ax = plt.subplots()
    ax.plot(f_m_x, f_n_x, label='Interaction Curve - About XX')
    ax.plot(m_star_x, n_star, 'ro', label='Applied Load')
    ax.plot([0, m_star_x], [0, n_star], 'b-', label='Load Line')
    ax.plot(f_m_y, f_n_y, label='Interaction Curve - About YY')
    ax.plot(m_star_y, n_star, 'ro', label='Applied Load')
    ax.plot([0, m_star_y], [0, n_star], 'b-', label='Load Line')
    if point_in_diagram_x:
        ax.plot([m_star_x, phi_m_x], [n_star, phi_n_x], 'b--', label='Extension Line')
    if point_in_diagram_y:
        ax.plot([m_star_y, phi_m_y], [n_star, phi_n_y], 'b--', label='Extension Line')
    ax.set_ylabel('Axial Load (kN)')
    ax.set_xlabel('Moment (kNm)')
    ax.set_title('Moment Interaction Diagram')
    ax.legend()
    ax.grid(True)

    results = (pass_fail, capacity_x, capacity_y, fig)
    return results


# input parameters
fc = 50
d = 2000
b = 250
h = 4200
v_bar_dia = 28
h_bar_dia = 20
cover = 30
v_bar_cts = 200
h_bar_cts = 200
n_star = 3000
m3_top = 3000
m3_bot = -1000
m2_top = -100
m2_bot = 100
v_star = 500
mu_sp = 2.6

# Check column interaction diagram
moment_interaction_results = moment_interaction_design(
    fc, cover, d, b, v_bar_dia, v_bar_cts, h_bar_dia, 'Unbraced', n_star, m3_top, m3_bot, m2_top, m2_bot
)

# Check column shear capacity
v_uc, v_us = column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts)

print('Column Design Results')
print('---------------------')
print('Column Dimensions: ' + str(b) + ' x ' + str(d) + ' mm')
print('Concrete Strength: ' + str(fc) + ' MPa')
print('Vertical Bar Diameter: ' + str(v_bar_dia) + ' mm')
print('Vertical Bar Spacing: ' + str(v_bar_cts) + ' mm')
print('Horizontal Bar Diameter: ' + str(h_bar_dia) + ' mm')
print('Horizontal Bar Spacing: ' + str(h_bar_cts) + ' mm')
print('Cover to Reinforcement: ' + str(cover) + ' mm')
print('Axial Load: ' + str(n_star) + ' kN')

print('Results X')
print('Result X: ' + moment_interaction_results[0][0]) # pass or fail
print('Phi N X: ' + str(round(moment_interaction_results[1][1], 1))) # intersection with N
print('Phi M X: ' + str(round(moment_interaction_results[1][0], 1))) # intersection with M

print('Results X')
print('Result Y: ' + moment_interaction_results[0][1]) # pass or fail
print('Phi N Y: ' + str(round(moment_interaction_results[2][1], 1))) # intersection with N
print('Phi M Y: ' + str(round(moment_interaction_results[2][0], 1))) # intersection with M

shear_pass_fail = 'Pass' if 0.7*v_uc + 0.7*v_us > v_star else 'Fail'
print('Shear Capacity Results: ' + shear_pass_fail)
print('Concrete Shear Capacity: ' + str(round(v_uc, 1)) + ' kN') # concrete shear capacity
print('Steel Shear Capacity: ' + str(round(v_us, 1)) + ' kN') # steel shear capacity

plt.show()
