import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sectionproperties.pre.library import concrete_rectangular_section

from concreteproperties import (
    Concrete,
    ConcreteLinear,
    ConcreteSection,
    RectangularStressBlock,
    SteelBar,
    SteelElasticPlastic,
)

from concreteproperties.design_codes import AS3600
from concreteproperties.results import MomentInteractionResults


def __round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def __column_buckling_load(d, v_bar_dia, h_bar_dia, cover, p_g, p_q):

    # Column buckling load - Cl 10.4.4
    k_u = 0.545
    phi = 0.65
    m_c = 0
    d_o = d - cover - h_bar_dia - v_bar_dia/2
    alpha_s = 0.6
    E_c = 40000
    beta_d = p_g / (p_g + p_q)
    n_c = (math.pi**2 * h**2) * ((182 * d_o * m_c) / (1 + beta_d))
    pass

def __moment_magnification_factor():
    pass

def moment_interaction_design(fc, cover, d, b, v_bar_dia, v_bar_cts, h_bar_dia, n_star, m_star):
    # calculations
    bar_area = math.pi * v_bar_dia**2 / 4
    alpha_1 = max(0.72, min(1 - 0.003 * fc, 0.85))
    alpha_2 = max(0.85 - 0.0015*fc, 0.67)
    alpha_options = [alpha_1, alpha_2]
    gamma = 0.97-0.0025*fc
    n_bars_x = 2
    n_bars_y = round((d - 2 * cover - 2 * h_bar_dia) / v_bar_cts, 0)

    # moment interaction diagram
    design_code = AS3600()
    concrete = design_code.create_concrete_material(compressive_strength=50)
    steel = design_code.create_steel_material()

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
        n_side=n_bars_y-2,
        c_side=cover + h_bar_dia,
        conc_mat=concrete,
        steel_mat=steel,
    )
    conc_sec = ConcreteSection(geom)
    design_code.assign_concrete_section(concrete_section=conc_sec)

    f_mi_res, mi_res, phis = design_code.moment_interaction_diagram(progress_bar=False, n_points=18)

    results_list = mi_res.get_results_lists(moment='m_x')
    n_x = results_list[0]
    m_x = results_list[1]

    df = pd.DataFrame({
        'N': [divmod(x, 1000)[0] for x in n_x], 
        'M': [divmod(x, 1000000)[0] for x in m_x],
        'phi': phis
    })

    df['phiN'] = round(df['N'] * df['phi'], 0)
    df['phiM'] = round(df['M'] * df['phi'], 0)

    point_in_diagram = f_mi_res.point_in_diagram(n_star, m_star)
    pass_fail = 'Pass' if point_in_diagram == True else 'Fail'

    return pass_fail


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
    print('dsy = ' + str(dsy) + ' mm')
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

# input parameters
fc = 50
d = 3000
b = 300
h = 4200
v_bar_dia = 24
h_bar_dia = 16
cover = 30
v_bar_cts = 200
h_bar_cts = 200
n_star = 5000
m_star = 200
v_star = 500
mu_sp = 2.6

v_uc, v_us = column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts)

print('D =' + str(d) + ' mm')
print('B =' + str(b) + ' mm')
print('Vuc = ' + str(v_uc) + ' kN')
print('Vus = ' + str(v_us) + ' kN')
