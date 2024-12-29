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

# input parameters
fc = 50
d = 1000
b = 350
h = 4200
v_bar_dia = 24
h_bar_dia = 16
cover = 30
v_bar_cts = 200
h_bar_cts = 200
n_star = 5000
m_star = 25000
v_star = 5500
mu_sp = 2.6

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


def calculate_effective_shear_depth(d, cover, v_bar_area, v_bar_cts, h_bar_dia):
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
    return dvy

def column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts):
    '''
    Determines shear capacity of column based on beam shear to Section 8

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

    fsy = 500 # MPa
    h_bar_area = math.pi * h_bar_dia**2 / 4 # mm2
    v_bar_area = math.pi * v_bar_dia**2 / 4 # mm2
    n_bars_v = round(2 * ( d - 2 * cover - 2 * h_bar_dia ) / v_bar_cts, 0)
    n_bars_h = round(2 * ( h - 2 * cover ) / h_bar_cts, 0)

    # Concrete contribution to shear cl 8.2.4
    strain_x = 0
    # Simplified method for kv and theta v cl 8.2.4.3
    theta_v = (29 + 7000 * strain_x)
    A_sv_min = h_bar_cts * 0.08 * math.sqrt(fc) * b / fsy
    A_sv = h_bar_area * n_bars_h

    kv = 200/(100 + 1.3*d) if A_sv/h_bar_cts < A_sv_min/h_bar_cts else 0.15
    dvy = calculate_effective_shear_depth(d, cover, v_bar_area, v_bar_cts, h_bar_dia)
    bvy = b
    v_uc = round(kv * bvy * dvy * max(math.sqrt(fc), 8) / 1000, 0)
    v_us = round(h_bar_area * fsy * dvy / h_bar_cts * (math.cos(theta_v)/math.sin(theta_v)) / 1000, 0)
    return v_uc, v_us



v_uc, v_us = column_shear(fc, cover, d, b, h, v_bar_dia, v_bar_cts, h_bar_dia, h_bar_cts, mu_sp, v_star, m_star, n_star)
print('Vuc = ' + str(v_uc))
print('Vus = ' + str(v_us))