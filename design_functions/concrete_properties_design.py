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
bar_dia = 24
tie_dia = 12
cover = 30
n_bars_b = 2
n_bars_d = 5
n_star = 5000
m_star = 25000
v_star = 5500
mu_sp = 2.6

def moment_interaction_design(fc, d, b, bar_dia, tie_dia, cover, n_bars_b, n_bars_d, n_star, m_star):
    # calculations
    bar_area = math.pi * bar_dia**2 / 4
    alpha_1 = max(0.72, min(1 - 0.003 * fc, 0.85))
    alpha_2 = max(0.85 - 0.0015*fc, 0.67)
    alpha_options = [alpha_1, alpha_2]
    gamma = 0.97-0.0025*fc

    # moment interaction diagram
    design_code = AS3600()
    concrete = design_code.create_concrete_material(compressive_strength=50)
    steel = design_code.create_steel_material()

    geom = concrete_rectangular_section(
        b=b,
        d=d,
        dia_top=bar_dia,
        area_top=bar_area,
        n_top=n_bars_b,
        c_top=cover + tie_dia,
        dia_bot=bar_dia,
        area_bot=bar_area,
        n_bot=n_bars_b,
        c_bot=cover + tie_dia,
        dia_side=bar_dia,
        area_side=bar_area,
        n_side=n_bars_d-2,
        c_side=cover + tie_dia,
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

