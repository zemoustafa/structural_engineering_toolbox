'''
Functions to design elements directly from an ETABS model

'''

# Import libraries and modules
import sys
sys.path.append(r'C:\_Github\structural_engineering_toolbox\design_functions')
from design_functions import as3600_column_design, as3600_wall_design
import pandas as pd
import math

def filter_bar_sizes(fc: float, story: str, story_names: list, phz_levels: list, vertical_spacing: float, b: float):
    # Define bar sizes
    bar_sizes = [12, 16, 20, 24, 28, 32, 36, 40]
    
    # Calculate minimum tension reinforcement
    rho_wv_crit, rho_wv_typ = as3600_wall_design.min_tension_reinforcement(fc, story, story_names, phz_levels)
    rho_max = max(rho_wv_crit, rho_wv_typ)
    
    # Filter bar sizes based on rho calculation
    filtered_bar_sizes = []
    for bar_size in bar_sizes:
        rho = 2 * (math.pi * (bar_size ** 2) / 4) * (1000 / vertical_spacing) / (1000 * b)
        if rho > rho_max:
            filtered_bar_sizes.append(bar_size)
    
    return filtered_bar_sizes

def design_etabs_pier_as_column(
        pier:dict, 
        eq_env_1:str, 
        eq_env_2:str, 
        wind_env:str, 
        vertical_spacing:float, 
        horizontal_spacing:float, 
        design_both_axes:bool
        ):
    '''
    Designs an ETABS pier as a column to AS3600

    Parameters:
    pier (dict): dict representing a single pier in ETABS
    eq_env_1 (str): name of unfactored ULS earthquake load combination envelope  
    eq_env_2 (str): name of factored ULS earthquake load combination envelope  
    wind_env (str): name of ULS wind load combination envelope  
    vertical_spacing (float): user defined spacing of vertical reinforcement
    horizontal_spacing (float): user defined spacing of horizontal reinforcement
    design_both_axes (bool): if True, design weak and strong axis. if False, design strong only
    
    Returns:
    designed_pier (dict): updated dict of designed pier

    '''
    
    # Extract pier dimensionsand story
    d = pier['Width Bot'] # Length of pier
    b = pier['Thickness Bot'] # Thickness of pier
    h = pier['Story Height'] # Height of pier
    fc = pier['fc'] # Concrete strength

    # Define bar sizes
    bar_sizes = [12, 16, 20, 24, 28, 32, 36, 40]

    # Determine design shear force
    design_v_star = max(
        abs(pier[eq_env_2 + ' Top Max']['v2']),
        abs(pier[eq_env_2 + ' Bottom Max']['v2']),
        abs(pier[eq_env_2 + ' Top Min']['v2']),
        abs(pier[eq_env_2 + ' Bottom Min']['v2']),
        abs(pier[wind_env + ' Top Max']['v2']),
        abs(pier[wind_env + ' Bottom Max']['v2']),
        abs(pier[wind_env + ' Top Min']['v2']),
        abs(pier[wind_env + ' Bottom Min']['v2'])
    ) / 1000

    # Determine design moments in strong and weak directions and corresponding moments at the other end
    design_moment_x_x = max(
        abs(pier[eq_env_1 + ' Top Max']['m3']),
        abs(pier[eq_env_1 + ' Bottom Max']['m3']),
        abs(pier[eq_env_1 + ' Top Min']['m3']),
        abs(pier[eq_env_1 + ' Bottom Min']['m3']),
        abs(pier[wind_env + ' Top Max']['m3']),
        abs(pier[wind_env + ' Bottom Max']['m3']),
        abs(pier[wind_env + ' Top Min']['m3']),
        abs(pier[wind_env + ' Bottom Min']['m3'])
    )

    def get_design_moments(moment, env, top_max, bot_max, top_min, bot_min, moment_key):
        if moment == abs(pier[env + top_max][moment_key]):
            return abs(pier[env + top_max][moment_key]), abs(pier[env + bot_max][moment_key])
        elif moment == abs(pier[env + bot_max][moment_key]):
            return abs(pier[env + bot_max][moment_key]), abs(pier[env + top_max][moment_key])
        elif moment == abs(pier[env + top_min][moment_key]):
            return abs(pier[env + top_min][moment_key]), abs(pier[env + bot_min][moment_key])
        elif moment == abs(pier[env + bot_min][moment_key]):
            return abs(pier[env + bot_min][moment_key]), abs(pier[env + top_min][moment_key])

    design_m_star_top_xx, design_m_star_bot_xx = get_design_moments(design_moment_x_x, eq_env_1, ' Top Max', ' Bottom Max', ' Top Min', ' Bottom Min', 'm3')
    design_m_star_top_yy = 0
    design_m_star_bot_yy = 0

    design_m_star_top_xx = design_m_star_top_xx / 1e6
    design_m_star_bot_xx = design_m_star_bot_xx / 1e6 

    # Design axial load
    p_max = max(pier[eq_env_1 + ' Top Max']['p'], pier[eq_env_1 + ' Bottom Max']['p']) / 1000
    p_max_tension = 0 if p_max < 0 else p_max
    p_max_compression = abs (min(pier[eq_env_1 + ' Top Min']['p'], pier[eq_env_1 + ' Bottom Min']['p']) / 1000 )

    # Design pier as a column under compression load
    v_bar_size_1 = 0

    # Initialise results of capacity
    phi_mu_comp_x_x = 0
    phi_nu_comp_x_x = 0
    phi_mu_tens_x_x = 0
    phi_nu_tens_x_x = 0
    
    phi_mu_comp_y_y = 0
    phi_nu_comp_y_y = 0
    phi_mu_tens_y_y = 0
    phi_nu_tens_y_y = 0

    # Check each bar size and if it works, end the for loop
    for bar_size in bar_sizes:
        # Moment interaction check with compression load
        moment_interaction_results_p_min = as3600_column_design.moment_interaction_design(
            fc=fc,
            cover=30,
            d=d,
            b=b,
            h=h,
            v_bar_dia=bar_size,
            v_bar_cts=vertical_spacing,
            h_bar_dia=12,
            bracing='Unbraced',
            n_star=p_max_compression,
            m_star_xx_top=design_m_star_top_xx,
            m_star_xx_bot=design_m_star_bot_xx,
            m_star_yy_top=design_m_star_top_yy,
            m_star_yy_bot=design_m_star_bot_yy,
            check_both_axes=design_both_axes
        )

        if design_both_axes == True:
            if moment_interaction_results_p_min[0][0] == 'Pass' and moment_interaction_results_p_min[0][1] == 'Pass':
                v_bar_size_1 = bar_size
                phi_mu_comp_x_x = moment_interaction_results_p_min[1][0]
                phi_nu_comp_x_x = moment_interaction_results_p_min[1][1]
                
                phi_mu_comp_y_y = moment_interaction_results_p_min[2][0]
                phi_nu_comp_y_y = moment_interaction_results_p_min[2][1]
                break
            break
        else:
            if moment_interaction_results_p_min[0] == 'Pass':
                v_bar_size_1 = bar_size
                phi_mu_comp_x_x = moment_interaction_results_p_min[1][0]
                phi_nu_comp_x_x = moment_interaction_results_p_min[1][1]
                break
            break


    # Determine pertentage capacity
    safety_factor_x_x = max(design_m_star_top_xx, design_m_star_bot_xx) / phi_mu_comp_x_x

    # WEAK AXIS
    # -------------------------------------------------------------------------------------------------
    # Repeat calculations if check both axes is true
    if design_both_axes == True:
        design_moment_y_y = max(
            abs(pier[eq_env_1 + ' Top Max']['m2']),
            abs(pier[eq_env_1 + ' Bottom Max']['m2']),
            abs(pier[eq_env_1 + ' Top Min']['m2']),
            abs(pier[eq_env_1 + ' Bottom Min']['m2']),
            abs(pier[wind_env + ' Top Max']['m2']),
            abs(pier[wind_env + ' Bottom Max']['m2']),
            abs(pier[wind_env + ' Top Min']['m2']),
            abs(pier[wind_env + ' Bottom Min']['m2'])
        )

        design_m_star_top_yy, design_m_star_bot_yy = get_design_moments(design_moment_y_y, eq_env_1, ' Top Max', ' Bottom Max', ' Top Min', ' Bottom Min', 'm2')
        design_m_star_top_yy = design_m_star_top_yy / 1e6
        design_m_star_bot_yy = design_m_star_bot_yy / 1e6
        safety_factor_y_y = max(design_m_star_top_yy, design_m_star_bot_yy) / phi_mu_comp_y_y        

    # -------------------------------------------------------------------------------------------------

    # Determine in_plane shear capacity of pier
    vuc = 0
    vus = 0
    phi_vu = 0
    for bar_size in bar_sizes:
        # Check if bar size works for in-plane shear
        vuc, vus = as3600_column_design.column_shear(
            fc=fc,
            cover=30,
            d=d,
            b=b,
            h=h,
            v_bar_dia=v_bar_size_1,
            v_bar_cts=vertical_spacing,
            h_bar_dia=bar_size,
            h_bar_cts=horizontal_spacing,
        )

        phi_vu = 0.65 * (vuc + vus)
        if phi_vu >= design_v_star:
            h_bar_size = bar_size
            break

    safety_factor_shear =  design_v_star / phi_vu

    # Add design data to designed pier dictionary
    designed_pier = {}
    designed_pier['Pier'] = pier['Pier Name']
    designed_pier['Story']= pier['Story Name']
    designed_pier['Lw (mm)'] = round(d, 0)
    designed_pier['tw (mm)'] = round(b, 0)
    designed_pier['Hw (mm)'] = round(h, 0)
    designed_pier['fc (MPa)'] = round(fc, 0)
    designed_pier['V* (kN)'] = round(design_v_star, 0)
    designed_pier['Phi Vu'] = round(phi_vu, 0)
    designed_pier['H Bar Size'] = h_bar_size
    designed_pier['Safety Factor Shear'] = round(safety_factor_shear, 2)
    designed_pier['P* Tension (kN)'] = round(p_max_tension, 0)
    designed_pier['P* Compression (kN)'] = round(p_max_compression, 0)
    designed_pier['V Bar Size'] = v_bar_size_1
    designed_pier['M* Top X-X (kNm)'] = round(design_m_star_top_xx, 0)
    designed_pier['M* Bot X-X (kNm)'] = round(design_m_star_bot_xx, 0)
    designed_pier['Phi Mu X-X'] = round(phi_mu_comp_x_x, 0)
    designed_pier['Phi Nu X-X'] = round(phi_nu_comp_x_x, 0)
    designed_pier['Safety Factor X-X'] = round(safety_factor_x_x, 2)

    if design_both_axes == True:
        designed_pier['M* Top Y-Y (kNm)'] = round(design_m_star_top_yy, 0)
        designed_pier['M* Bot Y-Y (kNm)'] = round(design_m_star_bot_yy, 0)
        designed_pier['Phi Mu Y-Y'] = round(phi_mu_comp_y_y, 0)
        designed_pier['Phi Nu Y-Y'] = round(phi_nu_comp_y_y, 0)
        designed_pier['Safety Factor Y-Y'] = round(safety_factor_y_y, 2)

    return designed_pier

def design_all_piers(
        piers:list[dict], 
        eq_env_1:str, 
        eq_env_2:str, 
        wind_env:str, 
        vertical_spacing:float, 
        horizontal_spacing:float,
        design_both_axes:bool
        ):
    '''
    Design all piers in an ETABS model in one go

    Parameters:
    pier (dict): dict representing a single pier in ETABS
    eq_env_1 (str): name of unfactored ULS earthquake load combination envelope  
    eq_env_2 (str): name of factored ULS earthquake load combination envelope  
    wind_env (str): name of ULS wind load combination envelope  
    vertical_spacing (float): user defined spacing of vertical reinforcement
    horizontal_spacing (float): user defined spacing of horizontal reinforcement

    Returns:
    designed_pier (dict): updated dict of designed pier

    '''

    # Initialise list to store designed piers
    designed_piers = []

    # Iterate through all piers and design them
    for pier in piers:
        designed_pier = design_etabs_pier_as_column(
            pier, 
            eq_env_1, 
            eq_env_2, 
            wind_env, 
            vertical_spacing, 
            horizontal_spacing,
            design_both_axes,
            )
        designed_piers.append(designed_pier)

    # Create a dataframe to store the results
    df = pd.DataFrame(designed_piers)
    return df