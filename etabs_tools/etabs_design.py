'''
Functions to design elements directly from an ETABS model

'''

# Import libraries and modules
import sys
sys.path.append(r'C:\_Github\structural_engineering_toolbox')
from as3600 import vertical_structure
import pandas as pd
import math

def filter_bar_sizes(fc: float, story: str, story_names: list, phz_levels: list, vertical_spacing: float, b: float):
    # Define bar sizes
    bar_sizes = [12, 16, 20, 24, 28, 32, 36, 40]
    
    # Calculate minimum tension reinforcement
    rho_wv_crit, rho_wv_typ = vertical_structure.min_tension_reinforcement(fc, story, story_names, phz_levels)
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
    if p_max is not None and p_max < 0:
        p_max_tension = 0
    else:
        p_max_tension = -1 * abs(p_max)

    p_max_compression = abs(min(pier[eq_env_1 + ' Top Min']['p'], pier[eq_env_1 + ' Bottom Min']['p']) / 1000 )

    # Design pier as a column under compression load
    v_bar_size_1 = 0

    # Initialise safety factors
    safety_factor_x_c = 0
    safety_factor_x_t = 0
    safety_factor_y_c = 0
    safety_factor_y_t = 0

    pier_section = vertical_structure.ColumnSection(
            section_type=vertical_structure.SectionType.WALL, fc=fc, d=d, b=b, h=h, cover=30, v_bar_dia=12, v_bar_cts=vertical_spacing, h_bar_dia=12, h_bar_cts=horizontal_spacing)
    loading = vertical_structure.Loading(
        n_star_compression=p_max_compression, 
        n_star_tension=p_max_tension,
        m_x_top=design_m_star_top_xx, 
        m_x_bot=design_m_star_bot_xx, 
        m_y_top=design_m_star_top_yy if design_both_axes == True else 0, 
        m_y_bot=design_m_star_bot_yy if design_both_axes == True else 0, 
        v_star=design_v_star
        )

    # Check each bar size and if it works, end the for loop
    for bar_size in bar_sizes:
        # Moment interaction check with compression load
        pier_section.v_bar_dia = bar_size
        moment_interaction_results = vertical_structure.moment_interaction_design(
            section=pier_section,
            bracing='Unbraced',
            loading=loading,
            check_both_axes=design_both_axes
        )

        if design_both_axes == False and loading.n_star_tension == 0:
            if moment_interaction_results.result_x_c == 'Pass':
                v_bar_size_1 = bar_size
                safety_factor_x_c = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.m_x_c
                break
            
        elif design_both_axes == False and loading.n_star_tension != 0:
            if moment_interaction_results.result_x_c == 'Pass' and moment_interaction_results.result_x_t == 'Pass':
                v_bar_size_1 = bar_size
                safety_factor_x_c = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.phi_m_x_c
                safety_factor_x_t = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.phi_m_x_t
                break

        elif design_both_axes == True and loading.n_star_tension == 0:
            if moment_interaction_results.result_x_c == 'Pass' and moment_interaction_results.result_y_c == 'Pass':
                v_bar_size_1 = bar_size
                safety_factor_x_c = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.phi_m_x_c
                safety_factor_y_c = max(design_m_star_top_yy, design_m_star_bot_yy) / moment_interaction_results.phi_m_y_c
                break

        elif design_both_axes == True and loading.n_star_tension != 0:
            if moment_interaction_results.result_x_c == 'Pass' and moment_interaction_results.result_x_t == 'Pass' and moment_interaction_results.result_y_c == 'Pass' and moment_interaction_results.result_y_t == 'Pass':
                v_bar_size_1 = bar_size
                safety_factor_x_c = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.phi_m_x_c
                safety_factor_x_t = max(design_m_star_top_xx, design_m_star_bot_xx) / moment_interaction_results.phi_m_x_t
                safety_factor_y_c = max(design_m_star_top_yy, design_m_star_bot_yy) / moment_interaction_results.phi_m_y_c
                safety_factor_y_t = max(design_m_star_top_yy, design_m_star_bot_yy) / moment_interaction_results.phi_m_y_t
                break

            break
        
    # -------------------------------------------------------------------------------------------------

    # Determine in_plane shear capacity of pier
    vuc = 0
    vus = 0
    phi_vu = 0
    h_bar_size = 12
    for bar_size in bar_sizes:
        # Check if bar size works for in-plane shear
        pier_section.h_bar_dia = bar_size
        vuc, vus = vertical_structure.column_shear(
            section=pier_section,
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
    designed_pier['Phi Mu X-X'] = round(moment_interaction_results.phi_m_x_c if safety_factor_x_c > safety_factor_x_t else moment_interaction_results.phi_m_x_t, 0)
    designed_pier['Phi Nu X-X'] = round(moment_interaction_results.phi_n_x_c if safety_factor_x_c > safety_factor_x_t else moment_interaction_results.phi_n_x_t, 0)
    designed_pier['Safety Factor X-X'] = round(max(safety_factor_x_c, safety_factor_x_t), 2)

    if design_both_axes == True:
        designed_pier['M* Top Y-Y (kNm)'] = round(design_m_star_top_yy, 0)
        designed_pier['M* Bot Y-Y (kNm)'] = round(design_m_star_bot_yy, 0)
        designed_pier['Phi Mu Y-Y'] = round(moment_interaction_results.phi_m_y_c if safety_factor_y_c > safety_factor_y_t else moment_interaction_results.phi_m_y_t, 0)
        designed_pier['Phi Nu Y-Y'] = round(moment_interaction_results.phi_n_y_c if safety_factor_y_c > safety_factor_y_t else moment_interaction_results.phi_n_y_t, 0)
        designed_pier['Safety Factor Y-Y'] = round(max(safety_factor_y_c, safety_factor_y_t), 2)

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