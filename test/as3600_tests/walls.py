import numpy as np
import pandas as pd
import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from design_functions import load_factors
import copy

'''
Wall overall design functions

These functions will be used to design many walls at once, where each wall is represented as a dict

'''
# DESIGN ALL WALLS

def full_wall_design(walls:list[dict], load_cases:list[str], vcts:list[float], hcts:list[float], story_names:list[str], phz_levels:list[str], wall_type:str, mu_sp:float) -> list[dict]:
    piers_as_walls = []
    for wall in walls:
        wall['Design as'] = design_as(wall['Thickness Bot'], wall['Width Bot'], wall['Story Height'])
        if wall['Design as'] == 'Design as wall':
            g = load_cases[0]
            q = load_cases[1]
            rs = load_cases[2]
            # wx = load_cases[3]
            # wy = load_cases[4]

            p_g = wall[g]['p']
            p_q = wall[q]['p']
            p_rs = wall[rs]['p']
            # p_wx = wall[wx]['p']
            # p_wy = wall[wy]['p']

            v2_g = wall[g]['v2']
            v2_q = wall[q]['v2']
            v2_rs = wall[rs]['v2']

            m3_g = wall[g]['m3']
            m3_q = wall[q]['m3']
            m3_rs = wall[rs]['m3']
            # m3_wx = wall[wx]['m3']
            # m3_wy = wall[wy]['m3']

            tw = wall['Thickness Bot']
            Lw = wall['Width Bot']
            hw = wall['Story Height']
            fc = wall['fc']

            # wall prelim checks
            wall['G+0.3Q (MPa)'] = (p_g + 0.3 * p_q) / (tw * Lw) # seismic weight
            wall['G+0.3Q+RS (C)(MPa)'] = load_factors.EQCompStress(p_g, p_q, p_rs, m3_g, m3_q, m3_rs, tw, Lw)
            wall['G+0.3Q-RS (T)(MPa)'] = load_factors.EQTensStress(p_g, p_q, p_rs, m3_g, m3_q, m3_rs, tw, Lw)
            wall['Axial Load Ratio'] = axial_load_ratio(tw, Lw, wall['G+0.3Q (MPa)']) # check axial load ratio
            wall['Slenderness Ratio'] = slenderness_ratio(tw, hw, mu_sp) # check slenderness ratio

            # wall reinforcement
            wall['Rho crit.'], wall['Rho typ.'] = min_tension_reinforcement(fc, wall['Story Name'], story_names, phz_levels) # minimum reo
            wall['db Vert'], wall['s Vert'] = tension_reinforcement(tw, Lw, wall['G+0.3Q-RS (T)(MPa)'], wall['Rho crit.'], wall['Rho typ.'], vcts, 500)

            # boundary element
            wall["0.15f'c"] = 0.15 * wall['fc']
            wall["0.2f'c"] = 0.2 * wall['fc']
            wall["0.585f'c"] = 0.585 * wall['fc']
            wall['BE Width'], wall['Lig Dia'], wall['Lig Cts'] = boundary_element(tw, Lw, fc, wall['G+0.3Q+RS (C)(MPa)'], wall['db Vert'])

            # shear reinforcement
            eq_shear = load_factors.EQShear(v2_g, v2_q, v2_rs)
            wall['EQ Shear'] = adjusted_eq_shear(wall['Story Name'], mu_sp, eq_shear, phz_levels)
            wall['Vuc'], wall['Vus'], wall['phiVu'], wall['db Horiz'], wall['s Horiz'] = shear_design(tw, Lw, hw, fc, hcts, 500, wall['EQ Shear'])

            piers_as_walls.append(wall)
    
    return piers_as_walls

def optimise_walls(walls:list[dict], stress_limit:float, load_cases:list[str], mu_sp:float) -> list[dict]:
    ''' Optimise wall design based on parameters

    :param walls: original wall design dataframe
    :param stress_limit: maximum allowable compressive stress
    :param load_case: load case names to design for
    :param mu_sp: ductility factor of building

    :type walls: list[dict]
    :type stress_limit: float
    :type load_case: list[str]
    :type mu_sp: float

    :return optimised_walls: pier dataframe with optimised design
    
    '''
    optimised_walls = []
    fc_range = [32, 40, 50, 65]
    tw_range = [200, 250, 300, 350]
    g = load_cases[0]
    q = load_cases[1]
    rs = load_cases[2]
    
    for wall in walls:
        current_wall = copy.deepcopy(wall)
        wall_tw = current_wall['Thickness Bot']
        tw_index = tw_range.index(round(wall_tw,0))
        wall_fc = current_wall['fc']
        fc_index = fc_range.index(round(wall_fc, 0))
        wall_Lw = current_wall['Width Bot']
        wall_hw = current_wall['Story Height']
        p_g = current_wall[g]['p']
        p_q = current_wall[q]['p']
        p_rs = current_wall[rs]['p']
        m3_g = current_wall[g]['m3']
        m3_q = current_wall[q]['m3']
        m3_rs = current_wall[rs]['m3']
        flag = False
        for i in range(tw_index):
            if flag:
                break
            eq_load = load_factors.EQCompStress(p_g, p_q, p_rs, m3_g, m3_q, m3_rs, tw_range[i], wall_Lw) # recalculate eq compressive stress
            slenderness = slenderness_ratio(wall_tw, wall_hw, mu_sp)
            for j in range(fc_index):
                max_fc = stress_limit * fc_range[j]
                if max_fc > eq_load and slenderness != 'Slender wall. Increase thickness':
                    current_wall['thickness_bot'] = tw_range[i]
                    current_wall['fc'] = fc_range[j]
                    current_wall['Slenderness Ratio'] = slenderness
                    flag = True
                    break
        optimised_walls.append(current_wall)

    return optimised_walls

def pier_as_columns(walls:list[dict], load_cases:list[str]) -> list[dict]:
    piers_as_columns = []
    for wall in walls:
        wall['Design as'] = design_as(wall['thickness_bot'], wall['width_bot'], wall['story_height'])
        if wall['Design as'] == 'Design as column':
            g = load_cases[0]
            q = load_cases[1]
            rs = load_cases[2]
            wx = load_cases[3]
            wy = load_cases[4]

            wall['N-G'] = wall[g]['p']
            wall['N-Q'] = wall[q]['p']
            wall['N-RS'] = wall[rs]['p']

            wall['M-G'] = wall[g]['m3']
            wall['M-Q'] = wall[q]['m3']
            wall['M-RS'] = wall[rs]['m3']

    return piers_as_columns

'''
Wall design individual functions

'''

def get_phz_stories(start_phz:str, above_phz:int, story_names:list[str]) -> list[str]:
    ''' Determine levels within and outside PHZ
    
    :param start_phz: lowest level within plastic hinge region
    :param above_phz: number of levels within plastic hinge region above lowest story
    :param story_names: list of all story names in building

    :type start_phz: str
    :type above_phz: int
    :type story_names: list[str]

    :return phz_list: list of story names within plastic hinge regions
    :type phz_list: list[str]
    '''
    phz_list = []
    index = story_names.index(start_phz)
    for story in story_names:
        if story_names.index(story) < index + above_phz + 1:
            phz_list.append(story)
        else:
            break
    return phz_list

def is_in_phz(story:str, phz_stories: list) -> bool:
    ''' Determines whether a story is within the plastic hinge region

    :param story: story to be checked if within plastic hinge region
    :param phz_stories: list of stories within plastic hinge region

    :type story: str
    :type phz_stories: list[str]

    :return in_phz: returns True or False whether story is within plastic hinge region
    :type in_phz: bool
    '''

    if story in phz_stories:
        in_phz = True
    else:
        in_phz = False
    return in_phz


def design_as(tw:float, Lw:float, hw:float):
    ''' Determines whether a pier should be designed as a wall or a column.
    Check with: 
    ACI Table R18.10.1 - Governing design provisions for vertical wall segments
    
    :param tw: thickness of pier
    :param Lw: length of pier
    :param hw: floor to floor height of pier

    :type tw: float
    :type Lw: float
    :type hw: float
    
    '''
    if hw / Lw < 2:
        design_as = "Design as wall"
    elif hw / Lw >= 2 and Lw / tw <= 6:
        design_as =  "Design as column"
    else:
        design_as =  "Design as wall"
    return design_as
    

def axial_load_ratio(tw:float, Lw:float, eq_load:float) -> str:
    ''' Checks axial load ratio requirement. Only applied to ductile walls.
    AS3600:2018 - Cl 14.4.4.3 - Axial load ratio
    
    :param tw: thickness of pier
    :param Lw: length of pier
    :param eq_load: earthquake mass - G+0.3Q (kN)

    :type tw: float
    :type Lw: float
    :type eq_load: float

    '''
    Ag = tw * Lw
    if eq_load*1000/Ag < 0.2:
        return "<0.2, OK"
    else:
        return "Increase f'c or wall thickness"

#14.6.5/14/7/2 Slenderness ratio
def slenderness_ratio(tw:float, hw:float, muSp) -> str:
    if muSp == 2.6:
        if hw/tw <= 20:
            return "<20, OK"
        else:
            return "Slender wall. Increase thickness"
    elif muSp == 4.5:
        if hw/tw <= 16:
            return "<16, OK"
        else:
            return "Slender wall. Increase thickness"

# -------------------------------------------------------------------------------------------------------------
# REINFORCEMENT
# Clause 14.6.7 - Reinforcement

# MINIMUM TENSION REINFORCEMENT WITHIN PLASTIC HINGE ZONE
def min_tension_reinforcement(fc:float, story:str, story_names:list, phz_levels:list):
    fsy = 500
    rho_wv_crit = (0.7 * fc**0.5) / fsy  # minimum vertical reinforcement ratio
    rho_wv_typ = (0.35 * fc**0.5) / fsy

    if story not in phz_levels:
        current_story_index = story_names.index(story)
        end_phz_index = story_names.index(phz_levels[-1])
        difference = current_story_index - end_phz_index
        rho_wv_crit = 0 # = not required
        for i in range(difference):
            rho_wv_typ = max(rho_wv_typ * 0.9, 0.0025)
    else:
        pass
    return rho_wv_crit, rho_wv_typ

# TENSION REINFORCEMENT DUE TO MAX TENSILE LOAD
def tension_reinforcement(tw:float, Lw:float, tens_stress:float, rho_wv_crit:float, rho_wv_typ:float, vcts:list[float], fsy:float):
    if isinstance(rho_wv_crit, float):
        rho_min = rho_wv_crit
    else:
        rho_min = rho_wv_typ
    
    T_star = tens_stress * 1000 * ((tw/1000) * (Lw/1000) / 2)
    max_cts = vcts[0]
    min_cts = vcts[1]
    bar_cts = list(range(max_cts, min_cts - 1, -50))
    bar_dia = [12, 16, 20, 24, 28, 32, 36]
    for db in bar_dia:
        for s in bar_cts:
            A_t = np.pi * db**2 / 4
            n_bars = round(2 * Lw/2 / s, 0) # divided by 2 because it number of bars in peak tension zone
            rho_t = 2 * n_bars * A_t / (tw * Lw) # calculate with all bars in wall
            phi_T = 0.8 * n_bars * A_t * fsy / 1000
            if phi_T > T_star and rho_t > rho_min:
                return db, s

# -------------------------------------------------------------------------------------------------------------
# Clause 14.6.2 - Boundary elements
def boundary_element(tw:float, Lw:float, fc:float, eq_comp_stress, vdia):
    be_width = max(0.15 * Lw, 1.5 * tw)

    if vdia > 28:
        lig_dia = 12
    else:
        lig_dia = 10
    
    if eq_comp_stress >= 0.2 * fc:
        # Clause 14.5.4 - Columns
        lig_cts = min(8 * vdia, 24 * lig_dia, 0.5 * tw, 300)

    elif 0.2 * fc >= eq_comp_stress > 0.15 * fc:
        # Clause 10.7.4.3 Diameter and spacing of fitments and helices
        lig_cts = min(tw, 15 * vdia)

    elif eq_comp_stress > 0.585 * fc:
        # Clause 10.7.4.3 Diameter and spacing of fitments and helices
        be_width = "Increase f'c or wall thickness"
        lig_cts = "FAIL"
        lig_dia = "FAIL"

    else:
        be_width = "Not required"
        lig_cts = "Not required"
        lig_dia = "Not required"

    return be_width, lig_dia, lig_cts

# -------------------------------------------------------------------------------------------------------------
# Clause 11.6 - DESIGN OF WALLS FOR IN-PLANE SHEAR FORCES
def adjusted_eq_shear(story:str, muSp:float, V_star_eq:float, phz_list:list[str]):
    if story in phz_list:
        V_star = V_star_eq * muSp
    else:
        V_star = V_star_eq
    return V_star

# Clause 11.6.2 - Strength in shear
def shear_design(tw:float, Lw:float, hw:float, fc:float, hcts:list[float], fsy:float, V_star:float):
    # Clause 14.6.6 - In-plane shear
    # if within plastic hinge zone, design for full V* (assume muSp = 1.3)
    # Clause 11.6.3 - Shear strength excluding wall reinforcement
    Vuc1 = ((0.66 * np.sqrt(fc) - 0.21 * hw / Lw * np.sqrt(fc))* 0.8* Lw* tw)/1000
    if hw / Lw <= 1:
        Vuc = Vuc1
    else:
        Vuc2 = min(Vuc1,(0.05 * np.sqrt(fc) + (0.1 * np.sqrt(fc)) / (hw / Lw - 1)) * 0.8 * Lw * tw)/1000
        Vuc3 = (0.17 * np.sqrt(fc) * 0.8 * Lw * tw) / 1000
        Vuc = max(Vuc2, Vuc3)

    # Clause 11.6.4 - Contribution to shear strength by wall reinforcement
    
    def ShearSteel(db, s):
        A_b = np.pi * db**2 / 4
        n_bars = round(2 * hw/2 / s, 0) # divided by 2 because it number of bars in peak shear zone
        rho_v = 2 * n_bars * A_b / (tw * hw) # calculate with all bars in wall
        Vus = (rho_v * fsy * (0.8 * Lw * tw))/1000
        return Vus
    
    bar_dia = [12, 16, 20, 24, 28, 32, 36]
    max_cts = hcts[0]
    min_cts = hcts[1]
    bar_cts = list(range(max_cts, min_cts - 1, -50))

    for db in bar_dia:
        for s in bar_cts:
            Vus = ShearSteel(db, s)
            phiVu = 0.65 * (Vuc + Vus)

            if phiVu > V_star:
                hdia = db
                hcts = s
                return Vuc, Vus, phiVu, hdia, hcts

    


# -------------------------------------------------------------------------------------------------------------
# PRECAST DOWELS
# only include these if wall type is "Precast"
def precast_dowels(tw:float, Lw:float, story:str, V_star:float, tens_stress:float, vdia:float, vcts:float, phz_levels:float, wall_type:str):
    if wall_type == "Precast":
        fsy = 500
        T_star = tens_stress * 1000 / (tw * Lw / 2)
        # vdia to be determined once wall design is complete, dowels to be designed second
        A_v_bars = ((1000 / vcts) * 2 * (np.pi * np.power(vdia, 2) / 4))  # area of vertical bars per m
        dowel_spacing = [600, 500, 400, 300, 250, 200]  # standard spacing of dowels
        bar_dia = [12, 16, 20, 24, 28, 32, 36]
        for db in bar_dia:
            for cts in dowel_spacing:
                A_dowels = ((1000 / cts) * 2 * (np.pi * np.power(db, 2) / 4))  # assumes one row of dowels central
                phi_V = 0.65 * A_dowels * fsy / 1000
                phi_T = 0.85 * A_dowels * Lw / 2 * fsy / 1000
                if story in phz_levels and A_dowels > A_v_bars and phi_V > V_star and phi_T > T_star:
                    dowel_dia = db
                    dowel_cts = cts
                    break
                elif A_dowels > 0.5 * A_v_bars and phi_V > V_star and phi_T > T_star:
                    dowel_dia = db
                    dowel_cts = cts
                    break
                else:
                    dowel_dia = "FAIL"
                    dowel_cts = "FAIL"
                    break
    else:
        dowel_dia = "N/A"
        dowel_cts = "N/A"
    return dowel_dia, dowel_cts


# -------------------------------------------------------------------------------------------------------------
# CREATE WALL DATAFRAME
# -------------------------------------------------------------------------------------------------------------

