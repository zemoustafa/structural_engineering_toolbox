'''
Functions to design elements directly from an ETABS model

'''

# Import libraries and modules
import sys
sys.path.append(r'C:\_Github\structural_engineering_toolbox')
from as3600 import vertical_structure_v1
import pandas as pd
import math

class PierColumnDesigner:
    def __init__(self, pier, eq_env_1, eq_env_2, wind_env, vertical_spacing, horizontal_spacing, design_both_axes, design_boundary_element):
        self.pier = pier
        self.eq_env_1 = eq_env_1
        self.eq_env_2 = eq_env_2
        self.wind_env = wind_env
        self.vertical_spacing = vertical_spacing
        self.horizontal_spacing = horizontal_spacing
        self.design_both_axes = design_both_axes
        self.design_boundary_element = design_boundary_element
        self.bar_sizes = [12, 16, 20, 24, 28, 32, 36, 40]
        self.extract_pier_properties()
    
    def extract_pier_properties(self):
        self.d = self.pier['Width Bot']
        self.b = self.pier['Thickness Bot']
        self.h = self.pier['Story Height']
        self.fc = self.pier['fc']

    def get_design_shear_force(self):
        return max(
            abs(self.pier[self.eq_env_2 + ' Top Max']['v2']),
            abs(self.pier[self.eq_env_2 + ' Bottom Max']['v2']),
            abs(self.pier[self.eq_env_2 + ' Top Min']['v2']),
            abs(self.pier[self.eq_env_2 + ' Bottom Min']['v2']),
            abs(self.pier[self.wind_env + ' Top Max']['v2']),
            abs(self.pier[self.wind_env + ' Bottom Max']['v2']),
            abs(self.pier[self.wind_env + ' Top Min']['v2']),
            abs(self.pier[self.wind_env + ' Bottom Min']['v2'])
        ) / 1000
    
    def get_design_moment(self, env, moment_key):
        '''
        Sorts out the design moment for the pier

        Parameters:
        env (str): name of earthquake combination envelope. Eg. '(88) RS ULS ENV'. Must be envelope used for flexure design
        moment_key (str): key to access the moment value in the pier dictionary. Eg. 'm3' for strong axis, 'm2' for weak axis

        Returns:
        design_m_star_top (float): design moment for top of pier
        design_m_star_bot (float): design moment for bottom of pier

        '''
        # Take max of earthquake and wind envelopes
        design_m_star_top = max(
            abs(self.pier[env + ' Top Max'][moment_key]),
            abs(self.pier[env + ' Top Min'][moment_key]),
            abs(self.pier[self.wind_env + ' Top Max'][moment_key]),
            abs(self.pier[self.wind_env + ' Top Min'][moment_key])
        ) / 1e6
        design_m_star_bot = max(
            abs(self.pier[env + ' Bottom Max'][moment_key]),
            abs(self.pier[env + ' Bottom Min'][moment_key]),
            abs(self.pier[self.wind_env + ' Bottom Max'][moment_key]),
            abs(self.pier[self.wind_env + ' Bottom Min'][moment_key])
        ) / 1e6
        return design_m_star_top, design_m_star_bot
    
    def get_design_axial_loads(self):
        p_min = -1 * max(self.pier[self.eq_env_1 + ' Top Max']['p'], self.pier[self.eq_env_1 + ' Bottom Max']['p']) / 1000
        p_max = abs(min(self.pier[self.eq_env_1 + ' Top Min']['p'], self.pier[self.eq_env_1 + ' Bottom Min']['p']) / 1000)
        return p_max, p_min
    
    def design_pier(self):
        design_v_star = self.get_design_shear_force()
        design_m_star_top_x, design_m_star_bot_x = self.get_design_moment(self.eq_env_1, 'm3')
        design_m_star_top_y, design_m_star_bot_y = (self.get_design_moment(self.eq_env_1, 'm2') if self.design_both_axes else (0, 0))

        p_max, p_min = self.get_design_axial_loads()

        pier_section = vertical_structure_v1.RectangularColumn(
            section_type=vertical_structure_v1.SectionType.WALL, 
            fc=self.fc, 
            d=self.d, 
            b=self.b, 
            h=self.h, 
            cover=30, 
            v_bar_dia=12, 
            v_bar_cts=self.vertical_spacing, 
            h_bar_dia=12,
            h_bar_cts=self.horizontal_spacing,
            bracing_x=vertical_structure_v1.BracingType.UNBRACED,
            bracing_y=vertical_structure_v1.BracingType.BRACED
        )

        loading = vertical_structure_v1.Loading(
            n_star_max=p_max, 
            n_star_min=p_min,
            m_x_top=design_m_star_top_x, 
            m_x_bot=design_m_star_bot_x, 
            m_y_top=design_m_star_top_y if self.design_both_axes else 0, 
            m_y_bot=design_m_star_bot_y if self.design_both_axes else 0, 
            v_star=design_v_star
        )

        # Design for flexure
        # Iterate through each bar size to find the first one that passes
        for bar_size in self.bar_sizes:
            pier_section.v_bar_dia = bar_size
            mi_results = vertical_structure_v1.rectangular_column_moment_interaction(
                section=pier_section,
                design_both_axes=self.design_both_axes
            )

            column_check_results = vertical_structure_v1.check_column_capacity(
                section=pier_section,
                loading=loading,
                mi_results=mi_results,
                )

            # Check if passes or fails
            # First, check if design_both_axes is True
            if self.design_both_axes: # If design_both_axes is True
                if column_check_results.result_x == 'Pass' and column_check_results.result_y == 'Pass' and column_check_results.buckling_x == False and column_check_results.buckling_y == False:
                    # Check utilisation of column capacity for both axes
                    utilisation_x_1 = max(design_m_star_top_x, design_m_star_bot_x) / column_check_results.phi_m_x_max # Check utilisation under compression load
                    utilisation_x_2 = max(design_m_star_top_x, design_m_star_bot_x) / column_check_results.phi_m_x_min if column_check_results.phi_m_x_min is not None else 0 # Check utilisation under tension load (if applicable)
                    utilisation_y_1 = max(design_m_star_top_y, design_m_star_bot_y) / column_check_results.phi_m_y_max # Check utilisation under compression load
                    utilisation_y_2 = max(design_m_star_top_y, design_m_star_bot_y) / column_check_results.phi_m_y_min if column_check_results.phi_m_x_min is not None else 0 # Check utilisation under tension load (if applicable)
                    break
            else: # If design_both_axes is False
                if column_check_results.result_x == 'Pass' and column_check_results.buckling_x == False:
                    # Check utilisation of column capacity for strong axis only
                    utilisation_x_1 = max(design_m_star_top_x, design_m_star_bot_x) / column_check_results.phi_m_x_max # Check utilisation under compression load
                    utilisation_x_2 = max(design_m_star_top_x, design_m_star_bot_x) / column_check_results.phi_m_x_min if column_check_results.phi_m_x_min is not None else 0 # Check utilisation under tension load (if applicable)
                    break

        # Design for shear
        # Iterate through each bar size to find the first one that passes
        for bar_size in self.bar_sizes:
            # Update bar size for current iteration
            pier_section.h_bar_dia = bar_size
            
            # Calculate shear capacity
            v_uc, v_us = vertical_structure_v1.column_shear_capacity(pier_section)
            
            # Calculate phiVu
            phi_vu = 0.65 * v_us

            # Check if passes or fails
            safety_factor_shear = loading.v_star / phi_vu
            if safety_factor_shear <= 1:
                break

        # Store designed pier results into a dictionary
        designed_pier = {
            'Pier': self.pier['Pier Name'],
            'Story': self.pier['Story Name'],
            'Lw (mm)': round(self.d, 0),
            'tw (mm)': round(self.b, 0),
            'Hw (mm)': round(self.h, 0),
            'fc (MPa)': round(self.fc, 0),
            'V* (kN)': round(design_v_star, 0),
            'phiVu (kN)': round(phi_vu, 0),
            'H bar size': pier_section.h_bar_dia,
            'Utilisation Shear': round(safety_factor_shear, 2),
            'P* Max (kN)': round(p_max, 0),
            'P* Min (kN)': round(p_min, 0),
            'M* X (kNm)': round(max(design_m_star_top_x, design_m_star_bot_x), 0),
            'phiMu X Max (kNm)': round(column_check_results.phi_m_x_max, 0),
            'phiMu X Min (kNm)': round(column_check_results.phi_m_x_min, 0) if column_check_results.phi_m_x_min is not None else 0,
            'V bar size': pier_section.v_bar_dia,
            'Utilisation X Max':round(utilisation_x_1, 2),
            'Utilisation X Min':round(utilisation_x_2, 2)
        }

        # If designing both axes, add additional results
        if self.design_both_axes:
            designed_pier['M* Y (kNm)'] = round(max( design_m_star_top_y, design_m_star_bot_y), 0)
            designed_pier['phiMu Y Max (kNm)'] = round(column_check_results.phi_m_y_max, 0)
            designed_pier['phiMu Y Min (kNm)'] = round(column_check_results.phi_m_y_min, 0)
            designed_pier['Utilisation Y Max'] = round(utilisation_y_1, 2)
            designed_pier['Utilisation Y Min'] = round(utilisation_y_2, 2)

        # Calculate boundary element (if required)
        if self.design_boundary_element:
            be_length, lig_size, lig_spacing = vertical_structure_v1.boundary_element(pier_section, loading)
            designed_pier['BE Length (mm)'] = be_length
            designed_pier['Lig Size (mm)'] = lig_size
            designed_pier['Lig Spacing (mm)'] = lig_spacing

        return designed_pier


def design_etabs_pier_as_column(pier, eq_env_1, eq_env_2, wind_env, vertical_spacing, horizontal_spacing, design_both_axes, design_boundary_element):
    designer = PierColumnDesigner(pier, eq_env_1, eq_env_2, wind_env, vertical_spacing, horizontal_spacing, design_both_axes, design_boundary_element)
    return designer.design_pier()


def design_all_piers(
        piers:list[dict], 
        eq_env_1:str, 
        eq_env_2:str, 
        wind_env:str, 
        vertical_spacing:float, 
        horizontal_spacing:float,
        design_both_axes:bool,
        design_boundary_element:bool
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
            design_boundary_element
            )
        designed_piers.append(designed_pier)

    # Create a dataframe to store the results
    df = pd.DataFrame(designed_piers)
    return df


def filter_bar_sizes(fc: float, story: str, story_names: list, phz_levels: list, vertical_spacing: float, b: float):

    # Define bar sizes
    bar_sizes = [12, 16, 20, 24, 28, 32, 36, 40]
    
    # Calculate minimum tension reinforcement
    rho_wv_crit, rho_wv_typ = vertical_structure_v1.min_tension_reinforcement(fc, story, story_names, phz_levels)
    rho_max = max(rho_wv_crit, rho_wv_typ)
    
    # Filter bar sizes based on rho calculation
    filtered_bar_sizes = []
    for bar_size in bar_sizes:
        rho = 2 * (math.pi * (bar_size ** 2) / 4) * (1000 / vertical_spacing) / (1000 * b)
        if rho > rho_max:
            filtered_bar_sizes.append(bar_size)
    
    return filtered_bar_sizes

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
