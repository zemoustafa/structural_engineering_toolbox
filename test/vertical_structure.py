'''
Module to design columns to Section 10 of AS3600:2018

'''

import math
from enum import Enum
from dataclasses import dataclass
from collections import namedtuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sectionproperties.pre.library import concrete_rectangular_section

from concreteproperties import Concrete, ConcreteLinear, ConcreteSection, RectangularStressBlock, SteelBar, SteelElasticPlastic

from concreteproperties.design_codes import AS3600
from concreteproperties.results import MomentInteractionResults
from shapely.geometry import LineString, Point, Polygon

"""
CONSTANTS

"""
DEFAULT_FSY = 500  # MPa
DEFUALT_DESIGN_CODE = AS3600()
DEFUALT_STEEL_MATERIAL = DEFUALT_DESIGN_CODE.create_steel_material()

"""
ENUMS

"""
class BracingType(Enum):
    BRACED = "Braced"
    UNBRACED = "Unbraced"

class SectionType(Enum):
    COLUMN = "Column"
    WALL = "Wall"

"""
DATACLASSES - DEFINE SECTION PROPERTIES AND LOADING

"""
@dataclass
class RectangularColumn:
    section_type: SectionType  # 'Column' or 'Wall'
    fc: float  # Concrete strength (MPa)
    d: float  # Depth (long direction for walls, width for columns) (mm)
    b: float  # Thickness (for walls) or width (for columns) (mm)
    h: float  # Height (mm)
    bracing_x: BracingType # 'Braced' or 'Unbraced'
    bracing_y: BracingType
    cover: float  # Cover to reinforcement (mm)
    v_bar_dia: float  # Vertical bar diameter (mm)
    h_bar_dia: float  # Horizontal bar diameter (mm), only used for walls
    v_bar_cts: float = None  # Vertical bar spacing (only for walls)
    h_bar_cts: float = None  # Horizontal bar spacing (only for walls)
    n_bars_x: int = None  # Vertical bar spacing (only for columns)
    n_bars_y: int = None # Horizontal bar spacing (only for columns)

    def __post_init__(self) -> None:
        """Initialize the number of bars in the x (b) and y (d) direction for walls."""
        if self.section_type == SectionType.WALL:
            if self.n_bars_x is None:
                self.n_bars_x = 2  # Fixed for walls
            if self.n_bars_y is None:
                self.n_bars_y = int(round((self.d - 2 * self.cover - 2 * self.h_bar_dia) / self.v_bar_cts, 0))
        
@dataclass
class Loading:
    n_star_compression: float  # Axial force in compression (+ve)
    n_star_tension: float      # Axial force in tension (+ve)
    m_x_top: float              # Major-axis moment at the top
    m_x_bot: float              # Major-axis moment at the bottom
    m_y_top: float              # Minor-axis moment at the top
    m_y_bot: float              # Minor-axis moment at the bottom
    v_star: float              # Shear force

"""
NAMED TUPLES - STORE RESULTS

"""
ColumnDesignResults = namedtuple(
    "ColumnDesignResults",
    ['result_x_c', 'phi_n_x_c', 'phi_m_x_c', 'result_x_t', 'phi_n_x_t', 'phi_m_x_t', 'result_y_c', 'phi_n_y_c', 'phi_m_y_c', 'result_y_t', 'phi_n_y_t', 'phi_m_y_t'],
    defaults=(None, None, None, None, None, None, None, None, None, None, None, None)
)

# Store results from the moment interaction analysis
ColumnMIResults = namedtuple('ColumnMIResults', ['f_m_x', 'f_n_x', 'f_m_y', 'f_n_y', 'balance_point_x', 'balance_point_y'], defaults=(None, None, None, None, None, None))

# Store results from column capacity check
ColumnCapacityResults = namedtuple('ColumnCapacity', ['phi_m_x_comp', 'phi_n_x_comp', 'phi_m_x_tens', 'phi_n_x_tens', 'phi_m_y_comp', 'phi_n_y_comp', 'phi_m_y_tens', 'phi_n_y_tens'], defaults=[None, None, None, None, None, None, None, None])

"""
FUNCTIONS

"""
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
    extended_line = LineString([(0, 0), (100000 * direction[0], 100000 * direction[1])])  # Extend by a factor of 1000 or more

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
    f_m_x_copy = f_m_x[:] + [f_m_x[0]]
    f_n_x_copy = f_n_x[:] + [f_n_x[0]]

    curve_coords = list(zip(f_m_x_copy, f_n_x_copy)) # Combine x and y coordinates into a list of tuples

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


def __calculate_magnified_moment(section:RectangularColumn, loading:Loading, slenderness, bracing, balance_point, m_star, r, axis):
    if slenderness == 'Slender':
        n_c = column_buckling_load(section.d if axis == 'x' else section.b, section.h, section.v_bar_dia, section.h_bar_dia, section.cover, balance_point[0], beta_d=0.5)
        delta = moment_magnification_factor(section.fc, bracing, section.d if axis == 'x' else section.b, section.b if axis == 'x' else section.d, n_c, loading.n_star_compression, loading.m_x_bot if axis == 'x' else loading.m_y_bot, loading.m_x_bot if axis == 'x' else loading.m_y_bot, section.h, r, beta_d=0.5)
        return m_star * delta
    else:
        return m_star


def __calculate_effective_shear_depth(d, cover, v_bar_dia, v_bar_cts, h_bar_dia):
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
    v_bar_area = math.pi * v_bar_dia**2 / 4 # mm2
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


def shear_induced_tension_reinforcement(d:float, m_star:float, n_star:float, v_star:float, phi_vuc:float):
    
    theta_v_rads = math.radians(36)
    f_tds = 0.5 * ( v_star + phi_vuc)*1000 * (1 / math.tan(theta_v_rads) )
    
    z = d / 4

    t_td = ( m_star * 1e6 / z ) + ( n_star*1000 / 2 ) + f_tds

    a_st = t_td * 0.85 / DEFAULT_FSY

    return a_st


def column_shear(section:RectangularColumn):
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

    Returns:
    v_uc (float): concrete contribution to shear strength
    v_us (float): steel contribution to shear strength

    '''
    # calculate reo area and reo rate in section
    fsy = 500 # MPa
    h_bar_area = math.pi * section.h_bar_dia**2 / 4 # mm2
    v_bar_area = math.pi * section.v_bar_dia**2 / 4 # mm2
    a_ct = section.d * section.b / 2 # mm2 - area of concrete in tensile region
    n_bars_v = __round_up_to_even(2 * ( section.d - 2 * section.cover - 2 * section.h_bar_dia ) / section.v_bar_cts) # no. vertical bars in tensile region
    n_bars_h = round(2 * ( section.h - 2 * section.cover ) / section.h_bar_cts, 0) # no. horizontal bars in section

    rho = ( n_bars_v * v_bar_area / (section.d * section.b) ) * 100 # reinforcement rate

    # Concrete contribution to shear cl 8.2.4
    dvy, ast_total = __calculate_effective_shear_depth(section.d, section.cover, section.v_bar_dia, section.v_bar_cts, section.h_bar_dia)
    dvs = section.d / 2 - dvy
    # strain_x = ( (m_star*1000000)/dvy + (v_star*1000) ) / ( 2*( dvs/ dvy * 200000 * ast_total))

    # Simple method for kv and theta v
    theta_v = 36 # angle of inclination of concrete strut
    theta_v_rads = math.radians(theta_v)
    A_sv_min = section.h_bar_cts * 0.08 * math.sqrt(section.fc) * section.b / fsy  # Minimum transverse shear reinforcement - Cl 8.2.1.7
    A_sv = h_bar_area * n_bars_h
    kv = 200/(100 + 1.3*section.d) if A_sv/section.h_bar_cts < A_sv_min/section.h_bar_cts else 0.15

    # concrete contribution to shear strength
    v_uc = round(kv * section.b * dvy * min(math.sqrt(section.fc), 8) / 1000, 0)

    # Transverse shear and torsion reinforcement contribution - Cl 8.2.5
    v_us = ((2 * h_bar_area) * fsy * dvy / section.h_bar_cts) * (1/math.tan(theta_v_rads)) / 1000

    # return values
    return v_uc, v_us


def moment_magnification_factor(fc:float, bracing:str, d:float, b:float, n_c:float, n_star:float, m_star_top:float, m_star_bot:float, l_e:float, r:float, beta_d:float=0.5):
    '''
    Determines moment amplification factor for slender columns based on Cl 10.4

    Parameters:
    bracing (str): bracing condition of column. must be either 'Braced' or 'Unbraced'
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


def column_buckling_load(d, h, cover, v_bar_dia, h_bar_dia, m_ub, beta_d=0.5):
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


def moment_interaction_design(
        section:RectangularColumn,
        bracing:str,
        loading:Loading,
        check_both_axes:float
        ):
    '''
    Checks column using moment interaction diagram to Section 10
    
    Parameters:
    section (RectangularColumn): ConcreteSection object
    bracing (str): bracing condition of column. Value must be 'Braced' or 'Unbraced'
    loading (Loading): ConcreteLoading object
    check_both_axes (bool): True if weak axis is to be checked, otherwise False

    Returns:
    results (MomentInteractionResults): Moment interaction results

    '''
    # Check tension load is negative
    if loading.n_star_tension is not None and loading.n_star_tension > 0:
        raise ValueError("Error: Tension must be a negative value.")

    def create_concrete_section(b, d, n_bars_x, n_bars_y):
        bar_area = math.pi * section.v_bar_dia**2 / 4
        geom = concrete_rectangular_section(
            b=b, d=d,
            dia_top=section.v_bar_dia, area_top=bar_area, n_top=n_bars_x, c_top=section.cover + section.h_bar_dia,
            dia_bot=section.v_bar_dia, area_bot=bar_area, n_bot=n_bars_x, c_bot=section.cover + section.h_bar_dia,
            dia_side=section.v_bar_dia, area_side=bar_area, n_side=n_bars_y - 2, c_side=section.cover + section.h_bar_dia,
            conc_mat=concrete, steel_mat=DEFUALT_STEEL_MATERIAL,
        )
        return ConcreteSection(geom)

    def get_moment_interaction_results(conc_sec):
        DEFUALT_DESIGN_CODE.assign_concrete_section(concrete_section=conc_sec)
        return DEFUALT_DESIGN_CODE.moment_interaction_diagram(progress_bar=False, n_points=18, control_points=[("fy", 1.0)])

    def process_results(f_mi_res):
        f_results_list = f_mi_res.get_results_lists(moment='m_x')
        f_n = [x / 1000 for x in f_results_list[0]]  # Convert to kN
        f_m = [x / 1000000 for x in f_results_list[1]]  # Convert to kNm
        return f_n, f_m

    def calculate_magnified_moment(slenderness, h, mi_res, m_star, r, axis):
        if slenderness == 'Slender':
            balance_point = __extract_balance_point(mi_res)
            n_c = column_buckling_load(section.d if axis == 'x' else section.b, h, section.v_bar_dia, section.h_bar_dia, section.cover, balance_point[0], beta_d=0.5)
            delta = moment_magnification_factor(section.fc, bracing, section.d if axis == 'x' else section.b, section.b if axis == 'x' else section.d, n_c, loading.n_star_compression, loading.m_x_bot if axis == 'x' else loading.m_y_bot, loading.m_x_bot if axis == 'x' else loading.m_y_bot, h, r, beta_d=0.5)
            return m_star * delta
        return m_star

    # Determine if column is short or slender - Cl 10.2
    rx = 0.3 * section.d # Radius of gyration - Cl 10.5.2
    slenderness_x = 'Short' if section.h / rx <= 22 else 'Slender'

    # Set up moment interaction diagram using concreteproperties library
    concrete = DEFUALT_DESIGN_CODE.create_concrete_material(compressive_strength=section.fc)

    # Create and process moment interaction results for both axes
    conc_sec_x = create_concrete_section(section.b, section.d, section.n_bars_x, section.n_bars_y)
    f_mi_res_xx, mi_res_xx, _ = get_moment_interaction_results(conc_sec_x)
    f_n_x, f_m_x = process_results(f_mi_res_xx)

    # Determine max applied moment in each direction
    m_star_x = max(abs(loading.m_x_top), abs(loading.m_x_bot))

    # Apply moment magnification if column is slender
    m_star_x = calculate_magnified_moment(slenderness_x, section.h, mi_res_xx, m_star_x, rx, 'x')

    # Check applied axial load and moment fall within diagram (compression)
    point_in_diagram_x_comp = __is_point_inside_curve(f_m_x, f_n_x, loading.n_star_compression, m_star_x)
    pass_fail_x_comp = 'Pass' if point_in_diagram_x_comp else 'Fail'
    
    # Find the intersection points
    phi_n_x_comp, phi_m_x_comp = __find_intersection(f_m_x, f_n_x, m_star_x, loading.n_star_compression)
    
    # If there is tension, find intersection
    if loading.n_star_tension == 0:
        capacity_x = (phi_n_x_comp, phi_m_x_comp)
        pass_fail_x = pass_fail_x_comp
    else:
        # Run same check but for tensile load (only if tensile load exists)
        point_in_diagram_x_tens = __is_point_inside_curve(f_m_x, f_n_x, loading.n_star_tension, m_star_x)
        pass_fail_x_tens = 'Pass' if point_in_diagram_x_tens else 'Fail'
    
        # Find the intersection points
        phi_n_x_tens, phi_m_x_tens = __find_intersection(f_m_x, f_n_x, m_star_x, loading.n_star_tension)
    
        capacity_x = (phi_n_x_comp, phi_m_x_comp, phi_n_x_tens, phi_m_x_tens)
        pass_fail_x = (pass_fail_x_comp, pass_fail_x_tens)


    # Repeat all calculations if check_weak_axis is checked
    if check_both_axes == True:
        ry = 0.3 * section.b  # Radius of gyration - Cl 10.5.2
        slenderness_y = 'Short' if section.h / ry <= 22 else 'Slender' # Determine if column is short or slender - Cl 10.2
        conc_sec_y = create_concrete_section(section.d, section.b, section.n_bars_y, section.n_bars_x) # Create concrete section about y axis
        f_mi_res_yy, mi_res_yy, _ = get_moment_interaction_results(conc_sec_y) # Moment interaction results about y axis
        f_n_y, f_m_y = process_results(f_mi_res_yy)
        m_star_y = max(abs(loading.m_y_top), abs(loading.m_y_bot)) # Grab max design moment about y axis
        m_star_y = calculate_magnified_moment(slenderness_y, section.h, mi_res_yy, m_star_y, ry, 'y') # Calculate magnified moment about y axis
        point_in_diagram_y_comp = __is_point_inside_curve(f_m_y, f_n_y, loading.n_star_compression, m_star_y) # Determine if point inside diagram
        pass_fail_y_comp = 'Pass' if point_in_diagram_y_comp else 'Fail'
        phi_n_y_comp, phi_m_y_comp = __find_intersection(f_m_y, f_n_y, m_star_y, loading.n_star_compression)
        # If tension > 0, find intersection
        if loading.n_star_tension == 0:

            results = ColumnDesignResults(
                result_x_c=pass_fail_x_comp,
                phi_n_x_c=phi_n_x_comp,
                phi_m_x_c=phi_m_x_comp,
                result_y_c=pass_fail_y_comp,
                phi_n_y_c=phi_n_y_comp,
                phi_m_y_c=phi_m_y_comp
            )
            
        else:
            # Run same check but for tensile load (only if tensile load exists)
            point_in_diagram_y_tens = __is_point_inside_curve(f_m_y, f_n_y, loading.n_star_tension, m_star_y)
            pass_fail_y_tens = 'Pass' if point_in_diagram_y_tens else 'Fail'
        
            # Find the intersection points
            phi_n_y_tens, phi_m_y_tens = __find_intersection(f_m_y, f_n_y, m_star_y, loading.n_star_tension)
        
            results = ColumnDesignResults(
                result_x_c=pass_fail_x_comp,
                phi_n_x_c=phi_n_x_comp,
                phi_m_x_c=phi_m_x_comp,
                result_x_t=pass_fail_x_tens,
                phi_n_x_t=phi_n_x_tens,
                phi_m_x_t=phi_m_x_tens,
                result_y_c=pass_fail_y_comp,
                phi_n_y_c=phi_n_y_comp,
                phi_m_y_c=phi_m_y_comp,
                result_y_t=pass_fail_y_tens,
                phi_n_y_t=phi_n_y_tens,
                phi_m_y_t=phi_m_y_tens
            )
    else:
        results = (pass_fail_x, capacity_x)
        results = ColumnDesignResults(
            result_x_c=pass_fail_x_comp,
            phi_n_x_c=phi_n_x_comp,
            phi_m_x_c=phi_m_x_comp,
            result_x_t=pass_fail_x_tens if loading.n_star_tension != 0 else None,
            phi_n_x_t=phi_n_x_tens if loading.n_star_tension != 0 else None,
            phi_m_x_t=phi_m_x_tens if loading.n_star_tension != 0 else None
        )

    return results
    # add fig_x and fig_y to results if you want to show plot again


def rectangular_column_moment_interaction(section:RectangularColumn, both_axes:bool):
    '''
    Calculates column interaction diagram for a rectangular column to AS3600:2018 Section 10

    Parameters:
    section (RectangularColumn): ConcreteSection object
    check_both_axes (bool): True if weak axis is to be checked, otherwise False

    Returns:
    results (MomentInteractionResults): Moment interaction results

    '''
    # Initialise results object
    results = ColumnMIResults()

    # Define concrete material property based on f'c defined in section
    concrete = DEFUALT_DESIGN_CODE.create_concrete_material(compressive_strength=section.fc)
    
    # Calculate vertical bar area
    bar_area = math.pi * section.v_bar_dia**2 / 4

    # Create concrete section
    def create_concrete_section(b, d, n_bars_x, n_bars_y):
        geom = concrete_rectangular_section(
            b=b, d=d,
            dia_top=section.v_bar_dia, area_top=bar_area, n_top=n_bars_x, c_top=section.cover + section.h_bar_dia,
            dia_bot=section.v_bar_dia, area_bot=bar_area, n_bot=n_bars_x, c_bot=section.cover + section.h_bar_dia,
            dia_side=section.v_bar_dia, area_side=bar_area, n_side=n_bars_y - 2, c_side=section.cover + section.h_bar_dia,
            conc_mat=concrete, steel_mat=DEFUALT_STEEL_MATERIAL,
        )
        return ConcreteSection(geom)
    
    # Function to get moment interaction results
    def get_moment_interaction_results(conc_sec):
        DEFUALT_DESIGN_CODE.assign_concrete_section(concrete_section=conc_sec)
        return DEFUALT_DESIGN_CODE.moment_interaction_diagram(progress_bar=False, n_points=18, control_points=[("fy", 1.0)])
    
    # Function to convert moment interaction results to kN and kNm
    def process_results(f_mi_res):
        f_results_list = f_mi_res.get_results_lists(moment='m_x')
        f_n = [x / 1000 for x in f_results_list[0]]  # Convert to kN
        f_m = [x / 1000000 for x in f_results_list[1]]  # Convert to kNm
        return f_n, f_m

    # Define concrete section (from concreteproperties library)
    conc_sec_x = create_concrete_section(section.b, section.d, section.n_bars_x, section.n_bars_y)

    # Create moment interaction diagram and convert to kN and kNm
    f_mi_res_x, mi_res_x, _ = get_moment_interaction_results(conc_sec_x)
    f_n_x, f_m_x = process_results(f_mi_res_x) # Convert units

    # Grab balance point
    balance_point_x = __extract_balance_point(mi_res_x)

    # Update results
    results = ColumnMIResults(f_m_x=f_m_x, f_n_x=f_n_x, balance_point_x=balance_point_x)

    # If both_axes is true, check and return moment interaciton results about y axis
    if both_axes == True:
        # Define concrete section (from concreteproperties library). Flip dimensions for y axis
        conc_sec_y = create_concrete_section(section.d, section.b, section.n_bars_y, section.n_bars_x)

        # Create moment interaction diagram and convert to kN and kNm
        f_mi_res_y, mi_res_y, _ = get_moment_interaction_results(conc_sec_y)
        f_n_y, f_m_y = process_results(f_mi_res_y) # Convert units
        
        # Grab balance point
        balance_point_y = __extract_balance_point(mi_res_y)

        # Update results
        results = ColumnMIResults(f_m_y=f_m_y, f_n_y=f_n_y, balance_point_y=balance_point_y)

    return results


def check_column_capacity(section:RectangularColumn, loading:Loading, mi_results:ColumnMIResults):
    '''
    Checks column capacity based on loading and moment interaction diagram.

    Parameters:
    section (RectangularColumn): ConcreteSection object
    loading (Loading): Loading object
    mi_results (ColumnMIResults): Results of moment interaction diagram

    Returns:
    capacity (ColumnCapacityResults): Named tuple containing column capacities
    
    '''

    # Initialise results object
    results = ColumnCapacityResults()
    
    # Determine if column is short or slender - Cl 10.2
    rx = 0.3 * section.d # Radius of gyration - Cl 10.5.2
    slenderness_x = 'Short' if section.h / rx <= 22 else 'Slender'

    # Determine max applied moment in each direction
    m_star_x = max(abs(loading.m_x_top), abs(loading.m_x_bot))

    # Apply moment magnification if column is slender
    m_star_x = __calculate_magnified_moment(section, loading, slenderness_x, section.bracing_x, mi_results.balance_point_x, m_star_x, rx,  'x')

    # Find the intersection points
    phi_n_x_comp, phi_m_x_comp = __find_intersection(mi_results.f_m_x, mi_results.f_n_x, m_star_x, loading.n_star_compression)
    
    # If there is tension, find intersection
    if loading.n_star_tension == 0:
        results.phi_m_x_comp = phi_m_x_comp
        results.phi_n_x_comp = phi_n_x_comp

    else:
        # Run same check but for tensile load (only if tensile load exists)
        phi_n_x_tens, phi_m_x_tens = __find_intersection(mi_results.f_m_x, mi_results.f_n_x, m_star_x, loading.n_star_tension)
        results.phi_m_x_tens = phi_m_x_tens
        results.phi_n_x_tens = phi_n_x_tens

    # If y direction loading exists, repeat calculations
    if hasattr(mi_results, 'f_m_y') and hasattr(loading, 'm_y_bot'):
        # Determine if column is short or slender - Cl 10.2
        ry = 0.3 * section.b # Radius of gyration - Cl 10.5.2
        slenderness_y = 'Short' if section.h / ry <= 22 else 'Slender'

        # Determine max applied moment in each direction
        m_star_y = max(abs(loading.m_y_top), abs(loading.m_y_bot))

        # Apply moment magnification if column is slender
        m_star_y = __calculate_magnified_moment(slenderness_y, section.h, mi_results.m_y, m_star_y, ry, 'y')

        # Find the intersection points
        phi_n_y_comp, phi_m_y_comp = __find_intersection(mi_results.f_m_y, mi_results.f_n_y, m_star_y, loading.n_star_compression)

        # If tension > 0, find intersection
        if loading.n_star_tension == 0:
            results.phi_m_y_comp = phi_m_y_comp
            results.phi_n_y_comp = phi_n_y_comp

        else:
            # Run same check but for tensile load (only if tensile load exists)
            phi_n_y_tens, phi_m_y_tens = __find_intersection(mi_results.f_m_y, mi_results.f_n_y, m_star_y, loading.n_star_tension)
            results.phi_m_y_tens = phi_m_y_tens
            results.phi_n_y_tens = phi_n_y_tens

    return results













# # # input parameters
# wall = RectangularColumn(
#     section_type=SectionType.WALL,
#     fc=50,
#     d=4000,
#     b=300,
#     h=4200,
#     cover=50,
#     v_bar_dia=20,
#     v_bar_cts=200,
#     h_bar_dia=16,
#     h_bar_cts=200
#     )

# loading = Loading(
#     n_star_compression=4000,
#     n_star_tension=0,
#     m_x_top=-3000,
#     m_x_bot=0,
#     m_y_top=0,
#     m_y_bot=400,
#     v_star=5000
# )

# check_both_axes = False
# # Check column interaction diagram
# results = moment_interaction_design(
#     section=wall,
#     bracing='Unbraced',
#     loading=loading,
#     check_both_axes=check_both_axes
# )

# # Check column shear capacity
# v_uc, v_us = wall_shear(wall)

# print('Column Design Results')
# print('---------------------')
# print('Column Dimensions: ' + str(wall.b) + ' x ' + str(wall.d) + ' mm')
# print('Concrete Strength: ' + str(wall.fc) + ' MPa')
# print('Vertical Bar Diameter: ' + str(wall.v_bar_dia) + ' mm')
# print('Vertical Bar Spacing: ' + str(wall.v_bar_cts) + ' mm')
# print('Horizontal Bar Diameter: ' + str(wall.h_bar_dia) + ' mm')
# print('Horizontal Bar Spacing: ' + str(wall.h_bar_cts) + ' mm')
# print('Numnber of bars x: ' + str(wall.n_bars_x))
# print('Numnber of bars y: ' + str(wall.n_bars_y))
# print('Cover to Reinforcement: ' + str(wall.cover) + ' mm')
# print('Axial Load (compression): ' + str(loading.n_star_compression) + ' kN')
# print('Axial Load (tension): ' + str(loading.n_star_tension) + ' kN')
# print('---------------------')
# print('Results X')
# print('X compression: '+ results.result_x_c) # pass or fail
# print('Phi N X C: ' + str(round(results.phi_n_x_c))) # intersection with N
# print('Phi M X C: ' + str(round(results.phi_m_x_c))) # intersection with M

# if loading.n_star_tension < 0:
#     print('X tension: '+ results.result_x_t) # pass or fail
#     print('Phi N X T: ' + str(round(results.phi_n_x_t))) # intersection with N
#     print('Phi M X T: ' + str(round(results.phi_m_x_t))) # intersection with M

# if check_both_axes == True:
#     print('---------------------')
#     print('Results Y')
#     print('Y compression: ' + results.result_y_c) # pass or fail
#     print('Y tension: ' + results.result_y_t) # pass or fail
#     print('Phi N Y C: ' + str(round(results.phi_n_y_c))) # intersection with N
#     print('Phi M Y C: ' + str(round(results.phi_m_y_c))) # intersection with M
# if check_both_axes == True and loading.n_star_tension < 0:
#     print('Phi N Y T: ' + str(round(results.phi_n_y_t))) # intersection with N
#     print('Phi M Y T: ' + str(round(results.phi_m_y_t))) # intersection with M

# shear_pass_fail = 'Pass' if 0.7*v_uc + 0.7*v_us > loading.v_star else 'Fail'
# print('Shear Capacity Results: ' + shear_pass_fail)
# print('Concrete Shear Capacity: ' + str(round(v_uc, 1)) + ' kN') # concrete shear capacity
# print('Steel Shear Capacity: ' + str(round(v_us, 1)) + ' kN') # steel shear capacity

# plt.show()
