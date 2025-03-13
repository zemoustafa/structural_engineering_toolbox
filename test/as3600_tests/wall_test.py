import sys
sys.path.append('C:\_Github\structural_engineering_toolbox')
from as3600.vertical_structure_v1 import *

wall = RectangularColumn(
    section_type=SectionType.WALL,
    fc=50,
    d=6000,
    b=300,
    h=4200,
    cover=30,
    v_bar_dia=20,
    v_bar_cts=200,
    h_bar_dia=16,
    h_bar_cts=200,
    bracing_x=BracingType.UNBRACED,
    bracing_y=BracingType.BRACED,
    )

loading = Loading(
    n_star_max=10000,
    n_star_min=4000,
    m_x_top=8000,
    m_x_bot=0,
    m_y_top=500,
    m_y_bot=0,
    v_star=5000
)

check_both_axes = False

# Create interaction diagram for column
results = rectangular_column_moment_interaction(section=wall, design_both_axes=check_both_axes)

# Calculate capacity
column_capacity = check_column_capacity(wall, loading, results)

# Wall shear capacity
v_uc, v_us = column_shear_capacity(wall)
shear_reo = shear_induced_tension_reinforcement(wall.d, m_star=loading.m_x_top, n_star=loading.n_star_min, v_star=loading.v_star, phi_vuc = 0.7*v_uc)

# Boundary element
be_length, lig_size, lig_cts = boundary_element(wall=wall, loading=loading)

print('Column Design Results')
print('---------------------')
print('Wall Dimensions: ' + str(wall.b) + ' x ' + str(wall.d) + ' mm')
print('Concrete Strength: ' + str(wall.fc) + ' MPa')
print('Vertical Bar Diameter: ' + str(wall.v_bar_dia) + ' mm')
print('Vertical Bar Spacing: ' + str(wall.v_bar_cts) + ' mm')
print('Horizontal Bar Diameter: ' + str(wall.h_bar_dia) + ' mm')
print('Horizontal Bar Spacing: ' + str(wall.h_bar_cts) + ' mm')
print('Numnber of bars x: ' + str(wall.n_bars_x))
print('Numnber of bars y: ' + str(wall.n_bars_y))
print('Cover to Reinforcement: ' + str(wall.cover) + ' mm')
# print('Axial Load (max): ' + str(loading.n_star_max) + ' kN')
# print('Axial Load (min): ' + str(loading.n_star_min) + ' kN')
print('---------------------')
print('Buckling Load: ' + str(round(column_capacity.n_c_x, 0)) + ' kN')
print('Is Buckling Load Exceeded: ' + str(column_capacity.buckling_x))

print('Results X')
print('Phi N X Max: ' + str(round(column_capacity.phi_n_x_max, 0))) # intersection with N
print('Phi M X Max: ' + str(round(column_capacity.phi_m_x_max, 0))) # intersection with M

if column_capacity.phi_m_x_min is not None:
    print('Phi N X Min: ' + str(round(column_capacity.phi_n_x_min))) # intersection with N
    print('Phi M X Min: ' + str(round(column_capacity.phi_m_x_min))) # intersection with M

if column_capacity.phi_m_y_max is not None:
    print('---------------------')
    print('Results Y')
    print('Buckling Load Y: ' + str(round(column_capacity.n_c_y, 0)) + ' kN')
    print('Phi N Y Max: ' + str(round(column_capacity.phi_n_y_max))) # intersection with N
    print('Phi M Y Max: ' + str(round(column_capacity.phi_m_y_max))) # intersection with M

if column_capacity.phi_m_y_min is not None:
    print('Phi N Y Min: ' + str(round(column_capacity.phi_n_y_min))) # intersection with N
    print('Phi M Y Min: ' + str(round(column_capacity.phi_m_y_min))) # intersection with M

shear_pass_fail = 'Pass' if 0.7*v_uc + 0.7*v_us > loading.v_star else 'Fail'
print('---------------------')
print('Shear Capacity Results: ' + shear_pass_fail)
print('Concrete Shear Capacity: ' + str(round(v_uc, 1)) + ' kN') # concrete shear capacity
print('Steel Shear Capacity: ' + str(round(v_us, 1)) + ' kN') # steel shear capacity
print('Addiitonal reo = ' + str(shear_reo) + ' mm2')

print('---------------------')
print('Boundary Element')
print('Length: ' + str(be_length) + ' mm')
print('Lig Size: ' + str(lig_size) + ' mm')
print('Lig Cts: ' + str(lig_cts) + ' mm')