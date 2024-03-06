import os
import sys
sys.path.append(r"C:\Program Files\Bentley\Engineering\RAM Concept CONNECT Edition\RAM Concept CONNECT Edition V8\python")
import ram_concept
import pandas as pd

from ram_concept.concept import Concept
from ram_concept.cad_manager import CadManager
from ram_concept.result_layers import ReactionContext
from ram_concept.point_2D import Point2D

# open model
reaction_concept = Concept.start_concept(headless=True)
reaction_model_path = r"C:\_Local Projects\ram_lrd_test\WD21009 2023-09-05 _679-683  Glen Huntly Rd_ L1.cpt"
reaction_model = reaction_concept.open_file(reaction_model_path)

# run model
reaction_model.calc_all()

# prepare layers
cad_manager = reaction_model.cad_manager
element_layer = cad_manager.element_layer

# loading and combo layers
loadings = cad_manager.force_loading_layers
load_combos = cad_manager.load_combo_layers

# combine the load combo layers and the loading layers, as we want reactions for both
loadings_and_combos = loadings + load_combos

dead_load_reaction_name = "All Dead LC"
live_load_reaction_name = "All Live Lc"
reaction_cases = [dead_load_reaction_name, live_load_reaction_name]

column_reactions = []
wall_reactions = []

# loop through all the loadings and combos and get their reactions
for loading in loadings_and_combos:
    if loading.name in reaction_cases:
        for column_element in element_layer.column_elements_below:
            reaction = loading.column_reaction(column_element, ReactionContext.STANDARD)
            column_reactions.append({
                'Column No.': column_element.number,
                'Location X': column_element.location.x,
                'Location Y': column_element.location.y,
                'Reaction Case': loading.name,
                'Fx': reaction.x,
                'Fy': reaction.y,
                'Fz': reaction.z,
                'Rx': reaction.rot_x,
                'Ry': reaction.rot_y
            })
        
        for wall_element in element_layer.wall_element_groups_below:
            reaction = loading.wall_group_reaction(wall_element, ReactionContext.STANDARD)
            wall_reactions.append({
                'Wall No.': wall_element.number,
                'Centroid X': wall_element.centroid.x,
                'Centroid Y': wall_element.centroid.y,
                'Reaction Case': loading.name,
                'Fx': reaction.x,
                'Fy': reaction.y,
                'Fz': reaction.z,
                'Rx': reaction.rot_x,
                'Ry': reaction.rot_y
            })

reaction_concept.shut_down() # close model

# add reactions to target model
target_concept = Concept.start_concept(headless=True)
target_model_path = r"C:\_Local Projects\ram_lrd_test\WD21009 2023-09-05 _679-683  Glen Huntly Rd_ GF.cpt"
target_model = target_concept.open_file(target_model_path)

cad_manager = target_model.cad_manager

dead_load_target_name = "Level Above DL"
live_load_target_name = "Level Above LL"

dead_loading_target_layer = cad_manager.force_loading_layer(dead_load_target_name)
live_loading_target_layer = cad_manager.force_loading_layer(live_load_target_name)
target_layers = [dead_loading_target_layer, live_loading_target_layer]

# add column point loads
for column_element in column_reactions:
    case_index = reaction_cases.index(column_element['Reaction Case'])
    target_layer = target_layers[case_index]
    point_load = target_layer.add_point_load(Point2D(column_element['Location X'], column_element['Location Y']))
    point_load.Fz = column_element['Fz']

# add column point loads
for wall_element in wall_reactions:
    case_index = reaction_cases.index(wall_element['Reaction Case'])
    target_layer = target_layers[case_index]
    point_load = target_layer.add_point_load(Point2D(wall_element['Centroid X'], wall_element['Centroid Y']))
    point_load.Fz = wall_element['Fz']

target_model.save_file(target_model_path)
target_concept.shut_down()