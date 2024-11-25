import os
import sys
import math
sys.path.append(r"C:\Program Files\Bentley\Engineering\RAM Concept\RAM Concept 2023\python")
import ram_concept

from ram_concept.concept import Concept
from ram_concept.cad_manager import CadManager
from ram_concept.result_layers import ReactionContext
from ram_concept.point_2D import Point2D
from ram_concept.line_segment_2D import LineSegment2D

def copy_object():
    pass

def full_load_rundown(file_path_list):
    """
    Complete a full load rundown, requires file path of every level and
    for each model to be aligned

    """
    for i in range(len(file_path_list) - 1):
        delete_exisitng_loads(file_path_list[i + 1])

        reaction_path = file_path_list[i]
        target_path = file_path_list[i + 1]
        ram_load_rundown(reaction_path, target_path)

def delete_exisitng_loads(model_path):
    '''
    Deletes all existing point and line loads that have a name
    
    '''
    concept = Concept.start_concept(headless=True)
    model = concept.open_file(model_path)

    dead_load_name = "Other Dead Loading"
    live_load_name = "Live (Unreducible) Loading"
    case_names = [dead_load_name, live_load_name]

    # prepare layers
    cad_manager = model.cad_manager
    element_layer = cad_manager.element_layer
    loadings = cad_manager.force_loading_layers

    for loading in loadings:
        if loading.name in case_names:
            line_loads = loading.line_loads
            point_loads = loading.point_loads
            for line_load in line_loads:
                if line_load.name != '':
                    line_load.delete()
            for point_load in point_loads:
                if point_load.name != '':
                    point_load.delete()

    model.save_file(model_path)
    concept.shut_down()

def ram_load_rundown(reaction_path, target_path):
    '''
    Grabs reaction from one model (reaction) and applies loads to another (target) 
    
    '''
    # open model
    reaction_concept = Concept.start_concept(headless=True)
    reaction_model_path = reaction_path
    reaction_model = reaction_concept.open_file(reaction_model_path)

    # run model
    reaction_model.generate_mesh()
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
    live_load_reaction_name = "Live (Unreducible) Loading"
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
                    'Length': wall_element.total_length,
                    'Angle': math.radians(wall_element.reaction_angle),
                    'Fx': reaction.x,
                    'Fy': reaction.y,
                    'Fz': reaction.z,
                    'Rx': reaction.rot_x,
                    'Ry': reaction.rot_y
                })

    reaction_concept.shut_down() # close model

    # add reactions to target model
    target_concept = Concept.start_concept(headless=True)
    target_model_path = target_path
    target_model = target_concept.open_file(target_model_path)

    cad_manager = target_model.cad_manager

    dead_load_target_name = "Other Dead Loading"
    live_load_target_name = "Live (Unreducible) Loading"

    dead_loading_target_layer = cad_manager.force_loading_layer(dead_load_target_name)
    live_loading_target_layer = cad_manager.force_loading_layer(live_load_target_name)
    target_layers = [dead_loading_target_layer, live_loading_target_layer]

    # add column point loads
    column_name_index = 0
    for column_element in column_reactions:
        case_index = reaction_cases.index(column_element['Reaction Case'])
        target_layer = target_layers[case_index]
        column_element['Point Load'] = target_layer.add_point_load(Point2D(column_element['Location X'], column_element['Location Y']))
        column_element['Point Load'].Fz = column_element['Fz']
        column_element['Point Load'].name = "column" + str(column_name_index)
        column_name_index = column_name_index + 1

    # add wall reactions as line loads
    wall_name_index = 0
    for wall_element in wall_reactions:
        case_index = reaction_cases.index(wall_element['Reaction Case'])
        target_layer = target_layers[case_index]
        x = wall_element['Length']/2 * math.cos(wall_element['Angle'])
        y = wall_element['Length']/2 * math.sin(wall_element['Angle'])
        start_point = Point2D(wall_element['Centroid X'] - x, wall_element['Centroid Y'] - y)
        end_point = Point2D(wall_element['Centroid X'] + x, wall_element['Centroid Y'] + y)
        line_segment = LineSegment2D(start_point, end_point)
        wall_element['Line Load'] = target_layer.add_line_load(line_segment)
        wall_element['Line Load'].Fz0 = wall_element['Fz'] / wall_element['Length']
        wall_element['Line Load'].Fz1 = wall_element['Fz'] / wall_element['Length']
        wall_element['Line Load'].name = "wall" + str(wall_name_index)
        wall_name_index = wall_name_index + 1

    target_model.save_file(target_model_path)
    target_concept.shut_down()



    # # add wall reactions as point loads
    # for wall_element in wall_reactions:
    #     case_index = reaction_cases.index(wall_element['Reaction Case'])
    #     target_layer = target_layers[case_index]
    #     point_load = target_layer.add_point_load(Point2D(wall_element['Centroid X'], wall_element['Centroid Y']))
    #     point_load.Fz = wall_element['Fz']