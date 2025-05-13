# Ensure these imports are at the top of your ram_api.py or ram_concept_tools.py file
import math
import os
import sys
sys.path.append(r"C:\Program Files\Bentley\Engineering\RAM Concept\RAM Concept 2023\python")
import ram_concept

from ram_concept.concept import Concept, Model
from ram_concept.column import Column, DefaultColumn
from ram_concept.wall import Wall, DefaultWall
from ram_concept.concrete import Concrete
from ram_concept.concretes import Concretes
from ram_concept.cad_manager import CadManager
from ram_concept.result_layers import ReactionContext
from ram_concept.point_2D import Point2D
from ram_concept.line_segment_2D import LineSegment2D

class RamConcept:
    """A class for interacting with the RAM Concept API."""

    def __init__(self):
        """Initializes the class with a Concept instance."""
        self.concept = None

    def create_concept(self) -> Concept:
        """Create instance of Concept class."""
        self.concept = Concept.start_concept(headless=True)

    def open_model(self, file_path: str) -> Model:
        """Open a RAM Concept model."""
        return self.concept.open_file(file_path)
    
    def shutdown_concept(self):
        """Shuts down the Concept instance."""
        self.concept.shut_down()

def ram_load_rundown(reaction_path:str, target_path:str, concept:Concept, progress_queue=None):
    '''
    Grabs reaction from one model (reaction) and applies loads to another (target).
    Sends progress updates to the provided queue.
    
    Args:
        reaction_path (str): Path to the RAM Concept file to get reactions from.
        target_path (str): Path to the RAM Concept file to apply loads to.
        progress_queue (queue.Queue, optional): Queue for sending progress messages.
                                                Defaults to None (prints to console).
    '''
    reaction_level = os.path.basename(reaction_path) # Extract names of level from path
    target_level = os.path.basename(target_path)    

    # Helper function to send progress messages
    def _send_progress(message):
        if progress_queue:
            progress_queue.put(message)
        else:
            print(message) 

    _send_progress(f"Initiating load transfer from reaction file to target file\n")

    try:
        reaction_model = concept.open_file(reaction_path) # Open reaction model
        _send_progress(f"  Reaction model opened: '{reaction_level}'\n")
        _send_progress(f"  Generating mesh for reaction model...\n")

        reaction_model.generate_mesh()
        _send_progress(f"  Mesh generated for reaction model. Calculating all...\n")
        reaction_model.calc_all()
        _send_progress(f"  Calculations complete for reaction model.\n")

        cad_manager = reaction_model.cad_manager
        element_layer = cad_manager.element_layer
        loadings = cad_manager.force_loading_layers
        load_combos = cad_manager.load_combo_layers
        loadings_and_combos = loadings + load_combos

        dead_load_reaction_name = "All Dead LC"  
        live_load_reaction_name = "Live (Unreducible) Loading" 
        reaction_cases = [dead_load_reaction_name, live_load_reaction_name]
        
        _send_progress(f"  Preparing to extract reactions for cases: {', '.join(reaction_cases)}...\n")

        column_reactions = []
        wall_reactions = []

        for loading in loadings_and_combos:
            if loading.name in reaction_cases:
                _send_progress(f"    Extracting reactions from '{loading.name}' in {reaction_level}...\n")
                num_cols_processed = 0
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
                    num_cols_processed +=1
                _send_progress(f"      Processed {num_cols_processed} column reactions for '{loading.name}'.\n")

                num_walls_processed = 0
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
                    num_walls_processed += 1
                _send_progress(f"      Processed {num_walls_processed} wall reactions for '{loading.name}'.\n")
        
        _send_progress(f"  Finished extracting all reactions from '{reaction_level}'. Total {len(column_reactions)} column, {len(wall_reactions)} wall reactions.\n")

    except Exception as e:
        _send_progress(f"  ERROR during reaction model processing for '{reaction_level}': {type(e).__name__} - {str(e)}\n")
        raise 
    finally:
        if concept:
            _send_progress(f"  Extracting loads from reaction model complete.\n")

    # --- Apply Reactions to Target Model ---
    try:
        _send_progress(f"  Opening file for target model: {target_level}...\n")

        target_model = concept.open_file(target_path)
        _send_progress(f"  Target model opened. Preparing to apply loads...\n")
        
        cad_manager = target_model.cad_manager

        dead_load_target_name = "Other Dead Loading" 
        live_load_target_name = "Live (Unreducible) Loading" 
        
        dead_loading_target_layer = cad_manager.force_loading_layer(dead_load_target_name)
        live_loading_target_layer = cad_manager.force_loading_layer(live_load_target_name)
        target_layers_map = {
            dead_load_reaction_name: dead_loading_target_layer,
            live_load_reaction_name: live_loading_target_layer
        }

        _send_progress(f"  Applying column loads to '{target_level}'...\n")
        column_name_index = 0
        for col_reaction_data in column_reactions:
            target_layer = target_layers_map.get(col_reaction_data['Reaction Case'])
            if not target_layer:
                _send_progress(f"    WARNING: No target layer mapping for reaction case '{col_reaction_data['Reaction Case']}'. Skipping column load.\n")
                continue
            
            point_location = Point2D(col_reaction_data['Location X'], col_reaction_data['Location Y']) # Use actual Point2D
            point_load = target_layer.add_point_load(point_location)
            point_load.Fz = col_reaction_data['Fz'] 
            point_load.name = f"RCol{col_reaction_data['Column No.']}_{column_name_index}"
            column_name_index += 1
        _send_progress(f"    Applied {column_name_index} column loads.\n")

        _send_progress(f"  Applying wall loads to '{target_level}'...\n")
        wall_name_index = 0
        for wall_reaction_data in wall_reactions:
            target_layer = target_layers_map.get(wall_reaction_data['Reaction Case'])
            if not target_layer:
                _send_progress(f"    WARNING: No target layer mapping for reaction case '{wall_reaction_data['Reaction Case']}'. Skipping wall load.\n")
                continue

            wall_len = wall_reaction_data['Length']
            if wall_len == 0:
                _send_progress(f"    WARNING: Wall {wall_reaction_data['Wall No.']} has zero length. Skipping line load.\n")
                continue

            angle = wall_reaction_data['Angle']
            half_len_x = (wall_len / 2.0) * math.cos(angle)
            half_len_y = (wall_len / 2.0) * math.sin(angle)
            
            start_point = Point2D(wall_reaction_data['Centroid X'] - half_len_x, wall_reaction_data['Centroid Y'] - half_len_y) # Use actual Point2D
            end_point = Point2D(wall_reaction_data['Centroid X'] + half_len_x, wall_reaction_data['Centroid Y'] + half_len_y) # Use actual Point2D
            line_segment = LineSegment2D(start_point, end_point) # Use actual LineSegment2D
            
            line_load = target_layer.add_line_load(line_segment)
            line_load.Fz0 = wall_reaction_data['Fz'] / wall_len 
            line_load.Fz1 = wall_reaction_data['Fz'] / wall_len
            line_load.name = f"RWall{wall_reaction_data['Wall No.']}_{wall_name_index}"
            wall_name_index += 1
        _send_progress(f"    Applied {wall_name_index} wall loads.\n")

        _send_progress(f"  All loads applied. Saving target model '{target_path}'...\n")
        target_model.save_file(target_path)
        _send_progress(f"  Target model '{target_path}' saved.\n")

    except Exception as e:
        _send_progress(f"  ERROR during target model processing for '{target_path}': {type(e).__name__} - {str(e)}\n")
        raise 
    finally:
        _send_progress(f"ram_load_rundown: Successfully completed for '{reaction_path}' -> '{target_path}'.\n")


def delete_existing_loads(model_path:str, concept:Concept, progress_queue=None):
    '''
    Deletes all existing point and line loads that have a name from specified layers.
    Sends progress updates to the provided queue.

    Args:
        model_path (str): Path to the RAM Concept file to modify.
        progress_queue (queue.Queue, optional): Queue for sending progress messages.
                                                Defaults to None (prints to console).
    '''
    
    model_level = os.path.basename(model_path) # Extract names of level from path

    # Helper function to send progress messages
    def _send_progress(message):
        if progress_queue:
            progress_queue.put(message)
        else:
            print(message)

    _send_progress(f"Deleting existing loads for model '{model_path}'. Opening file... \n")
    try:
        model = concept.open_file(model_path)
        _send_progress(f"  Model opened.\n")

        dead_load_name = "Other Dead Loading"
        live_load_name = "Live (Unreducible) Loading"
        case_names_to_clear = [dead_load_name, live_load_name]

        cad_manager = model.cad_manager
        loadings = cad_manager.force_loading_layers # Get all force loading layers

        loads_deleted_count = 0
        found_layers_to_process = False

        for loading_layer in loadings: # Iterate through all available force loading layers
            if loading_layer.name in case_names_to_clear:
                found_layers_to_process = True
                _send_progress(f"    Processing layer '{loading_layer.name}' for load deletion...\n")
                
                # Process Point Loads
                if hasattr(loading_layer, 'point_loads'): # Check if the layer object has point_loads
                    for point_load in loading_layer.point_loads:
                        if hasattr(point_load, 'name') and point_load.name != '':
                            point_load.delete()
                            loads_deleted_count += 1

                # Process Line Loads
                # Create a list of line loads to delete to avoid issues with modifying list while iterating
                if hasattr(loading_layer, 'line_loads'): # Check if the layer object has line_loads
                    for line_load in loading_layer.line_loads:
                        if hasattr(line_load, 'name') and line_load.name != '':
                            line_load.delete()
                            loads_deleted_count += 1
                
        if not found_layers_to_process:
            _send_progress(f"  WARNING: None of the specified layers ({', '.join(case_names_to_clear)}) were found in the model's force loading layers.\n")

        if loads_deleted_count > 0:
            _send_progress(f"  Total {loads_deleted_count} named loads deleted. Saving model '{model_level}'...\n")
            model.save_file(model_path)
            _send_progress(f"  Model '{model_level}' saved after load deletion.\n")
        else:
            _send_progress(f"  No named loads were found to delete in the specified layers, or specified layers not found. Model not re-saved.\n")
        
        _send_progress(f"delete_exisitng_loads: Finished for model '{model_level}'.\n")

    except Exception as e:
        _send_progress(f"  ERROR during delete_exisitng_loads for '{model_level}': {type(e).__name__} - {str(e)}\n")
        raise # Re-raise to be caught by the UI's worker thread
    finally:
        if concept:
            _send_progress(f"  Load deletion complete.\n")


def create_column(model:Model, column_data:dict=None, load_case_names:list[str]=None, progress_queue=None):
    '''
    Creates a column in the target model based on column data and loading from ETABS.
    Sends progress updates to the provided queue.

    Args:
        model_path (str): Path to the RAM Concept file to modify.
        progress_queue (queue.Queue, optional): Queue for sending progress messages.
    '''
    # Helper function to send progress messages
    def _send_progress(message):
        if progress_queue:
            progress_queue.put(message)
        else:
            print(message)

    # Step 1 - Open StructureLayer (mesh input in UI)
    _send_progress('Step 1 - Open StructureLayer (mesh input in UI)\n')
    cad_manager = model.cad_manager
    structure_layer = cad_manager.structure_layer

    # Step 2 - Set default properties of column before placing
    _send_progress('Step 2 - Set default properties of column before placing\n')
    # for convenience we'll use the concrete mix we defined earlier
    concretes = model.concretes
    concrete_40 = concretes.concrete("40 MPa")
    
    # Now set remaining properties
    default_column = cad_manager.default_column
    default_column.angle = 90 + column_data['Angle']
    default_column.below_slab = False
    default_column.compressible = True
    default_column.concrete = concrete_40
    default_column.b = column_data['B']
    default_column.d = column_data['D']
    default_column.height = column_data['H'] / 1000
    default_column.i_factor = 1.0
    default_column.fixed_near = True
    default_column.fixed_far = True
    default_column.roller = False
    default_column.use_specified_LLR_parameters = False
    
    # Step 3 - Place column in location with Point2D
    _send_progress('Step 3 - Place column in location with Point2D\n')
    _send_progress (f'Location: {column_data["Point 2 X"]/1000}, {column_data["Point 2 Y"]/1000}\n')
    column_location = Point2D(column_data['Point 2 X']/1000, column_data['Point 2 Y']/1000)
    structure_layer.add_column(column_location)

    # Step 4 - Add loading to each layer, starting with Other Dead Loading
    _send_progress('Step 4 - Add loading to each layer, starting with Other Dead Loading\n')
    dead_loading_layer = cad_manager.force_loading_layer("Other Dead Loading")
    dead_load_case_name = load_case_names[0] # Use first load case name

    dead_point_load = dead_loading_layer.add_point_load(column_location) # Add point load to Point2D location
    dead_point_load.Fz = max(abs(load) for load in column_data[dead_load_case_name]['P']) / 1000 # Set the load value to the value from the load case data

    # Step 5 - Now add load to Unreducible Live Load layer
    _send_progress('Step 5 - Now add load to Unreducible Live Load layer\n')
    live_loading_layer = cad_manager.force_loading_layer("Live (Unreducible) Loading")
    live_load_case_name = load_case_names[1] # Use second load case name

    live_point_load = live_loading_layer.add_point_load(column_location) # Add point load to Point2D location
    live_point_load.Fz = max(abs(load) for load in column_data[live_load_case_name]['P']) / 1000 # Set the load value to the value from the load case data

    _send_progress('Column creation complete.\n')
    # Done