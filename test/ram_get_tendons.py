import os
import sys
import math

sys.path.append(r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\python")
import ram_concept


from ram_concept.concept import Concept
from ram_concept.cad_manager import CadManager
from ram_concept.enums import GeneratedBy
from ram_concept.line_segment_2D import LineSegment2D
from ram_concept.tendon_layer import TendonLayer
from ram_concept.tendon_segment import TendonSegment
from ram_concept.model import Model
from ram_concept.point_2D import Point2D

# link to the model path
model_path = r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\2nd floor slab design - IFC.cpt"

# pick tendon direction
tendon_direction = 'Latitude'

def __get_tendons_and_jacks(model_path, tendon_direction):
    """
    Private method.
    Gets tendon segments, nodes and jacks from RAM model.
    
    Args:
        model_path (str): Path to existing RAM model. Must be 100% working and checked before use.
        tendon_direction(str): User to select tendon direction to grab tendons from desired layer.
    
    Returns:
        segments (list): List of TendonSegment objects.
        jacks (list): List of of Jack objects.
    """

    # create instance of Concept class
    concept = Concept.start_concept(headless=True)

    try:
        # create instance of Concept class
        concept = Concept.start_concept(headless=True)

        # Open the file in the RAM Concept process and create a Model representing it.
        model = concept.open_file(model_path)
    except:
        print("Invalid model path or RAM Concept not installed.")

    # access layers in the models
    cad_manager = model.cad_manager

    # get all the tendon_layers
    tendon_layers = cad_manager.tendon_layers

    # select layer based on desired direction
    if tendon_direction == 'Latitude':
        selected_tendon_layer = tendon_layers[2] # latitude manual tendons
    elif tendon_direction == 'Longitude':
        selected_tendon_layer = tendon_layers[3]  # longitude manual tendons
    else:
        print("Error: Invalid tendon direction.")

    # grab segments and jacks from model
    segments = selected_tendon_layer.tendon_segments
    jacks = selected_tendon_layer.jacks

    return segments, jacks


def __chain_tendon_segments(all_segments):
    """
    Private method.
    Chains together all tendon segments into lists of connected segments.

    Args:
        all_segments (List[TendonSegment]): A list of all tendon segments in the model.

    Returns:
        List[List[TendonSegment]]: A list of chains, where each chain is a list of connected tendon segments.
    """
    visited_segments = set()  # To keep track of segments we've already processed
    chains = []  # List of chains

    def follow_chain(segment, chain):
        """Recursively follow a chain of tendon segments."""
        if segment.uid in visited_segments:
            return
        visited_segments.add(segment.uid)
        chain.append(segment)

        # Explore the next connected segments
        for node in [segment.node_1, segment.node_2]:
            connected_segments = node.connected_tendon_segments_except(segment)
            for next_segment in connected_segments:
                follow_chain(next_segment, chain)

    # Iterate over all segments and chain them together
    for segment in all_segments:
        if segment.uid not in visited_segments:
            chain = []
            follow_chain(segment, chain)
            chains.append(chain)

    return chains

def __convert_to_dicts(segment_chains, jacks, tendon_direction):
    """
    Private method.
    Converts the TendonSegment objects to dicts for future use and includes jack info.
    
    """
    # iterate through each segment
    chains_as_dicts = []
    current_chain = []

    for segment_chain in segment_chains:
        for segment in segment_chain:
            node_1 = segment.node_1
            node_2 = segment.node_2
            if tendon_direction == 'Latitude':
                start_node = node_1 if node_1.location.x < node_2.location.x else node_2
                end_node = segment.other_node(start_node)
                length = end_node.location.x - start_node.location.x
            elif tendon_direction == 'Longitude':
                start_node = node_1 if node_1.location.y < node_2.location.y else node_2
                end_node = segment.other_node(start_node)
                length = end_node.location.y - start_node.location.y
            current_chain.append({
                'UID': segment.uid,
                'Strand Number': segment.strand_count,
                'Start X': start_node.location.x,
                'Start Y': start_node.location.y,
                'End X': end_node.location.x,
                'End Y': end_node.location.y,
                'Start Elevation': int(round(segment.elevation_value_1, 0)),
                'End Elevation': int(round(segment.elevation_value_2, 0)),
                'Start Reference': segment.elevation_reference_1,
                'End Reference': segment.elevation_reference_2,
                'Intermediate Elevations': segment.elevations_along_segment(fractional_lengths = [i / length for i in range(0, int(length) + 1)] if length > 0 else [0]),
                'Type': 'Internal' # Internal is default value
            })
        chains_as_dicts.append(current_chain)
        current_chain = []
    
    
    checked_jacks = set() # keep track of which jacks have been used
    # add jacks and whether segment is internal or end to dicts
    for segment_chain in chains_as_dicts:
        start_segment = segment_chain[0] # grab first segment in chain
        end_segment = segment_chain[-1] # grab last segment in chain
        chain_start_point = (start_segment['Start X'], start_segment['Start Y'])
        chain_end_point = (end_segment['End X'], end_segment['End Y'])
        for jack in jacks:
            if jack.uid not in checked_jacks:
                jack_location = (jack.location.x, jack.location.y)
                if jack_location == chain_start_point:
                    start_segment['Type'] = 'Start ' + 'Live End' # TO UPDATE IN FUTURE
                    end_segment['Type'] = 'End ' + 'Dead End'
                    checked_jacks.add(jack.uid) 
                    break
                elif jack_location == chain_end_point:
                    end_segment['Type'] = 'End ' + 'Live End' # TO UPDATE IN FUTURE
                    start_segment['Type'] = 'Start ' + 'Dead End'
                    checked_jacks.add(jack.uid)   
                    break
        
    return chains_as_dicts

def get_ram_tendons(model_path, tendon_direction):
    segments, jacks = __get_tendons_and_jacks(model_path, tendon_direction)
    tendon_chains = __chain_tendon_segments(segments)
    tendon_chains_as_dicts = __convert_to_dicts(tendon_chains, jacks, tendon_direction)
    return tendon_chains_as_dicts

tendons = get_ram_tendons(model_path, tendon_direction)
print(tendons)


