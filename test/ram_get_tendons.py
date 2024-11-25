import os
import sys
import pandas as pd
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
model_path = r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\241119 23-064 Level 2 Flat Plate Option 1.6 DL 260 thk ZM.cpt"

# create instance of Concept class
concept = Concept.start_concept(headless=True)

# Open the file in the RAM Concept process and create a Model representing it.
model = concept.open_file(model_path)

# access layers in the models
cad_manager = model.cad_manager

# get all the tendon_layers
tendon_layers = cad_manager.tendon_layers

# separate into manual latitude and longitude layers
latitude_tendon_layer = tendon_layers[2]
longitude_tendon_layer = tendon_layers[3]

lat_segments = []

# grab a segment
for segment in latitude_tendon_layer.tendon_segments:
    lat_segments.append({   
        'Strand Number': segment.strand_count,
        'Start X': segment.location.start_point.x,
        'Start Y': segment.location.start_point.y,
        'End X': segment.location.end_point.x,
        'End Y': segment.location.end_point.y,
        'Start Elevation': segment.elevation_value_1,
        'End Elevation': segment.elevation_value_2,
        'Start Reference': segment.elevation_reference_1,
        'End Reference': segment.elevation_reference_2,
    })

concept.shut_down()

pass