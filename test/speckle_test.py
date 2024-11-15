import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from speckle_tools import speckle_client
from etabs_tools import etabs_api

# inputs
private_token = "50278a3a662962a9426c798ec642ea164abb4f0914"
model_url = "https://app.speckle.systems/projects/9439578ade/models/166f646beb@fa5bf97c18"

"""
STEP BY STEP HOW TO GRAB UNIQUE LEVELS

"""

# step 1 - create speckle_client class
speckle_client = speckle_client.speckle_client(model_url, private_token)

# step 2 - create instance of gql client
gql_client = speckle_client.gql_client()

# step 3 - grab dict of unique levels from speckle model
unique_levels = speckle_client.get_unique_levels()
story_elevations = list(unique_levels.values()) # create list of just elevations
story_names = list(unique_levels.keys()) # create list of just story names

# step 4 - serever transport to receive model data
speckle_client.server_transport()

# step 5 - get revit floors from model
revit_floors = speckle_client.get_revit_floors()

# step 6 - get revit columns from model
revit_columns = speckle_client.get_revit_columns(unique_levels)

# step 7 - get revit walls from model
revit_walls = speckle_client.get_revit_walls(unique_levels)

"""
STEP BY STEP HOW TO DRAW MODEL IN ETABS

"""

# step 1 - create instance of sap model object from etabs api class
sap_model = etabs_api.etabs_api()

# step 2 - input unique levels in etabs model
ret = sap_model.set_stories(unique_levels)

# step 3 - draw floors in etabs model
for floor in revit_floors:
    ret = sap_model.sap_model.AreaObj.AddByCoord(
    NumberPoints = len(floor['X']),
    X = floor['X'],
    Y = floor['Y'],
    Z = floor['Z'],
    Name = floor['Name']
)
 
# step 4 - draw columns
for column in revit_columns:
    ret = sap_model.sap_model.FrameObj.AddByCoord(
        column['Bottom X'],
        column['Bottom Y'],
        column['Bottom Z'],
        column['Top X'],
        column['Top Y'],
        column['Top Z'],
        column['Name']
    )

# step 5 - draw walls
for wall in revit_walls:
    X_coord = [wall['Start X'], wall['Start X'], wall['End X'], wall['End X'],]
    Y_coord = [wall['Start Y'], wall['Start Y'], wall['End Y'], wall['End Y']]
    Z_coord = [wall['Bottom Z'], wall['Top Z'], wall['Top Z'], wall['Bottom Z']]
    ret = ret = sap_model.sap_model.AreaObj.AddByCoord(
        4,
        X_coord,
        Y_coord,
        Z_coord,
        wall['Name']
    )

pass






# iterate within range of start and end level
# in_range = False
# for index, level in enumerate(story_names):
#     if level == baseLevel:
#         in_range = True
#     if in_range:
#         # extract start and end coordinates of wall
#         startX = round(start.x, 0)
#         startY = round(start.y, 0)
#         startZ = story_elevations[index]
#         endX = round(end.x, 0)
#         endY = round(end.y, 0)
#         endZ = story_elevations[index + 1]

#         # DRAW WALL IN ETABS
#         X_coord = [startX, startX, endX, endX]
#         Y_coord = [startY, startY, endY, endY]
#         Z_coord = [startZ, endZ, endZ, startZ]

#         Name = " "
#         ret = SapModel.AreaObj.AddByCoord(4, X_coord, Y_coord, Z_coord, Name, "CW-200-C40")

#     if story_names[index + 1] == topLevel:
#         in_range = False
#         break

# DRAW COLUMNS THIS WAY IF MODELLED OVER SEVERAL LEVELS
# iterate within range of start and end level
# for column in revit_columns:
#     in_range = False
#     for index, level in enumerate(story_names):
#         if level == revit_columns:
#             in_range = True
#         if in_range:
#             # extract start and end coordinates of column
#             startX = round(column['Bottom X'], 0)
#             startY = round(column['Bottom Y'], 0)
#             startZ = story_elevations[index]
#             endX = round(column['Top X'], 0)
#             endY = round(column['Top Y'], 0)
#             endZ = story_elevations[index + 1]
#             # add column with ETABS API
#             ret = sap_model.SapModel.FrameObj.AddByCoord(startX, startY, startZ, endX, endY, endZ)
#         if story_names[index + 1] == column['Top Level']:
#             in_range = False
#             break