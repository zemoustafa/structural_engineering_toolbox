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

# step 4 - serever transport to receive model data
speckle_client.server_transport()

# step 5 - get revit floors from model
revit_floors = speckle_client.get_revit_floors()


pass


# """
# STEP BY STEP HOW TO DRAW MODEL IN ETABS

# """

# # step 1 - create instance of sap model object from etabs api class
# sap_model = etabs_api.etabs_api()

# # step 2 - input unique levels in etabs model
# ret = sap_model.set_stories(unique_levels)

# # step 3 - draw floors in etabs model
# revit_floors = speckle_client.get_revit_floors()
    
# for floor in revit_floors:
#     sap_model.sap_model.AreaObj.AddByCoord(
#     NumberPoints = len(floor['X']),
#     X = floor['X'],
#     Y = floor['Y'],
#     Z = floor['Z'],
#     Name = floor['Name']
# )
 
# step 4 - draw columns
revit_columns = speckle_client.get_revit_columns()
pass