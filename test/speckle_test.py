import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from speckle_tools import speckle_client
from etabs_tools import etabs_api

# inputs
private_token = "a17392654d309903769de972c15296ec3dc4b8a7c5"
model_url = "https://app.speckle.systems/projects/a5f5d0cbdc/models/8084b76bba@1e1e6adbb9"

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
# print(ret)
