import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from speckle_tools import speckle_client

# inputs
private_token = "a17392654d309903769de972c15296ec3dc4b8a7c5"
model_url = "https://app.speckle.systems/projects/a5f5d0cbdc/models/8084b76bba@1e1e6adbb9"

"""
STEP BY STEP HOW TO GRAB UNIQUE LEVELS

"""

# step 1 - create speckle_client class
speckle_client = speckle_client.speckle_client(model_url, private_token)

gql_client = speckle_client.gql_client()
unique_levels = speckle_client.get_unique_levels()

print(unique_levels)