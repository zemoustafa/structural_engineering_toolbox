import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from speckle_tools import speckle_client
from etabs_tools import etabs_api

"""
Create instance of speckle_client class and grab unique levels

"""
private_token = "a17392654d309903769de972c15296ec3dc4b8a7c5"
model_url = "https://app.speckle.systems/projects/a5f5d0cbdc/models/8084b76bba@1e1e6adbb9"

speckle_client = speckle_client.speckle_client(model_url, private_token)
gql_client = speckle_client.gql_client()

unique_levels = speckle_client.get_unique_levels()

"""
Create instance of etabs_api class

"""
etabs_api = etabs_api.etabs_api()

"""
Set levels and elevations in etabs model from data grabbed from speckle via graphql

"""
set_stories = etabs_api.set_stories(story_data=unique_levels)

