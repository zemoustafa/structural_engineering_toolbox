import comtypes.client
import sys
sys.path.append('C:\_Github\structural_engineering_toolbox')
from etabs_tools import etabs_api, etabs_design 
from as3600 import vertical_structure  
import pandas as pd

etabs_api = etabs_api.etabs_api()
sap_model = etabs_api.sap_model

eq_env_1 = '(88) RS ULS ENV' # earthquake envelope unfactored EQ for moment design
eq_env_2 = '(88) RS ULS ENV SHEAR' # earthquake envelope factored EQ for amplified shear design
wind_env = '(88) WIND ULS ENV' # wind envelope
load_cases = [eq_env_1, eq_env_2, wind_env]

piers = etabs_api.get_piers(load_cases=load_cases)
designed_piers_df = etabs_design.design_all_piers(
    piers=piers,
    eq_env_1=eq_env_1,
    eq_env_2=eq_env_2,
    wind_env=wind_env,
    vertical_spacing=200,
    horizontal_spacing=200,
    design_both_axes=False
)

pass
# """
# Creates instance of sap model object to access data in etabs model

# """
# helper = comtypes.client.CreateObject("ETABSv1.Helper")
# helper = helper.QueryInterface(comtypes.gen.ETABSv1.cHelper)

# try:
#     etabs_object = helper.GetObject("CSI.ETABS.API.ETABSObject")
# except (OSError, comtypes.COMError):
#     raise Exception(
#         "No running instance of the program found or failed to attach."
#     )

# sap_model = etabs_object.SapModel

# story_data = {
#     'GROUND': 1000,
#     'LEVEL 1': 4000,
#     'ROOF': 7000
# }

# first_pair = next(iter(story_data.items())) # grab first key value pair
# base_elev = first_pair[1]
# num_stories = len(story_data)
# story_names = list(story_data.keys())
# story_elevations = list(story_data.values())
# story_heights = [story_elevations[0]] + [story_elevations[i] - story_elevations[i-1] for i in range(1, len(story_elevations))]
# is_master = [False] * num_stories
# similar_story = ["None"] * num_stories
# splice_above = [False] * num_stories
# splice_h = [0] * num_stories

# sap_model.Story.SetStories_2(
#     base_elev,
#     num_stories-1,
#     story_names[1:], 
#     story_heights[1:], 
#     is_master[1:], 
#     similar_story[1:], 
#     splice_above[1:], 
#     splice_h[1:]
#     )