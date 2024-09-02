import comtypes.client

"""
Creates instance of sap model object to access data in etabs model

"""
helper = comtypes.client.CreateObject("ETABSv1.Helper")
helper = helper.QueryInterface(comtypes.gen.ETABSv1.cHelper)

try:
    etabs_object = helper.GetObject("CSI.ETABS.API.ETABSObject")
except (OSError, comtypes.COMError):
    raise Exception(
        "No running instance of the program found or failed to attach."
    )

sap_model = etabs_object.SapModel

story_data = {
    'GROUND': 1000,
    'LEVEL 1': 4000,
    'ROOF': 7000
}

first_pair = next(iter(story_data.items())) # grab first key value pair
base_elev = first_pair[1]
num_stories = len(story_data)
story_names = list(story_data.keys())
story_elevations = list(story_data.values())
story_heights = [story_elevations[0]] + [story_elevations[i] - story_elevations[i-1] for i in range(1, len(story_elevations))]
is_master = [False] * num_stories
similar_story = ["None"] * num_stories
splice_above = [False] * num_stories
splice_h = [0] * num_stories

sap_model.Story.SetStories_2(
    base_elev,
    num_stories-1,
    story_names[1:], 
    story_heights[1:], 
    is_master[1:], 
    similar_story[1:], 
    splice_above[1:], 
    splice_h[1:]
    )