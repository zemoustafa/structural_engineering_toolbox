import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox\etabs_tools")
import etabs_api
import pandas as pd

etabs_api = etabs_api.etabs_api()
 
sap_model = etabs_api.sap_model

def get_columns_on_story(story_name):
    # Step 1 - Grab all frame objects on a certain story
    frame_names_on_story = sap_model.FrameObj.GetNameListOnStory(story_name) # Story name will be an input
    frame_names = frame_names_on_story[1] # Frame object names

    # Step 2 - Filter out beams, only keep columns
    columns = []
    for frame_name in frame_names:
        # Get the two joints either side of the frame
        points = sap_model.FrameObj.GetPoints(frame_name)
        point_1 = points[0]
        point_2 = points[1]
        point_1_coords = sap_model.PointObj.GetCoordCartesian(point_1)
        point_2_coords = sap_model.PointObj.GetCoordCartesian(point_2)

        # If the two points have the same x and y coordinate, it's a column
        if point_1_coords[0] == point_2_coords[0] and point_1_coords[1] == point_2_coords[1]:
            column = {
                'Name': frame_name,
                'Point 1 X': point_1_coords[0],
                'Point 1 Y': point_1_coords[1],
                'Point 1 Z': point_1_coords[2],
                'Point 2 X': point_2_coords[0],
                'Point 2 Y': point_2_coords[1],
                'Point 2 Z': point_2_coords[2],
                'H': point_2_coords[2] - point_1_coords[2], # Height (mm)
                'Story': story_name,
            }
            columns.append(column)

    # Step 3 - Grab the properties of each column
    for column in columns:
        # Grab the section property name of the column
        column_section = sap_model.FrameObj.GetSection(column['Name'])

        # Try to get concrete column properties as a rectangular section
        column_section_name = column_section[0]
        column_section = sap_model.PropFrame.GetRectangle(column_section_name)

        # If the section is not rectangular, try as circle
        if column_section[-1] == 1:
            column_section = sap_model.PropFrame.GetCircle(column['Name'])

        # Add section properties to column dict
        column['fc'] = column_section[1] # Concrete strength (MPa)
        column['D'] = column_section[2] # Depth (mm)
        column['B'] = column_section[3] # Width (mm)


    # Step 4 - Add loading to column dicts
    dead_load = 'V:G'
    live_load = 'V:Q'
    load_case_names = [dead_load, live_load]

    # Input for this funciton will be a list of load case names

    for column in columns:
        # Iterate through each load case name
        for load_case_name in load_case_names:
            # Deselect all cases and combos for output
            sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
            
            # Try as load case first, otherwise try as combo
            output = sap_model.Results.Setup.SetCaseSelectedForOutput(load_case_name)
            if output == 0:
                continue
            else:
                output = sap_model.Results.Setup.SetComboSelectedForOutput(load_case_name)

            # Grab frame forces for current load case and current column
            frame_forces = sap_model.Results.FrameForce(column['Name'], 1)

            # Create dict of forces for current load case
            frame_force_current_load_case_dict = {
                'NumberResults': frame_forces[0],
                'Obj': frame_forces[1],
                'ObjSta': frame_forces[2],
                'Elm': frame_forces[3],
                'ElmSta': frame_forces[4],
                'LoadCase': frame_forces[5],
                'StepType': frame_forces[6],
                'StepNum': frame_forces[7],
                'P': frame_forces[8],
                'V2': frame_forces[9],
                'V3': frame_forces[10],
                'T': frame_forces[11],
                'M2': frame_forces[12],
                'M3': frame_forces[13],
            }
            
            # Add force dict to column dict
            column[load_case_name] = frame_force_current_load_case_dict

    return columns

pd.DataFrame(columns).to_csv('columns.csv')