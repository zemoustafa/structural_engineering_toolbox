"""
Module to visualise data from ETABS model using pandas dataframes
and plotly.graph_objects. 

"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go


'''
DATAFRAMES

'''
def piers_as_walls_dataframe(walls:list[dict]) -> pd.DataFrame:
    ''' Create dataframe containing wall design.

    :param walls: list of dicts containing designed walls. must have run full_wall_design() function.
    :type walls: list[dicts]

    :return styled_df: styled dataframe containing wall design 
    :type df: pd.Dataframe
    '''
    selected_keys = [
            'Pier Name',
            'Story Name',
            'Story Height',
            'Thickness Bot',
            'Width Bot',
            'fc',
            'G+0.3Q (MPa)',
            'G+0.3Q+RS (C)(MPa)',
            'G+0.3Q-RS (T)(MPa)',
            'Axial Load Ratio',
            'Slenderness Ratio',
            'Rho crit.',
            'Rho typ.',
            'db Vert',
            's Vert',
            "0.15f'c",
            "0.2f'c",
            "0.585f'c",
            'BE Width',
            'Lig Dia',
            'Lig Cts',
            'EQ Shear',
            'Vuc',
            'Vus',
            'phiVu',
            'db Horiz',
            's Horiz'
            ]
    df = pd.DataFrame([{key: wall[key] for key in selected_keys} for wall in walls])

    return df

def pier_forces(pier_forces:list[dict]) -> pd.DataFrame:
    '''
    
    '''
    

    
    pass

'''
FIGURES

'''

def plot_loading_plans(floor_objs:list[dict], wall_objs:list[dict]) -> dict:
    """ Plot

    :param floor_objs: floor area objects from etabs model
    :type floor_objs: list[dict]
    :return figures_by_story: dict where key is story name and value is plotly figure
    :type figures_by_story: dict
    """
    figures_by_story = {}

    # Iterate through the floors data
    for floor in floor_objs:
        story_label = floor['Story']

        # If the story label is not in the dictionary, create a new figure for it
        if story_label not in figures_by_story:
            figures_by_story[story_label] = go.Figure()

        # Add a trace for each floor's data to the corresponding figure
        figures_by_story[story_label].add_trace(
            go.Scatter(
                x=floor['Point X'],
                y=floor['Point Y'],
                mode='lines',
                fill='tozeroy',
                name=f"{floor['Property']}-{floor['Diaphragm']}"
            )
        )

    for wall in wall_objs:
        story_label = wall['Story']
        # Add a trace for each floor's data to the corresponding figure
        figures_by_story[story_label].add_trace(
            go.Scatter(
                x=wall['X'],
                y=wall['Y'],
                mode='lines',
                marker={
                    'color': 'black', 
                    'line': {'width': 6}
                    }
                )
            )

    # Show or save the figures
    for story_label, figure in figures_by_story.items():
        figure.update_layout(title=f"Story: {story_label}")
        figure.update_yaxes(scaleanchor = "x", scaleratio = 1)

    return figures_by_story
