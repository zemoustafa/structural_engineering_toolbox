"""
Module to visualise data from ETABS model using pandas dataframes
and plotly.graph_objects. 

"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

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
