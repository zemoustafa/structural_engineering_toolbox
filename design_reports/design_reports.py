"""

Module to create graphs and plots and then export them as readable files.

"""

import os
from decimal import Decimal
import pandas as pd
import tkinter as tk
from tkinter import filedialog

def dataframe_to_xlsx(dataframe: pd.DataFrame, title: str):
    # Open file dialog to select folder
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    folder_path = filedialog.askdirectory(title="Select Folder to Save Excel File")

    if not folder_path:  # If user cancels, return without saving
        print("No folder selected. Operation cancelled.")
        return

    # Define the full file path
    file_path = os.path.join(folder_path, f"{title}.xlsx")

    # Get number of rows of table
    df_shape = dataframe.shape
    num_rows = str(df_shape[0] + 1)

    # Create writer object
    writer = pd.ExcelWriter(file_path, engine='xlsxwriter')
    workbook = writer.book

    # Create worksheet from dataframe
    dataframe.to_excel(writer, sheet_name=title, index=False)
    worksheet1 = writer.sheets[title]

    # Create format for column titles
    header_format = workbook.add_format({
        'bold': True,
        'text_wrap': False,
        'valign': 'center',
        'fg_color': '#0E70A6',
        'border': 1
    })

    # Add format to columns in worksheet
    for col_num, value in enumerate(dataframe.columns.values):
        worksheet1.write(0, col_num, value, header_format)  # Adjusted to match column index

    # Autofit column widths
    worksheet1.autofit()

    # Export xlsx file
    writer.close()

    print(f"Excel file saved at: {file_path}")


'''
DATAFRAMES

'''
# def piers_as_walls_dataframe(walls:list[dict]) -> pd.DataFrame:
#     ''' Create dataframe containing wall design.

#     :param walls: list of dicts containing designed walls. must have run full_wall_design() function.
#     :type walls: list[dicts]

#     :return styled_df: styled dataframe containing wall design 
#     :type df: pd.Dataframe
#     '''
#     selected_keys = [
#             'Pier Name',
#             'Story Name',
#             'Story Height',
#             'Thickness Bot',
#             'Width Bot',
#             'fc',
#             'G+0.3Q (MPa)',
#             'G+0.3Q+RS (C)(MPa)',
#             'G+0.3Q-RS (T)(MPa)',
#             'Axial Load Ratio',
#             'Slenderness Ratio',
#             'Rho crit.',
#             'Rho typ.',
#             'db Vert',
#             's Vert',
#             "0.15f'c",
#             "0.2f'c",
#             "0.585f'c",
#             'BE Width',
#             'Lig Dia',
#             'Lig Cts',
#             'EQ Shear',
#             'Vuc',
#             'Vus',
#             'phiVu',
#             'db Horiz',
#             's Horiz'
#             ]
#     df = pd.DataFrame([{key: wall[key] for key in selected_keys} for wall in walls])

#     return df



'''
FIGURES

'''

# def plot_loading_plans(floor_objs:list[dict], wall_objs:list[dict]) -> dict:
    # """ Plot

    # :param floor_objs: floor area objects from etabs model
    # :type floor_objs: list[dict]
    # :return figures_by_story: dict where key is story name and value is plotly figure
    # :type figures_by_story: dict
    # """
    # figures_by_story = {}

    # # Iterate through the floors data
    # for floor in floor_objs:
    #     story_label = floor['Story']

    #     # If the story label is not in the dictionary, create a new figure for it
    #     if story_label not in figures_by_story:
    #         figures_by_story[story_label] = go.Figure()

    #     # Add a trace for each floor's data to the corresponding figure
    #     figures_by_story[story_label].add_trace(
    #         go.Scatter(
    #             x=floor['Point X'],
    #             y=floor['Point Y'],
    #             mode='lines',
    #             fill='tozeroy',
    #             name=f"{floor['Property']}-{floor['Diaphragm']}"
    #         )
    #     )

    # for wall in wall_objs:
    #     story_label = wall['Story']
    #     # Add a trace for each floor's data to the corresponding figure
    #     figures_by_story[story_label].add_trace(
    #         go.Scatter(
    #             x=wall['X'],
    #             y=wall['Y'],
    #             mode='lines',
    #             marker={
    #                 'color': 'black', 
    #                 'line': {'width': 6}
    #                 }
    #             )
    #         )

    # # Show or save the figures
    # for story_label, figure in figures_by_story.items():
    #     figure.update_layout(title=f"Story: {story_label}")
    #     figure.update_yaxes(scaleanchor = "x", scaleratio = 1)

    # return figures_by_story
