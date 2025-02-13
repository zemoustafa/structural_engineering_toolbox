"""

Module to create graphs and plots and then export them as readable files.

"""

import os
from decimal import Decimal
import pandas as pd
import tkinter as tk
from tkinter import filedialog

def dataframe_to_xlsx(dataframe: pd.DataFrame, default_filename: str = "output.xlsx"):
    """
    Exports a Pandas DataFrame to an XLSX file, allowing the user to specify the filename.

    Args:
        dataframe (pd.DataFrame): The DataFrame to export.
        default_filename (str): The default filename for the "Save As" dialog.
    """
    # Initialize the Tkinter root window (required for file dialogs)
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Open "Save As" dialog
    file_path = filedialog.asksaveasfilename(
        defaultextension=".xlsx",  # Ensure the extension is added
        filetypes=[("Excel Files", "*.xlsx")],  # Filter for Excel files
        title="Save As",
        initialfile=default_filename  # Set a default filename
    )

    if not file_path:
        print("No file selected. Operation cancelled.")
        return

    try:
        # Use ExcelWriter to save the DataFrame to an Excel file
        with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
            # Write the DataFrame to the Excel file
            dataframe.to_excel(writer, sheet_name="Sheet1", index=False)

            # Access the workbook and worksheet objects for formatting
            workbook = writer.book
            worksheet = writer.sheets["Sheet1"]

            # Define a format for the header row
            header_format = workbook.add_format({
                'bold': True,
                'text_wrap': False,
                'valign': 'center',
                'fg_color': '#0E70A6',
                'border': 1
            })

            # Apply the header format to the first row
            for col_num, value in enumerate(dataframe.columns):
                worksheet.write(0, col_num, value, header_format)

            # Autofit the columns
            worksheet.autofit()

        print(f"Excel file saved at: {file_path}")

    except Exception as e:
        print(f"An error occurred during export: {e}")
        return


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


    # # Plot the interaction curves and applied loads for x-x axis
    # fig_x, ax_x = plt.subplots()
    # ax_x.plot(f_m_x, f_n_x, label='Interaction Curve - About XX')
    # ax_x.plot(m_star_x, loading.n_star_compression, 'ro', label='Applied Load (Compression)')
    # ax_x.plot([0, m_star_x], [0, loading.n_star_compression], 'b-', label='Load Line (Compression)')
    # if loading.n_star_tension != 0:
    #     ax_x.plot(m_star_x, loading.n_star_tension, 'go', label='Applied Load (Tension)')
    #     ax_x.plot([0, m_star_x], [0, loading.n_star_tension], 'g-', label='Load Line (Tension)')
    # if point_in_diagram_x_comp:
    #     ax_x.plot([m_star_x, phi_m_x_comp], [loading.n_star_compression, phi_n_x_comp], 'b--', label='Extension Line (Compression)')
    # if loading.n_star_tension != 0 and point_in_diagram_x_tens:
    #     ax_x.plot([m_star_x, phi_m_x_tens], [loading.n_star_tension, phi_n_x_tens], 'g--', label='Extension Line (Tension)')
    # ax_x.set_ylabel('Axial Load (kN)')
    # ax_x.set_xlabel('Moment (kNm)')
    # ax_x.set_title('Moment Interaction Diagram - About XX')
    # ax_x.legend()
    # ax_x.grid(True)

    # # Plot the interaction curves and applied loads for y-y axis
    # if check_both_axes:
    #     fig_y, ax_y = plt.subplots()
    #     ax_y.plot(f_m_y, f_n_y, label='Interaction Curve - About YY')
    #     ax_y.plot(m_star_y, loading.n_star_compression, 'ro', label='Applied Load (Compression)')
    #     ax_y.plot([0, m_star_y], [0, loading.n_star_compression], 'b-', label='Load Line (Compression)')
    #     if loading.n_star_tension != 0:
    #         ax_y.plot(m_star_y, loading.n_star_tension, 'go', label='Applied Load (Tension)')
    #         ax_y.plot([0, m_star_y], [0, loading.n_star_tension], 'g-', label='Load Line (Tension)')
    #     if point_in_diagram_y_comp:
    #         ax_y.plot([m_star_y, phi_m_y_comp], [loading.n_star_compression, phi_n_y_comp], 'b--', label='Extension Line (Compression)')
    #     if loading.n_star_tension != 0 and point_in_diagram_y_tens:
    #         ax_y.plot([m_star_y, phi_m_y_tens], [loading.n_star_tension, phi_n_y_tens], 'g--', label='Extension Line (Tension)')
    #     ax_y.set_ylabel('Axial Load (kN)')
    #     ax_y.set_xlabel('Moment (kNm)')
    #     ax_y.set_title('Moment Interaction Diagram - About YY')
    #     ax_y.legend()
    #     ax_y.grid(True)
    # else:
    #     fig_y = None

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
