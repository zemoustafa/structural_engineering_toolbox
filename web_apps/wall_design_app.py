import numpy as np
import pandas as pd
import streamlit as st
import xlsxwriter
import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from modules import as3600_wall_design, etabs_api, design_reports

# -------------------------------------------------------------------------------------------------------------------
# SET UP PAGE
st.set_page_config(
    page_title="ETABS Wall Design Tool",
    page_icon=":book:",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# create header
st.title(body="ETABS Wall Design Tool")
# -------------------------------------------------------------------------------------------------------------------
# SET UP CONTAINERS
# create container that holds all items
setup_col, results_col = st.columns([1, 4])
# -------------------------------------------------------------------------------------------------------------------
# SET UP BUTTON SESSION STATES

if "load_cases" not in st.session_state:
    st.session_state.load_cases = []

if "selected_load_cases" not in st.session_state:
    st.session_state.selected_load_cases = []
    
if "v_cts" not in st.session_state:
    st.session_state.v_cts = []

if "h_cts" not in st.session_state:
    st.session_state.h_cts = []

if "phz_levels" not in st.session_state:
    st.session_state.phz_levels = []

if "story_data" not in st.session_state:
    st.session_state.story_data = []

if "story_names" not in st.session_state:
    st.session_state.story_names = []

if 'walls_dicts' not in st.session_state:
    st.session_state.walls_dicts = []

if "walls_df" not in st.session_state:
    st.session_state.walls_df = pd.DataFrame()

if 'optimised_walls_dicts' not in st.session_state:
    st.session_state.optimised_walls_dicts = []

if "optimised_walls_df" not in st.session_state:
    st.session_state.optimised_walls_df = pd.DataFrame()
    
if "columns_df" not in st.session_state:
    st.session_state.columns_df = pd.DataFrame()

if 'index' not in st.session_state:
    st.session_state.index = 0

if 'pier_forces' not in st.session_state:
    st.session_state.pier_forces = []
    
# -------------------------------------------------------------------------------------------------------------------
# SET UP MAIN CONTAINER

# create instance of class object and grab SapModel
etabsAPI = etabs_api.etabs_api()

# SET UP COLUMN 1 - CONTAINS SET UP
with setup_col:
    # add buttons
    get_load_cases_button = st.button(
        label="Get Load Cases", 
        help="Make sure the model is open",
        use_container_width=True
        )

    # once button clicked, make visible the selectboxes to select case names
    if get_load_cases_button:
        if etabsAPI:
            st.session_state.load_cases = etabsAPI.get_load_cases()
            st.session_state.story_data = etabsAPI.get_story_data()
            st.session_state.story_names = etabsAPI.story_names
            
        else:
            st.error("No ETABS model detected. Please open an ETABS model.")

    form = st.form('input_form')
    with form:
        form.g = st.selectbox(label="Dead Load", options=st.session_state.load_cases)
        form.q = st.selectbox(label="Live Load", options=st.session_state.load_cases)
        form.rs = st.selectbox(label="Response Spectrum", options=st.session_state.load_cases)
        form.wall_type = st.radio(label="Wall type", options=["In-situ", "Precast"])
        st.write("<style>div.row-widget.stRadio>div{flex-direction:row;}</style",unsafe_allow_html=True,)
        form.mu_sp = st.selectbox(label="Ductility factor", options=[1.3, 2.6, 4.5], index=1)

        form.start_phz = st.selectbox(label="Start of plastic hinge zone", options=st.session_state.story_names)
        form.num_phz_levels = st.number_input(label="Number of levels within PHZ", step=1)
        
        form.h_cts_max = st.selectbox(label="Max horiz. bar spacing", options=[150, 200, 250, 300, 400], index=3)
        form.h_cts_min = st.selectbox(label="Min horiz. bar spacing", options=[150, 200, 250, 300, 400], index=1)
        
        form.v_cts_max = st.selectbox(label="Max vert. bar spacing", options=[150, 200, 250, 300, 400], index=3)
        form.v_cts_min = st.selectbox(label="Min vert. bar spacing", options=[150, 200, 250, 300, 400], index=1)
        
        # form submit button
        submit_button = form.form_submit_button('Get pier forces', use_container_width=True)

    st.session_state.selected_load_cases = [form.g, form.q, form.rs]
    st.session_state.h_cts = [form.h_cts_max, form.h_cts_min]
    st.session_state.v_cts = [form.v_cts_max, form.v_cts_min]
    
    if submit_button:
        del st.session_state.pier_forces[:]
        for case in st.session_state.selected_load_cases:        
            st.session_state.pier_forces.append(etabsAPI.get_pier_forces(case))
        st.session_state.phz_levels = as3600_wall_design.get_phz_stories(form.start_phz, form.num_phz_levels, st.session_state.story_names)
        st.session_state.walls_dicts = etabsAPI.get_piers(st.session_state.selected_load_cases)  
        piers_as_walls = as3600_wall_design.full_wall_design(
                    st.session_state.walls_dicts, 
                    st.session_state.selected_load_cases,
                    st.session_state.v_cts,
                    st.session_state.h_cts,
                    st.session_state.story_names,
                    st.session_state.phz_levels,
                    form.wall_type,
                    form.mu_sp
                    )
        st.session_state.walls_df = design_reports.piers_as_walls_dataframe(piers_as_walls)

    export_xlsx_button = st.button(label='Export XLSX', use_container_width=True)

    # EXPORT XLSX - ADD TO FUNCTION LATER 
    if export_xlsx_button:
        # get number of rows of table
        df_shape = st.session_state.walls_df.shape
        num_rows = str(df_shape[0]+1)

        # create wrtier object
        writer = pd.ExcelWriter('pandas_simple.xlsx', engine='xlsxwriter')
        workbook  = writer.book

        # create worksheet1 from dataframe
        st.session_state.walls_df.to_excel(writer, sheet_name="Wall Design")
        worksheet1 = writer.sheets['Wall Design']        
        
        # create format for column titles
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': False,
            'valign': 'center',
            'fg_color': '#0E70A6',
            'border': 1
            })
        
        # add format to columns in worksheet11
        for col_num, value in enumerate(st.session_state.walls_df.columns.values):
            worksheet1.write(0, col_num + 1, value, header_format)

        # add conditional format to G+0.3Q+RS column
        compression_format = worksheet1.conditional_format('I1:I'+num_rows, options={
            'type': '3_color_scale',
            'min_color': '#8EE80E',
            'mid_color': '#FFF400',
            'max_color': '#D61600'
            }
            )
        
        # add conditional format to G+0.3Q-RS column
        tension_format = worksheet1.conditional_format('J1:J'+num_rows, options={
            'type': '3_color_scale',
            'min_color': '#D61600',
            'mid_color': '#FFF400',
            'max_color': '#8EE80E'
            }
            )
        
        # create and add formats for boundary elements
        format1 = workbook.add_format({"bg_color": "#FF4000"}) # red
        format2 = workbook.add_format({"bg_color": "#8EE80E"}) # green
        be_format_pass = worksheet1.conditional_format('Q2:S'+num_rows, options={
            'type': 'cell',
            'criteria': '>=',
            'value': '$I$2:$I$1048576',
            'format': format2
        })
        be_format_fail = worksheet1.conditional_format('Q2:S'+num_rows, options={
            'type': 'cell',
            'criteria': '<',
            'value': '$I$2:$I$1048576',
            'format': format1
        })

        # autofit column widths
        worksheet1.autofit()

        # export xlsx file
        writer.close()


# -------------------------------------------------------------------------------------------------------------------
# SET UP COLUMN 2 - CONTAINS TABS WITH RESULTS
with results_col:
    # set up results tabs and start page tabs
    wall_tab, column_tab = st.tabs(["Piers as walls", "Piers as columns"])

    
       
    # -------------------------------------------------------------------------------------------------------------------
    with wall_tab:
        st.dataframe(st.session_state.walls_df, use_container_width=True, height=1000)
        # wall_options_col, wall_df_col = st.columns([1, 4])

        # with wall_options_col:
        #     st.subheader("Optimise walls")
        #     st.text_input(label="Max. % of f'c", value='0.2', key='max_fc')
        #     optimise_button = st.button('Optimise walls')
        #     show_optimised_walls_checkbox = st.checkbox(label='Show optimised walls', key='optimised_walls_checkbox')
        #     if optimise_button:
        #         stress_limit = float(st.session_state.max_fc)
        #         st.session_state.optimised_walls_dicts = as3600_wall_design.optimise_walls(
        #             st.session_state.walls_dicts, 
        #             stress_limit, 
        #             st.session_state.selected_load_cases, 
        #             st.session_state.mu_sp
        #             )
        #         st.session_state.optimised_walls_redesigned = as3600_wall_design.full_as3600_wall_design(
        #             st.session_state.optimised_walls_dicts, 
        #             st.session_state.selected_load_cases,
        #             st.session_state.v_cts,
        #             st.session_state.h_cts,
        #             st.session_state.story_names,
        #             st.session_state.phz_levels,
        #             st.session_state.wall_type,
        #             st.session_state.mu_sp
        #             )
        #         st.session_state.optimised_walls_df = as3600_wall_design.piers_as_walls_daframe(st.session_state.optimised_walls_redesigned)

        #     update_etabs_button = st.button(label='Update ETABS Model')

        # with wall_df_col:
        #     if show_optimised_walls_checkbox:
        #         st.dataframe(st.session_state.optimised_walls_df, use_container_width=True, height=1000)
        #     else:
        #         st.dataframe(st.session_state.walls_df, use_container_width=True, height=1000)

        
    # -------------------------------------------------------------------------------------------------------------------
    with column_tab:
        st.dataframe(st.session_state.columns_df, use_container_width=True, height=1000)
        

        


# SIDEBAR TO DEBUG APP
with st.sidebar:
    st.write("Selected load cases")
    st.session_state.selected_load_cases
    st.write(len(st.session_state.selected_load_cases))

    st.write("Pier forces list")
    st.session_state.pier_forces

    st.write('Story names')
    st.session_state.story_names

    st.write('Start PHZ')
    st.write(form.start_phz)

    st.write('Stories within PHZ')
    st.session_state.phz_levels

    