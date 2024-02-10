import numpy as np
import pandas as pd
import streamlit as st
import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from modules import as3600_wall_design, etabs_api

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
# -------------------------------------------------------------------------------------------------------------------
# SET UP MAIN CONTAINER

# create instance of class object and grab SapModel
etabsAPI = etabs_api.etabs_api()

# SET UP COLUMN 1 - CONTAINS SET UP
with setup_col:
    # add buttons
    buttonCol1, buttonCol2 = st.columns([0.5, 0.5])
    with buttonCol1:
        get_load_cases_button = st.button(
            label="Get Load Cases", 
            help="Make sure the model is open"
            )
    with buttonCol2:
        get_pier_forces_button = st.button(
            "Get Pier Forces",
            help="Please make sure your ETABS model is correct, all piers labels have been assigned and the model has been run.",
        )

    # label - select load case names
    st.info("Select names of load cases from ETABS model.")

    # once button clicked, make visible the selectboxes to select case names
    if get_load_cases_button:
        if etabsAPI:
            st.session_state.load_cases = etabsAPI.get_load_cases()
            st.session_state.story_names = etabsAPI.get_story_names()
        else:
            st.error("No ETABS model detected. Please open an ETABS model.")

    g = st.selectbox(label="Dead Load", options=st.session_state.load_cases, key='g')
    q = st.selectbox(label="Live Load", options=st.session_state.load_cases, key='q')
    rs = st.selectbox(label="Response Spectrum", options=st.session_state.load_cases, key='rs')
    wx = st.selectbox(label="Wind X", options=st.session_state.load_cases, key='wx')
    wy = st.selectbox(label="Wind Y", options=st.session_state.load_cases, key='wy')
    st.session_state.selected_load_cases = [st.session_state.g, st.session_state.q, st.session_state.rs, st.session_state.wx, st.session_state.wy,]

    wall_type = st.radio(label="Wall type", options=["In-situ", "Precast"], key='wall_type')
    st.write("<style>div.row-widget.stRadio>div{flex-direction:row;}</style",unsafe_allow_html=True,)
    mu_sp = st.selectbox(label="Ductility factor", options=[1.3, 2.6, 4.5], index=1, key='mu_sp')

    start_phz = st.selectbox(label="Start of plastic hinge zone", options=st.session_state.story_names, key='start_phz')
    num_phz_levels = st.number_input(label="Number of levels within PHZ", step=1, key='num_phz_levels')
    
    h_cts_max = st.selectbox(label="Max horiz. bar spacing", options=[150, 200, 250, 300, 400], index=3)
    h_cts_min = st.selectbox(label="Min horiz. bar spacing", options=[150, 200, 250, 300, 400], index=1)
    st.session_state.h_cts = [h_cts_max, h_cts_min]

    v_cts_max = st.selectbox(label="Max vert. bar spacing", options=[150, 200, 250, 300, 400], index=3)
    v_cts_min = st.selectbox(label="Min vert. bar spacing", options=[150, 200, 250, 300, 400], index=1)
    st.session_state.v_cts = [v_cts_max, v_cts_min]

# -------------------------------------------------------------------------------------------------------------------
# SET UP COLUMN 2 - CONTAINS TABS WITH RESULTS
with results_col:
    # set up results tabs and start page tabs
    wall_tab, column_tab = st.tabs(["Piers as walls", "Piers as columns"])

    piers_as_columns = []
    
    if get_pier_forces_button:

        st.session_state.phz_levels = as3600_wall_design.get_phz_stories(st.session_state.start_phz, st.session_state.num_phz_levels, st.session_state.story_names)
        
        st.session_state.walls_dicts = etabsAPI.get_piers(st.session_state.selected_load_cases)
        
        piers_as_walls = as3600_wall_design.full_wall_design(
                    st.session_state.walls_dicts, 
                    st.session_state.selected_load_cases,
                    st.session_state.v_cts,
                    st.session_state.h_cts,
                    st.session_state.story_names,
                    st.session_state.phz_levels,
                    st.session_state.wall_type,
                    st.session_state.mu_sp
                    )
        
        st.session_state.walls_df = as3600_wall_design.piers_as_walls_dataframe(piers_as_walls)
       
    # -------------------------------------------------------------------------------------------------------------------
    with wall_tab:
        wall_options_col, wall_df_col = st.columns([1, 4])

        with wall_options_col:
            st.subheader("Optimise walls")
            st.text_input(label="Max. % of f'c", value='0.2', key='max_fc')
            optimise_button = st.button('Optimise walls')
            show_optimised_walls_checkbox = st.checkbox(label='Show optimised walls', key='optimised_walls_checkbox')
            if optimise_button:
                stress_limit = float(st.session_state.max_fc)
                st.session_state.optimised_walls_dicts = as3600_wall_design.optimise_walls(
                    st.session_state.walls_dicts, 
                    stress_limit, 
                    st.session_state.selected_load_cases, 
                    st.session_state.mu_sp
                    )
                st.session_state.optimised_walls_redesigned = as3600_wall_design.full_as3600_wall_design(
                    st.session_state.optimised_walls_dicts, 
                    st.session_state.selected_load_cases,
                    st.session_state.v_cts,
                    st.session_state.h_cts,
                    st.session_state.story_names,
                    st.session_state.phz_levels,
                    st.session_state.wall_type,
                    st.session_state.mu_sp
                    )
                st.session_state.optimised_walls_df = as3600_wall_design.piers_as_walls_daframe(st.session_state.optimised_walls_redesigned)

            update_etabs_button = st.button(label='Update ETABS Model')

        with wall_df_col:
            if show_optimised_walls_checkbox:
                st.dataframe(st.session_state.optimised_walls_df, use_container_width=True, height=1000)
            else:
                st.dataframe(st.session_state.walls_df, use_container_width=True, height=1000)

        
    # -------------------------------------------------------------------------------------------------------------------
    with column_tab:
        st.dataframe(st.session_state.columns_df, use_container_width=True, height=1000)
        

        
