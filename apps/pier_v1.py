import streamlit as st
import pandas as pd
import sys

# Add the project directory to the Python path
sys.path.insert(0, '.')

# Import your modules
from etabs_tools import etabs_api, etabs_design_v1
from design_reports import design_reports

# Default load case names
default_eq_flex = "EqX_envelope"
default_eq_shear = "EqX_envelope"
default_wind = "Wind_envelope"

# Initialize session state variables
if 'connected' not in st.session_state:
    st.session_state.connected = False
if 'analysis_run' not in st.session_state:
    st.session_state.analysis_run = False
if 'designed_piers_df' not in st.session_state:
    st.session_state.designed_piers_df = None
if 'load_case_names' not in st.session_state:
    st.session_state.load_case_names = None


st.title("ETABS Pier Design")

# Step 1 - Connect to ETABS button
if st.button("Connect to ETABS"):
    try:
        st.session_state.load_case_names = etabs_api.get_load_cases_and_combinations()
        st.session_state.connected = True
        st.write("Successfully connected to ETABS!")
    except Exception as e:
        st.write(f"Error connecting to ETABS: {e}")

# Input fields for load case names
eq_flex = st.selectbox("Earthquake Flexure Load Case", st.session_state.load_case_names, index=st.session_state.load_case_names.index(default_eq_flex))
eq_shear = st.selectbox("Earthquake Shear Load Case", st.session_state.load_case_names, index=st.session_state.load_case_names.index(default_eq_shear))
wind = st.selectbox("Wind Load Case", st.session_state.load_case_names, index=st.session_state.load_case_names.index(default_wind))



# Run Analysis button (enabled only after connection)
if st.button("Run Analysis", disabled=not st.session_state.connected):
    try:
        piers = etabs_api.get_piers(eq_flex, eq_shear, wind)
        if piers:
            st.session_state.designed_piers_df = etabs_design_v1.design_all_piers(
                piers,
                conc_grade="32 MPa",
                steel_grade="500 MPa",
                fsu_grade="1860 MPa",
                cover=30,
                AS3600_2018_cl11_6_4=True,
                phi_os=0.85,
                k_u=0.36,
            )
            st.session_state.analysis_run = True
            st.write("Analysis completed successfully!")
            st.dataframe(st.session_state.designed_piers_df)
        else:
            st.write("No piers found in the ETABS model.")

    except Exception as e:
        st.write(f"Error during analysis: {e}")


# Generate Report button (enabled only after analysis)
if st.button("Generate Report", disabled=not st.session_state.analysis_run):
    try:
        design_reports.dataframe_to_xlsx(st.session_state.designed_piers_df)
        st.write("Report generated successfully!")
    except Exception as e:
        st.write(f"Error generating report: {e}")