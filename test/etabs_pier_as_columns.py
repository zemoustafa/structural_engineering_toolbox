# This script is used to design piers in ETABS as columns
# The script will design the piers for the selected load combinations

# Import libraries
import sys
sys.path.append('C:\_Github\structural_engineering_toolbox')
from etabs_tools import etabs_api
from design_functions import as3600_column_design

import pandas as pd

# Select names of load combinations to design for
eq_env_1 = '(88) RS ULS ENV' # earthquake envelope unfactored EQ for moment design
eq_env_2 = '(88) RS ULS ENV SHEAR' # earthquake envelope factored EQ for amplified shear design
wind_env = '(88) WIND ULS ENV' # wind envelope
gravity_env = ' ' # gravity envelope

load_cases = [eq_env_1, eq_env_2, wind_env]

# Select spacing for vertical and horizontal bars
vertical_spacing = 200
horizontal_spacing = 200

# Create ETABS API object
etabs_api = etabs_api.etabs_api()
sap_model = etabs_api.sap_model

# Get piers with pier forces from ETABS
piers = etabs_api.get_piers(load_cases=load_cases)

# Design all piers as columns
designed_piers = []
for pier in piers:

    designed_pier = as3600_column_design.design_etabs_pier_as_column(pier, eq_env_1, eq_env_2, wind_env, vertical_spacing, horizontal_spacing, False)
    designed_piers.append(designed_pier)

df = pd.DataFrame(designed_piers)
pass
