import sys
sys.path.append('C:\_Github\structural_engineering_toolbox')
from etabs_tools import etabs_api, etabs_design_v1
from design_reports import design_reports

# Button 1 - connect to API and grab load case names
# Call function etabs_api.get_load_cases_and_combinations() to get load combination names as a list of dic

etabs_api = etabs_api.etabs_api()

# Once we have the load case names, let the user select the relevant case for each envelope from a list or dropdown
# Envelopes need to be Earthquake Flexure, Earthquake Shear, and Wind
eq_env_1 = '(88) RS ULS ENV' # Earthquake Flexure -> earthquake envelope unfactored EQ for moment design
eq_env_2 = '(88) RS ULS ENV mu=1' # Earthquake Shear -> earthquake envelope factored EQ for amplified shear design
wind_env = '(88) WIND ULS ENV' # Wind -> wind envelope

# Once load cases are selected:
# Button 2 - grab pier properties and forces - and run the analysis
try:
    piers = etabs_api.get_piers(load_cases=[eq_env_1, eq_env_2, wind_env])
except Exception as e:
    # raise error if piers are not found
    raise ValueError(f'Invalid load case names: {e}')

designed_piers_df = etabs_design_v1.design_all_piers(
    piers=piers,
    eq_env_1=eq_env_1,
    eq_env_2=eq_env_2,
    wind_env=wind_env,
    vertical_spacing=200,
    horizontal_spacing=200,
    design_both_axes=False,
    design_boundary_element=True
)

# Once it has been run, display the summary table of pier design results

# Button 3 - generate report
design_reports.dataframe_to_xlsx(designed_piers_df)

# Next steps:
# View indiviudal pier design diagrams when clicking on a row in the summary table
# Add summary section, eg number of walls, number of walls pass/fdailed, etc
# Locations of walls on plan?
# Power BI dashboard?