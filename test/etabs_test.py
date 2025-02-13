import sys
sys.path.append('C:\_Github\structural_engineering_toolbox')
from etabs_tools import etabs_api, etabs_design 
from design_reports import design_reports

etabs_api = etabs_api.etabs_api()

eq_env_1 = '(88) RS ULS ENV' # earthquake envelope unfactored EQ for moment design
eq_env_2 = '(88) RS ULS ENV SHEAR' # earthquake envelope factored EQ for amplified shear design
wind_env = '(88) WIND ULS ENV' # wind envelope

piers = etabs_api.get_piers(load_cases=[eq_env_1, eq_env_2, wind_env])
designed_piers_df = etabs_design.design_all_piers(
    piers=piers,
    eq_env_1=eq_env_1,
    eq_env_2=eq_env_2,
    wind_env=wind_env,
    vertical_spacing=200,
    horizontal_spacing=200,
    design_both_axes=False
)

design_reports.dataframe_to_xlsx(designed_piers_df)
