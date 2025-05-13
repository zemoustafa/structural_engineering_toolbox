import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from etabs_tools import etabs_api
from ram_concept_tools import ram_api_v1 as ram_api

# Setup
model_path = r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Raft Test.cpt"
load_case_names = ['V:G', 'V:Q'] # ETABS load case names

# Open etabs api and grab columns from lowest story
etabs_api = etabs_api.etabs_api()
columns = etabs_api.get_columns_on_story('Story1', load_case_names)

# Create instance of Concept and open a model
ram_concept = ram_api.RamConcept()
ram_concept.create_concept()
print('Concept created')

model = ram_concept.open_model(model_path)
print('Model opened')

print('Drawing columns in RAM Concept')
for column in columns:
    print('Column: ' + column['Name'])
    ram_api.create_column(model, column, load_case_names)

# Save file
model.save_file(model_path)
print('File saved')

# Shutdown concept
ram_concept.shutdown_concept()
print('Concept shut down')
