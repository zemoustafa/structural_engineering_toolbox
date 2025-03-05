import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from ram_concept_tools import ram_api

# file paths

file_path_list = [
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 8.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 7.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 6.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 5.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 4.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 3.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250226 - Level 2.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250110 - Level 1 Podium.cpt",
r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Documents\RAM API Tests\Load Rundown\250304 - Upper Ground.cpt"
]

for i in range(len(file_path_list) - 1):
    ram_api.delete_exisitng_loads(file_path_list[i + 1])

    reaction_path = file_path_list[i]
    target_path = file_path_list[i + 1]
    ram_api.ram_load_rundown(reaction_path, target_path)