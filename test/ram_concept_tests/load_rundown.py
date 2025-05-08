import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from ram_concept_tools import ram_api_v1 as ram_api

# file paths

# file_path_list = [
# r"C:\Temp\Load Rundown\Level 8.cpt",
# r"C:\Temp\Load Rundown\Level 7.cpt",
# r"C:\Temp\Load Rundown\Level 6.cpt",
# r"C:\Temp\Load Rundown\Level 5.cpt",
# r"C:\Temp\Load Rundown\Level 4.cpt",
# r"C:\Temp\Load Rundown\Level 3.cpt",
# r"C:\Temp\Load Rundown\Level 2.cpt",
# r"C:\Temp\Load Rundown\Level 1.cpt",
# r"C:\Temp\Load Rundown\Upper Ground.cpt"
# ]

ram_api.delete_existing_loads(r"C:\Temp\Load Rundown\Level 7.cpt")
pass
# for i in range(len(file_path_list) - 1):
#     reaction_path = file_path_list[i]
#     target_path = file_path_list[i + 1]

#     print("Opening file " + reaction_path) # opens reaction path
#     ram_api.delete_exisitng_loads(target_path) # deletes loads in target path

#     print("Applying loads to file " + target_path)
#     ram_api.ram_load_rundown(reaction_path, target_path)
#     print("Finished file " + target_path)
#     print('---------------------------------------------------')