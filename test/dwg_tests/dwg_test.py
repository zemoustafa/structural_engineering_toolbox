import sys
sys.path.append(r"C:\_Github\structural_engineering_toolbox")
from dwg_tools import read_dwg

dxf_path = r"C:\Users\zeyad.moustafa\OneDrive - Pritchard Francis\Desktop\Automation\DXF tests\StructuralPlan-LEVEL05_W.dxf"
target_layer = "250 - STRUCTURAL COLUMNS HIDDEN LINES"

entities = read_dwg.extract_entities_from_layer(dxf_path, target_layer)
pass