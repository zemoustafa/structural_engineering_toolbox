import ezdxf
import os

def extract_entities_from_layer(dxf_file_path: str, target_layer: str) -> list[ezdxf.entities.DXFEntity]:
    """
    Extracts all DXF entities found on a specific layer in the modelspace.

    Args:
        dxf_file_path (str): The full path to the DXF file.
        target_layer (str): The name of the layer to filter entities from.

    Returns:
        list[DXFEntity]: A list containing the raw ezdxf entity objects found on the target layer.
    """
    entities_on_layer = [] # Initialize an empty list to store matching entities

    # --- 1. Check if file exists ---
    if not os.path.exists(dxf_file_path):
        print(f"Error: DXF file not found at {dxf_file_path}")
        return entities_on_layer # Return empty list

    # --- 2. Open the DXF file and access modelspace ---
    try:
        doc = ezdxf.readfile(dxf_file_path)
        msp = doc.modelspace()
        print(f"Successfully opened {os.path.basename(dxf_file_path)}. Searching layer '{target_layer}'...")
    except Exception as e:
        # Basic error handling for file reading issues
        print(f"Error reading DXF file {dxf_file_path}: {e}")
        return entities_on_layer # Return empty list on error

    # --- 3. Loop through entities and filter by layer ---
    entity_count = 0
    for entity in msp:
        # Check if the entity has a 'layer' attribute (most common entities do)
        if hasattr(entity.dxf, 'layer'):
            # Check if the entity's layer matches the target layer (case-sensitive)
            if entity.dxf.layer == target_layer:
                # If it matches, add the entire entity object to the list
                entities_on_layer.append(entity)
                entity_count += 1

    print(f"Found {entity_count} entities on layer '{target_layer}'.")
    # --- 4. Return the list of found entities ---
    return entities_on_layer

