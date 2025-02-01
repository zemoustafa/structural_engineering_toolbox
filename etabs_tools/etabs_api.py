"""
This module gets or sets data in an ETABS model using the API.
"""

import comtypes.client

class etabs_api:
    def __init__(self):
        """
        Creates instance of sap model object to access data in etabs model

        """
        self.helper = comtypes.client.CreateObject("ETABSv1.Helper")
        self.helper = self.helper.QueryInterface(comtypes.gen.ETABSv1.cHelper)

        try:
            self.etabs_object = self.helper.GetObject("CSI.ETABS.API.ETABSObject")
        except (OSError, comtypes.COMError):
            raise Exception(
                "No running instance of the program found or failed to attach."
            )

        self.sap_model = self.etabs_object.SapModel

        self.pier_section_properties = self.get_pier_section_properties()
        self.concrete_material_properties = self.get_concrete_material_properties()
        self.load_cases = self.get_load_cases()
        self.story_data = self.get_story_data()
        self.story_names = [story_data["Story Name"] for story_data in self.story_data]

    """
    GET FUNCTIONS
    
    Retrieve properties, assignments and analysis results from etabs model
    
    """

    def get_piers(self, load_cases: list) -> list[dict]:
        """ Get wall/pier property assignments and pier forces for selected load cases
        
        :param load_cases: list of load case names in etabs model
        :return walls: list of dicts. each dict contains wall properties, geometry and forces
        """
        piers = []
        pier_forces = self.get_pier_forces(load_cases)

        # Create a dictionary for quick lookup of story heights and material properties
        story_heights = {story["Story Name"]: round(story["Story Height"], 0) for story in self.story_data}
        material_fc = {material["Mat Name"]: material["fc"] for material in self.concrete_material_properties}

        for section_property in self.pier_section_properties:
            section_property["Story Height"] = story_heights.get(section_property["Story Name"], 0)
            section_property["fc"] = material_fc.get(section_property["Mat Prop"], 0)
            pier_force = {}
            for pier_force in pier_forces:
                if section_property["Story Name"] == pier_force["Story Name"] and section_property["Pier Name"] == pier_force["Pier Name"]:
                    key = f"{pier_force['Load Case']} {pier_force['Location']} {pier_force['Step Type']}"
                    section_property[key] = pier_force

            piers.append(section_property)

        return piers

    def get_pier_forces(self, load_cases: list) -> list[dict]:
        """ Get pier property forces for selected load cases

        :param load_cases: list of load case names
        :return walls: list of dicts containing pier forces for each load case in list
        """
        pier_force_list = []
        for case in load_cases:
            self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
            result = self.sap_model.Results.Setup.SetCaseSelectedForOutput(case)
            if result == 1:
                self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
                result = self.sap_model.Results.Setup.SetComboSelectedForOutput(case)
            pier_force = self.sap_model.Results.PierForce()

            max_min_order = ['Max', 'Max', 'Min', 'Min']

            number_iterations = pier_force[0]  # Total number of iterations
            number_piers = number_iterations // len(max_min_order)  # Number of unique piers

            for j in range(number_piers):  # Loop over each unique pier
                for k in range(len(max_min_order)):  # Loop over Max/Min order
                    i = j * len(max_min_order) + k  # Correct index calculation
                    pier_force_dict = {
                        "Story Name": pier_force[1][i],
                        "Pier Name": pier_force[2][i],
                        "Load Case": pier_force[3][i],
                        "Step Type": max_min_order[k],
                        "Location": pier_force[4][i],
                        "p": pier_force[5][i],
                        "v2": pier_force[6][i],
                        "v3": pier_force[7][i],
                        "t": pier_force[8][i],
                        "m2": pier_force[9][i],
                        "m3": pier_force[10][i],
                    }
                    pier_force_list.append(pier_force_dict)
        return pier_force_list

    def get_story_names(self) -> list[str]:
        """ Get list of story names

        :return story_data: list of dicts containing story data
        """
        stories = self.sap_model.Story.GetStories_2()
        story_names = stories[2]
        return story_names
    
    def get_story_data(self) -> list[dict]:
        """ Get list of story names, elevations and heights

        :return story_data: list of dicts containing story data
        """
        story_data = []
        stories = self.sap_model.Story.GetStories_2()
        story_names = stories[2]
        story_elevations = stories[3]
        story_heights = stories[4]
        number_stories = stories[1]
        for i in range(number_stories):
            story_data.append(
                {
                "Story Name": story_names[i],
                "Story Elevation":story_elevations[i], 
                "Story Height": story_heights[i]}
            )
        return story_data

    def get_pier_section_properties(self) -> list[dict]:
        """ Get section properties of piers

        :return section_properties: list dicts containing section properties of all defined piers in model
        """
        section_properties = []
        name_list = self.sap_model.PierLabel.GetNameList()[1]
        for name in name_list:
            piers = self.sap_model.PierLabel.GetSectionProperties(name)
            number_stories = piers[0]
            for i in range(number_stories):
                section_properties.append(
                    {
                        "Pier Name": name,
                        "Story Name": piers[1][i],
                        "Axis Angle": piers[2][i],
                        "Num Area Objs": piers[3][i],
                        "Num Line Objs": piers[4][i],
                        "Width Bot": round(piers[5][i], 0),
                        "Thickness Bot": piers[6][i],
                        "Width Top": round(piers[7][i], 0),
                        "Thickness Top": piers[8][i],
                        "Mat Prop": piers[9][i],
                        "CG Bot X": piers[10][i],
                        "CG Bot Y": piers[11][i],
                        "CG Bot Z": piers[12][i],
                        "CG Top X": piers[13][i],
                        "CG Top Y": piers[14][i],
                        "CG Top Z": piers[15][i],
                    }
                )
        return section_properties

    def get_load_case_list(self) -> list[str]:
        """ Get list of load case names  (excludes load combinations)

        :return load_cases_list: list of load cases
        """
        load_cases_list = []
        load_cases = self.sap_model.LoadCases.GetNameList()[1]
        for load_case in load_cases:
            load_cases_list.append(load_case)
        return load_cases_list

    def get_load_cases(self) -> list[str]:
        """ Get list of load case names (including load combinations)

        :return load_cases: list of load cases
        """
        load_cases_list = []
        load_cases = self.sap_model.LoadCases.GetNameList()[1]
        for load_case in load_cases:
            load_cases_list.append(load_case)
        load_cases_combo = self.sap_model.RespCombo.GetNameList()[1]
        for load_case_combo in load_cases_combo:
            load_cases_list.append(load_case_combo)
        return load_cases_list

    def get_concrete_material_properties(self) -> list[dict]:
        """ Get concrete material properties

        :return conc_mat: list of dicts containing conc material properties
        """
        conc_mat = []
        mat_names = self.sap_model.PropMaterial.GetNameList()[1]
        # mat type for concrete = 2
        for mat_name in mat_names:
            material = self.sap_model.PropMaterial.GetOConcrete_1(mat_name)
            number_mat = len(mat_names)
            conc_mat.append(
                {
                    "Mat Name": mat_name,
                    "fc": material[0],
                    "Is Lightweight": material[1],
                    "fcs Factor": material[2],
                    "SS Type": material[3],
                    "SS Hys Type": material[4],
                    "Strain at fc": material[5],
                    "Strain Ult": material[6],
                    "Final SLope": material[7],
                    "Friction Angle": material[8],
                    "Dilation Angle": material[9],
                }
            )
        return conc_mat

    def get_modal_results(self, modal_case:str) -> dict:
        """ Get modal analysis results

        :params modal_case: name of modal load case
        :return load_cases_list: list of dicts containing modal results
        """
        self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
        self.sap_model.Results.Setup.SetCaseSelectedForOutput(modal_case)
        modal_results = self.sap_model.Results.ModalParticipatingMassRatios()
        modes = [f'Mode {i}' for i in range(1, modal_results[0] + 1)]
        modal_data = {
            'Modes': modes,
            'Period':modal_results[4], 
            'UX':modal_results[5], 
            'UY':modal_results[6], 
            'UZ':modal_results[7], 
            'SumUX':modal_results[8],
            'SumUY':modal_results[9],
            'SumUZ':modal_results[10],
            'RX':modal_results[11], 
            'RY':modal_results[12], 
            'RZ':modal_results[13], 
            'SumRX':modal_results[14],
            'SumRY':modal_results[15],
            'SumRZ':modal_results[16],
        }
        return modal_data

    def get_story_drifts(self, load_case:str) -> list[dict]:
        """ Get modal analysis results

        :params load_case: name of load case
        :return story_drifts_data: list of dicts containing story drifts for selected load case
        """
        story_drifts_data = []
        self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
        self.sap_model.Results.Setup.SetCaseSelectedForOutput(load_case)
        story_drifts = self.sap_model.Results.StoryDrifts()

        for i in range(story_drifts[0]):
            story_drifts_data.append({
                'Story': story_drifts[1][i],
                'Load Case': story_drifts[2][i],
                'Step Type': story_drifts[3][i],
                'Step Num': story_drifts[4][i],
                'Direction': story_drifts[5][i],
                'Drift': story_drifts[6][i],
                'Label': story_drifts[7][i],
                'Disp X': story_drifts[8][i]/1000,
                'Disp Y': story_drifts[9][i]/1000,
                'Disp Z': story_drifts[10][i]
            })
        
        stories = self.get_story_data()
        
        for drift in story_drifts_data:
            for story in stories:
                if drift["Story"] == story["story_name"]:
                    drift["Height"] = story["story_height"]
                    drift["Elevation"] = story["story_elevation"]

        return story_drifts_data

    def get_story_forces(self, load_case:str) -> list[dict]:
        # story_forces = []
        # self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
        # self.sap_model.Results.Setup.SetCaseSelectedForOutput(load_case)
        # reactions = self.sap_model.Results.JointReact("All", 2)
        # # for reaction in reactions:
        # return reactions
        pass

    def get_base_reactions(self, load_case:str):
        """ Get base reactions
        :params load_case: name of load case
        :return base_reactions_data: base reactions of selected load case
        """
        self.sap_model.Results.Setup.DeselectAllCasesAndCombosForOutput()
        self.sap_model.Results.Setup.SetCaseSelectedForOutput(load_case)
        base_reactions = self.sap_model.Results.BaseReact()
        base_reactions_data = {
            'Load Case': base_reactions[1],
            'Step Type': base_reactions[2],
            'Step Num': base_reactions[3],
            'FX': base_reactions[4],
            'FY': base_reactions[5],
            'FZ': base_reactions[6],
            'MX': base_reactions[7],
            'ParamMy': base_reactions[8],
            'Mz': base_reactions[9],
            'FZ': base_reactions[10],
            'FZ': base_reactions[11],
        }
        return base_reactions_data
    
    def get_floor_objs(self) -> list[dict]:
        """ Get floor geometry and properties

        :return floor_objs: floor geometry and assignments
        """
        area_objs_labels = self.sap_model.AreaObj.GetLabelNameList()
        floor_objs = []
        for i in range(area_objs_labels[0]):
            if area_objs_labels[2][i].startswith('F'):
                floor_objs.append({
                    'Name': area_objs_labels[1][i],
                    'Label': area_objs_labels[2][i],
                    'Story': area_objs_labels[3][i]
                })  
        
        for floor in floor_objs:
            # add floor properties
            floor_prop = self.sap_model.AreaObj.GetProperty(floor['Name'])
            floor['Property'] = floor_prop[0]

            # add diaphragm assignment
            floor_diaphragm = self.sap_model.AreaObj.GetDiaphragm(floor['Name'])
            floor['Diaphragm'] = floor_diaphragm[0]

            # add if floor is opening or not
            is_opening = self.sap_model.AreaObj.GetOpening(floor['Name'])
            floor['Opening'] = is_opening[0]

            # add floor points and x, y and z coordinates
            floor_points = self.sap_model.AreaObj.GetPoints(floor['Name'])
            floor['Num Points'] = floor_points[0]
            floor['Points'] = floor_points[1]
            point_coords_x = []
            point_coords_y = []
            point_coords_z = []
            for point in floor['Points']:
                coords = self.sap_model.PointObj.GetCoordCartesian(point)
                point_coords_x.append(coords[0])
                point_coords_y.append(coords[1])
                point_coords_z.append(coords[2])
            floor['Point X'] = point_coords_x
            floor['Point Y'] = point_coords_y
            floor['Point Z'] = point_coords_z
            
            # put the first coordinate at the end of the list to close the loop
            floor['Point X'].append(floor['Point X'][0])
            floor['Point Y'].append(floor['Point Y'][0])
            floor['Point Z'].append(floor['Point Z'][0])

        return floor_objs

    def get_wall_objs(self) -> list[dict]:
        """ Get wall geometry

        :return wall_objs
        """
        area_objs_labels = self.sap_model.AreaObj.GetLabelNameList()
        wall_objs = []
        for i in range(area_objs_labels[0]):
            if area_objs_labels[2][i].startswith('W'):
                wall_objs.append({
                    'Name': area_objs_labels[1][i],
                    'Label': area_objs_labels[2][i],
                    'Story': area_objs_labels[3][i]
                })  
        
        for wall in wall_objs:
            # add floor points and x, y and z coordinates
            wall_points = self.sap_model.AreaObj.GetPoints(wall['Name'])
            wall['Num Points'] = wall_points[0]
            wall['Points'] = wall_points[1]
            point_coords_x = []
            point_coords_y = []
            for point in wall['Points']:
                coords = self.sap_model.PointObj.GetCoordCartesian(point)
                point_coords_x.append(coords[0])
                point_coords_y.append(coords[1])

            wall['X'] = [min(point_coords_x), max(point_coords_x)]
            wall['Y'] = [min(point_coords_y), max(point_coords_y)]

        return wall_objs
    
    def get_column_objs(self) -> list[dict]:
        """ Get column locations on each story
        
        """


        pass

    """
    CREATE/DRAW FUNCTIONS
    Functions that will create properties or draw objects in etabs model using their geometry and property assignments
 
    """
    def create_concrete_material_property():
        """ Create new concrete material property
        
        """
        pass

    def create_concrete_column_property(self, width:float, depth:float, bars_x:int, bars_y:int, db:int, material:str, reo_mat:str) -> int:
        """ Create new concrete column frame property
        
        :param width: width or column section (shorter dimension). if depth = 0, width = diameter of column
        :param depth: depth of column section (longer dimension). if depth = 0, create circular column
        :param bars_x: number of bars along the x (3) direction
        :param bars_y: number of bars along the y (2) direction
        :param db: diameter of vertical bars
        :param material: name of material property to assign to column
        :param reo_mat: name of reinforcement material property to assign to column bars

        :type width: float
        :type depth: float
        :type bars_x: int
        :type bars_y: int
        :type db: int
        :type material: str
        :type reo_mat: str

        :return result: returns 0 if column successfully created, otherwise returns 1
        """
        material_props = self.sap_model.PropMaterial.GetOConcrete_1(material)
        fc = round(material_props[0] / 1000)

        # constants
        cover = 50
        tie_dia = 12
        tie_cts = 300
        bars_circle = 0
        num_2_dir_bars, num_3_dir_bars = 0, 0

        if depth == 0: # circular column
            name = "C"+str(width)+"C"+str(fc)
            pattern = 2
            confine_type = 2
            conc_column = self.sap_model.PropFrame.SetCircle(name, material, width/1000)
            bars_circle = bars_x
            bars_x, bars_y = 0, 0

        else: # rectangular column
            name = "C"+str(width)+"X"+str(depth)+"C"+str(fc)
            pattern = 1
            confine_type = 1
            conc_column = self.sap_model.PropFrame.SetRectangle(name, material, width/1000, depth/1000)
            num_2_dir_bars = bars_y
            num_3_dir_bars = bars_x

        set_reo = self.sap_model.PropFrame.SetRebarColumn(
            name, 
            reo_mat, 
            reo_mat, 
            pattern, 
            confine_type, 
            cover/1000, 
            bars_circle, 
            bars_x, 
            bars_y, 
            "N"+str(db), 
            "N"+str(tie_dia),
            tie_cts/1000,
            num_2_dir_bars,
            num_3_dir_bars,
            False
            )
        
        result = 0 if set_reo == 0 and conc_column == 0 else 1

        return result

    def draw_floors_by_points():
        pass

    def draw_walls_by_points():
        pass

    def draw_column():
        pass

    def draw_beam():
        pass

    """
    SET FUNCTIONS
    Functions that will set options or properties in etabs model

    """
    

    def set_stories(self, story_data:dict) -> int:
        """ Define story names and elevations in etabs model

        :param story_data: dict where each key is the story name and each value is the story elevation

        :return result: returns 0 if column successfully created, otherwise returns 1
        """
        first_pair = next(iter(story_data.items())) # grab first key value pair
        base_elev = first_pair[1]
        num_stories = len(story_data)
        story_names = list(story_data.keys())
        story_elevations = list(story_data.values())
        story_heights = [story_elevations[0]] + [story_elevations[i] - story_elevations[i-1] for i in range(1, len(story_elevations))]
        is_master = [False] * num_stories
        similar_story = ["None"] * num_stories
        splice_above = [False] * num_stories
        splice_h = [0] * num_stories

        result = self.sap_model.Story.SetStories_2(
        base_elev,
        num_stories-1,
        story_names[1:], 
        story_heights[1:], 
        is_master[1:], 
        similar_story[1:], 
        splice_above[1:], 
        splice_h[1:]
        )

        return result[-1]


    def set_units():
        pass
