"""
This module gets or sets data in an ETABS model using the API.
"""

import comtypes.client

class etabs_api:
    def __init__(self):
        """
        Creates instance of sap model object

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
        self.story_names = self.get_story_names()
        self.story_data = self.get_story_data()

    def get_piers(self, load_cases: list) -> list[dict]:
        """ Get wall/pier property assignments and pier forces for selected load cases
        
        :param load_cases: list of load case names in etabs model
        :return walls: list of dicts. each dict contains wall properties, geometry and forces
        """
        piers = []
        pier_forces = self.get_pier_forces(load_cases)

        for section_property in self.pier_section_properties:
            for story in self.story_data:
                if section_property["story_name"] == story["story_name"]:
                    section_property["story_height"] = round(story["story_height"] * 1000, 0)

            for material_prop in self.concrete_material_properties:
                if section_property["mat_prop"] == material_prop["mat_name"]:
                    section_property["fc"] = material_prop["fc"]

            for pier_force in pier_forces:
                if (
                    section_property["story_name"] == pier_force["story_name"]
                    and section_property["pier_name"] == pier_force["pier_name"]
                ):
                    section_property[pier_force["load_case"]] = pier_force

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
            number_results = pier_force[0]
            for i in range(number_results):
                location = pier_force[4][i]
                if location == "Bottom":
                    pier_force_dict = {
                        "story_name": pier_force[1][i],
                        "pier_name": pier_force[2][i],
                        "load_case": pier_force[3][i],
                        "location": pier_force[4][i],
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

        :return story_names: list of story names defined in model
        """
        story_names = []
        name_list = self.sap_model.Story.GetNameList()[1]
        for name in name_list:
            story_names.append(name)
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
                "story_name": story_names[i],
                "story_elevation":story_elevations[i], 
                "story_height": story_heights[i]}
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
                        "pier_name": name,
                        "story_name": piers[1][i],
                        "axis_angle": piers[2][i],
                        "num_area_objs": piers[3][i],
                        "num_line_objs": piers[4][i],
                        "width_bot": round(piers[5][i] * 1000, 0),
                        "thickness_bot": piers[6][i] * 1000,
                        "width_top": round(piers[7][i] * 1000, 0),
                        "thickness_top": piers[8][i] * 1000,
                        "mat_prop": piers[9][i],
                        "cg_bot_x": piers[10][i],
                        "cg_bot_y": piers[11][i],
                        "cg_bot_z": piers[12][i],
                        "cg_top_x": piers[13][i],
                        "cg_top_y": piers[14][i],
                        "cg_top_z": piers[15][i],
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
                    "mat_name": mat_name,
                    "fc": material[0] / 1000,
                    "is_lightweight": material[1],
                    "fcs_factor": material[2],
                    "ss_type": material[3],
                    "ss_hys_type": material[4],
                    "strain_at_fc": material[5],
                    "strain_ult": material[6],
                    "final_slope": material[7],
                    "friction_angle": material[8],
                    "dilation_angle": material[9],
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

    def set_wall_modifiers():
        pass