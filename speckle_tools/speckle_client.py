import re
from urllib.parse import urlparse, urljoin
from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport
from specklepy.api import operations
from specklepy.api.client import SpeckleClient
from specklepy.transports.server import ServerTransport

class speckle_client:
    def __init__(self, url:str, private_token:str):
        """
        Creates instance of speckle client with commit info ready to go from model url

        :param url: url of speckle model on latest server
        :param private_token: private token from speckle account
        """
        self.host_server, self.project_id, self.model_id, self.commit_id = extract_ids_from_url(url=url)
        self.private_token = private_token
        self.client = SpeckleClient(host=self.host_server) # app.speckle.systems is for new web app (not experimental server)
        self.client.authenticate_with_token(token=private_token)
        self.commit = self.client.commit.get(stream_id=self.project_id, commit_id=self.commit_id)
        self.commit_referenced_object = self.commit.referencedObject


    def server_transport(self):
        """
        Create server transport and extract elements from speckle model
        
        """
        # create an authenticated server transport from the client and receive the commit obj
        self.transport = ServerTransport(client=self.client, stream_id=self.project_id)

        # extract data from model
        self.model_data = operations.receive(self.commit_referenced_object, remote_transport=self.transport)
        self.elements = self.model_data["elements"]

    def get_revit_floors(self):
        """
        Extract floor coordinates 
        
        :return revit_floors: list of dicts representing each floor. each X, Y and Z in dict is a list of coordinates
        """
        revit_floors = None  # Initialize to None
        for element in self.elements:
            if element['name'] == 'Floors':
                revit_floors = element['elements']
                break  # Exit the loop once found
        
        floors = []

        for index, floor in enumerate(revit_floors):
            if floor.speckle_type == "Objects.BuiltElements.Floor:Objects.BuiltElements.Revit.RevitFloor":
                # grab outline of floor 
                segments = floor.outline.segments
                num_segments = len(segments)

                # intialise coorindates list
                X = [] 
                Y = []
                Z = []

                for i, segment in enumerate(segments):
                    # grab start and end points from current segment
                    # remember - start of segment 2 = end of segment 1 
                    
                    # find if segment is line or arc
                    if  segment.speckle_type == "Objects.Geometry.Line":
                        start = segment.start
                    else:
                        start = segment.startPoint

                    # check if on the final iteration of the segments loop or not
                    if i == num_segments:
                        # if on the final iteration, point = start of first segment
                        x = X[0]
                        y = Y[0]
                        z = Z[0]
                    else:
                        # all other iterations = grab start point of the segment
                        x = round(start.x, 0)
                        y = round(start.y, 0)
                        z = round(floor.level.elevation, 0)

                    X.append(x)
                    Y.append(y)
                    Z.append(z)

                floors.append({
                    'X': X,
                    'Y': Y,
                    'Z': Z,
                    'Name': "floor {index}"
                })

        return floors

    def get_revit_columns(self):
        """
        Extract floor coordinates 
        
        :return revit_floors: list of dicts representing each column

        """

        # grab revit columns from self.elements
        revit_columns = None  # Initialize to None
        for element in self.elements:
            if element['name'] == 'Structural Columns':
                revit_columns = element['elements']
                break  # Exit the loop once found

        for index, column in enumerate(revit_columns):   
            # grab start and end points of column
            start = column.baseLine.start # grab x, y and z at bottom of column
            end = column.baseLine.end # grab x, y and z at top of column

            # grab column base and top level names
            baseLevel = column.level.name
            topLevel = column.topLevel.name

            # iterate within range of start and end level
            in_range = False
            for index, level in enumerate(story_names):
                if level == baseLevel:
                    in_range = True
                if in_range:
                    # extract start and end coordinates of column
                    startX = round(start.x, 0)
                    startY = round(start.y, 0)
                    startZ = story_elevations[index]
                    endX = np.round(end.x, 0)
                    endY = np.round(end.y, 0)
                    endZ = story_elevations[index + 1]
                    # add column with ETABS API
                    ret = SapModel.FrameObj.AddByCoord(startX, startY, startZ, endX, endY, endZ)
                if story_names[index + 1] == topLevel:
                    in_range = False
                    break

    def gql_client(self):
        self.gql_client = Client(
            transport=RequestsHTTPTransport( url=f"{self.host_server}graphql" )
        )


    def get_unique_levels(self):
        query = gql(
            """ 
            query Commit( $stream_id: String!, $commit_referenced_object: String! ) { 
                stream(id: $stream_id) { 
                    object(id: $commit_referenced_object) {
                        children( limit:1000, 
                            select: [
                                "level.name", 
                                "level.parameters.LEVEL_ELEV.value"
                            ]
                        ) { objects { data } } } } } """
        )

        params = {
            "stream_id": self.project_id, 
            "commit_referenced_object": self.commit_referenced_object
        }

        received_data = self.gql_client.execute(query, variable_values=params)
        
        # Initialize an empty set to store unique levels
        unique_levels = {}

        # Navigate through the nested structure to extract level names and elevations
        if (received_data and 'stream' in received_data and 
                'object' in received_data['stream']):
            objects = received_data['stream']['object']\
                .get('children', {}).get('objects', [])
            for obj in objects:
                level_data = obj.get('data', {}).get('level', {})
                name = level_data.get('name', 'Unknown')
                elevation = level_data.get('parameters', {})\
                    .get('LEVEL_ELEV', {}).get('value', 'Unknown')
                unique_levels[name] = round(elevation, 0)

        sorted_unique_levels = dict(sorted(unique_levels.items(), key=lambda item: item[1]))
        return sorted_unique_levels


def extract_ids_from_url(url):
    """ Extract project ID, model ID and commit ID straight from URL of Speckle model on new web app
    
    """
    # Parse the URL to extract the host and path
    parsed_url = urlparse(url)
    host_server = f"{parsed_url.scheme}://{parsed_url.netloc}/"
    
    # Define the regex pattern to capture the project_id, model_id, and commit_id
    pattern = r"/projects/([a-f0-9]{10})/models/([a-f0-9]{10})@([a-f0-9]{10})"
    
    # Use regular expression to find matches
    match = re.search(pattern, parsed_url.path)
    
    if match:
        project_id = match.group(1)
        model_id = match.group(2)
        commit_id = match.group(3)
        return host_server, project_id, model_id, commit_id
    else:
        raise ValueError("URL format is incorrect or does not match the expected pattern")
