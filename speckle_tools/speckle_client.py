import re
from urllib.parse import urlparse, urljoin
from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport
from specklepy.api.client import SpeckleClient

class speckle_client:
    def __init__(self, url:str, private_token:str):
        """
        Creates instance of speckle client with commit info ready to go from model url

        """
        self.host_server, self.project_id, self.model_id, self.commit_id = extract_ids_from_url(url=url)
        self.private_token = private_token
        self.client = SpeckleClient(host=self.host_server) # app.speckle.systems is for new web app (not experimental server)
        self.client.authenticate_with_token(token=private_token)
        self.commit = self.client.commit.get(stream_id=self.project_id, commit_id=self.commit_id)
        self.commit_referenced_object = self.commit.referencedObject

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
