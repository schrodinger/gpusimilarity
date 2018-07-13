import requests
from multicorn import ForeignDataWrapper

FASTSIM_ENDPOINT = '{0}:{1}/similarity_search_json_{2}'


class GPUSimilarityFDW(ForeignDataWrapper):

    def __init__(self, options, columns):
        super(GPUSimilarityFDW, self).__init__(options, columns)
        self.columns = columns
        input_smiles = options['smiles']
        # We can't rely on 'query' for this, as execute() isn't given count
        self.max_results = int(options['max_results'])
        input_db = 'all'
        if 'db_name' in options:
            input_db = options['db_name']
        server = options['server']
        port = options['port']

        data = {'smiles': input_smiles, 'return_count': self.max_results}

        endpoint = FASTSIM_ENDPOINT.format(server, port, input_db)
        response = requests.post(endpoint, data)
        if not response.ok:
            raise RuntimeError("Server connection failed")
        self.data = response.json()

    def execute(self, quals, columns):
        for line_data in self.data:
            line = {}
            line['id'] = line_data[0]
            line['smiles'] = line_data[1]
            line['similarity'] = line_data[2]
            yield line
