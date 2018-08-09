import requests
from multicorn import ForeignDataWrapper

GPUSIM_ENDPOINT = '{0}:{1}/similarity_search_json_{2}'


class GPUSimilarityFDW(ForeignDataWrapper):

    def __init__(self, options, columns):
        super(GPUSimilarityFDW, self).__init__(options, columns)
        self.columns = columns
        # We can't rely on 'query' for this, as execute() isn't given count
        self.max_results = int(options['max_results'])
        input_db = 'all'
        if 'db_name' in options:
            input_db = options['db_name']
        server = options['server']
        port = options['port']

        self.query = {'return_count': self.max_results}
        self.endpoint = GPUSIM_ENDPOINT.format(server, port, input_db)
        self._last_query_smiles = None

    def execute(self, quals, columns):
        smiles = None
        for qual in quals:
            if qual.field_name == 'query' and qual.operator == '=':
                smiles = qual.value
                break
        if smiles is None:
            print("No query smiles given, so no results can be given.")
            raise StopIteration

        if smiles != self._last_query_smiles:
            self.query['smiles'] = smiles
            response = requests.post(self.endpoint, self.query)
            if not response.ok:
                raise RuntimeError("Server connection failed")
            self.data = response.json()
            self._last_query_smiles = smiles
        for line_data in self.data:
            line = {}
            line['id'] = line_data[0]
            line['query'] = smiles
            line['smiles'] = line_data[1]
            line['similarity'] = line_data[2]
            yield line
