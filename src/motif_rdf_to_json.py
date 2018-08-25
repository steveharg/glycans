from SPARQLWrapper import SPARQLWrapper, JSON
import json
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    sparql_endpoint_url = config['DEFAULT']['SparqlEndpointUrl']
    sparql_query = config['DEFAULT']['MotifIdAndLabelSparqlQuery']
    motif_ids_and_labels_json_file = config['DEFAULT']['MotifIdsAndLabelsJsonFile']

    sparql = SPARQLWrapper(sparql_endpoint_url)
    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    # write json results to file
    with open(motif_ids_and_labels_json_file, 'w') as f:
        json.dump(results, f, ensure_ascii=False)

