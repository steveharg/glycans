from SPARQLWrapper import SPARQLWrapper, JSON
import json
import configparser

config = configparser.ConfigParser(allow_no_value=True)
config.read('config.ini')

if __name__ == "__main__":

    sparql_endpoint_url = config['DEFAULT']['SparqlEndpointUrl']
    sparql_query = config['DEFAULT']['GlycanAndMotifSparqlQuery']
    glycans_json_file = config['DEFAULT']['GlycansJsonFile']

    sparql = SPARQLWrapper(sparql_endpoint_url)
    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        if "MotifPrimaryId" in result:
            print(result["Sequence"]["value"], " ", result["MotifPrimaryId"]["value"])
        else:
            print(result["Sequence"]["value"])

    # write json results to file
    with open(glycans_json_file, 'w') as f:
        json.dump(results, f, ensure_ascii=False)