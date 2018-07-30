from SPARQLWrapper import SPARQLWrapper, JSON
import json

sparql = SPARQLWrapper("https://ts.glytoucan.org/sparql")
sparql.setQuery("""
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>

    SELECT DISTINCT ?Motif ?MotifPrimaryId ?MotifLabel
    FROM <http://rdf.glytoucan.org/core>
    FROM <http://rdf.glytoucan.org/sequence/wurcs>
    FROM <http://rdf.glytoucan.org/motif>
    WHERE {
            ?Saccharide glycan:has_motif ?Motif .
            ?Motif glytoucan:has_primary_id ?MotifPrimaryId.
            ?Motif rdfs:label ?MotifLabel
    }
    ORDER BY ?MotifPrimaryId
""")
sparql.setReturnFormat(JSON)
# sparql.setReturnFormat(N3)
results = sparql.query().convert()

# write json results to file
with open('motif_ids_and_labels.json', 'w') as f:
    json.dump(results, f, ensure_ascii=False)

