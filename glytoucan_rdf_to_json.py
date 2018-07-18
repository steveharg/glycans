from SPARQLWrapper import SPARQLWrapper, JSON, N3
import json

sparql = SPARQLWrapper("https://ts.glytoucan.org/sparql")
sparql.setQuery("""
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>
    
    SELECT DISTINCT ?Saccharide ?PrimaryId ?Sequence ?Motif ?MotifPrimaryId
    FROM <http://rdf.glytoucan.org/core>
    FROM <http://rdf.glytoucan.org/sequence/wurcs>
    FROM <http://rdf.glytoucan.org/motif>
    WHERE {
        ?Saccharide glytoucan:has_primary_id ?PrimaryId .
        ?Saccharide glycan:has_glycosequence ?GlycoSequence .
        ?GlycoSequence glycan:has_sequence ?Sequence .
        ?GlycoSequence glycan:in_carbohydrate_format glycan:carbohydrate_format_wurcs.
        OPTIONAL { ?Saccharide glycan:has_motif ?Motif .
                   ?Motif glytoucan:has_primary_id ?MotifPrimaryId } .
    }
    ORDER BY ?PrimaryId
    # limit 100
""")
sparql.setReturnFormat(JSON)
# sparql.setReturnFormat(N3)
results = sparql.query().convert()

for result in results["results"]["bindings"]:
    if "MotifPrimaryId" in result:
        print(result["Sequence"]["value"], " ", result["MotifPrimaryId"]["value"])
    else:
        print(result["Sequence"]["value"])

# write json results to file
with open('glycans.json', 'w') as f:
  json.dump(results, f, ensure_ascii=False)

# Explore this query - breaks the WURCS string down into its component parts:
# PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
# PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>
#
# SELECT DISTINCT ?Saccharide ?p ?o ?gp ?go
# FROM <http://rdf.glytoucan.org/core>
# FROM <http://rdf.glytoucan.org/sequence/wurcs>
# WHERE {
#     ?Saccharide glytoucan:has_primary_id 'G00005DE'.
#     ?Saccharide glycan:has_glycosequence ?GlycoSequence .
# #    ?GlycoSequence glycan:has_sequence ?Sequence .
#     ?GlycoSequence ?gp ?go .
# #    ?GlycoSequence glytoucan:has_motif ?motif .
#     ?Saccharide ?p ?o
# }
# order by ?Saccharide
# #limit 10

# and this:
# PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
# PREFIX glytoucan:  <http://www.glytoucan.org/glyco/owl/glytoucan#>
#
# SELECT DISTINCT ?Saccharide ?p ?o ?gp ?go ?p2 ?o2 ?p3 ?o3
# FROM <http://rdf.glytoucan.org/core>
# FROM <http://rdf.glytoucan.org/sequence/wurcs>
# WHERE {
# #    ?Saccharide glytoucan:has_primary_id ?PrimaryID.
#     ?Saccharide glytoucan:has_primary_id 'G00009BX'.
#     ?Saccharide glycan:has_glycosequence ?GlycoSequence .
# #    ?GlycoSequence glycan:has_sequence ?Sequence .
#     ?GlycoSequence ?gp ?go .
# #    ?Saccharide glycan:has_motif ?motif .
#     ?Saccharide ?p ?o.
#     ?go ?p2 ?o2 .
#     ?o2 ?p3 ?o3 .
# }
# order by ?p, ?o, ?gp, ?go, ?p2, ?o2, ?p3, ?o3
# #limit 100
