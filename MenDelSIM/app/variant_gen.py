from GUD.ORM import Gene
import json
from . import establish_GUD_session
session = establish_GUD_session()

def get_genes(name: str) -> list:
    gene = Gene()
    genes = gene.select_by_name(session, name, True)
    response = []
    if len(genes) > 0:
        for g in genes:
            response.append({'uid': int(g.qualifiers["uid"]),
                            'accession' str(g.qualifiers["name"]):,
                            'start': int(g.start),
                            'end': int(g.end)})
    return response