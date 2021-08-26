import sqlite3
import pandas as pd
import py4cytoscape as p4c
import owlready2 as owl
from owlready2 import ThingClass


def owl_tags():
    onto = owl.get_ontology("file://data/go.owl")
    onto.load()
    # mito = onto.search(comment='*mitochon*')
    cc = onto.search(label='*cell*cycle*')
    return [o.name.replace('_', ':') for o in cc if type(o) == ThingClass]


def interactors(db_sqlite):
    go_tags = owl_tags()
    db = sqlite3.connect(f"file:{db_sqlite}", uri=True)
    cur = db.cursor()

    csep_str = "', '"
    fensmbl_qstr = f"""
            SELECT DISTINCT Ensembl as FilteredGenes
            FROM annotations 
            WHERE goEvidence IN ('EXP','IDA','IPI','IMP','IGI','IEP', 'HTP','HDA','HMP','HGI','HEP')
            AND goId IN ('{csep_str.join(go_tags)}');
        """

    cur.execute(fensmbl_qstr)
    filt_ensembl = [r[0] for r in cur]

    sql_filtered_genes = f"""
    SELECT p1, p2 FROM interactome 
    WHERE p1 IN ('{csep_str.join(filt_ensembl)}') OR p2 IN ('{csep_str.join(filt_ensembl)}')
    """
    cur.execute(sql_filtered_genes)
    out = list(cur)
    db.close()
    return out


if __name__ == '__main__':
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'

    print(dir(p4c))
    print(p4c.cytoscape_ping())
    print(p4c.cytoscape_version_info())

    interactions = interactors(db_sqlite)
    src = [i[0] for i in interactions]
    trg = [i[1] for i in interactions]
    nodes = pd.DataFrame(data={'id': sorted(set(src + trg))})
    edges = pd.DataFrame(data={
        'source':      src,
        'target':      trg,
        'interaction': ["interacts"] * len(interactions),
        'weight':      [1.0] * len(interactions)})

    p4c.create_network_from_data_frames(nodes, edges, title="my first network", collection="HuRI")
    # p4c.load_table_data(all_genes, data_key_column='Ensembl', table='node', table_key_column='name')
    p4c.close_session()
