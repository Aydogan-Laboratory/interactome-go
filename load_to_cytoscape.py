import sqlite3
import pandas as pd
import py4cytoscape as p4c
import owlready2 as owl
from owlready2 import ThingClass


def owl_tags():
    onto = owl.get_ontology("file://data/go.owl")
    onto.load()
    # flt = onto.search(comment='*mitochon*')
    flt = onto.search(label='*cell*cycle*')
    return [o.name.replace('_', ':') for o in flt if type(o) == ThingClass]


def interactors(db_sqlite):
    go_tags = owl_tags()
    print(go_tags)
    db = sqlite3.connect(f"file:{db_sqlite}", uri=True)
    cur = db.cursor()

    cur.execute('DROP TABLE IF EXISTS temp.filtered;')
    cur.execute(
        'CREATE VIRTUAL TABLE temp.filtered using fts3("UniProtKB_AC" TEXT,  "Ensembl" TEXT,  "geneName" TEXT);')
    csep_str = "', '"
    fensmbl_qstr = f"""
        INSERT INTO temp.filtered
        SELECT M.UniProtKB_AC, M.Ensembl, RM.ID 
            FROM mapping AS M
            INNER JOIN rawmap RM ON RM.UniProtKB_AC = M.UniProtKB_AC
            WHERE M.Ensembl IN (
                SELECT DISTINCT Ensembl as FilteredGenes
                FROM annotations 
                WHERE goEvidence IN ('EXP','IDA','IPI','IMP','IGI','IEP', 'HTP','HDA','HMP','HGI','HEP')
                AND goId IN ('{csep_str.join(go_tags)}')
        );
        """
    cur.execute(fensmbl_qstr)

    sql_filtered_genes = f"""
    SELECT I.p1 AS p1e, 
        I.p2 AS p2e, 
        F1.geneName as p1gn, 
        F2.geneName as p2gn
        --F1.UniProtKB_AC as p1u, 
        --F2.UniProtKB_AC as p2u, 
    FROM interactome as I
    INNER JOIN filtered F1 ON F1.Ensembl = I.p1
    INNER JOIN filtered F2 ON F2.Ensembl = I.p2;
    """
    cur.execute(sql_filtered_genes)
    df = pd.DataFrame(cur, columns=['id1', 'id2', 'GeneName1', 'GeneName2'])
    db.close()
    return df


if __name__ == '__main__':
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'

    print(dir(p4c))
    print(p4c.cytoscape_ping())
    print(p4c.cytoscape_version_info())

    interactions = interactors(db_sqlite)
    nodes = interactions[['id1', 'GeneName1']].rename(columns={'id1': 'id', 'GeneName1': 'name'}).append(
        interactions[['id2', 'GeneName2']].rename(columns={'id2': 'id', 'GeneName2': 'name'})).drop_duplicates()

    edges = pd.DataFrame(data={
        'source':      interactions['id1'],
        'target':      interactions['id2'],
        'interaction': ["interacts"] * len(interactions),
        'weight':      [1.0] * len(interactions)})

    # print(interactions)
    # print(nodes)
    # print(edges)
    p4c.create_network_from_data_frames(nodes, edges,
                                        title="my first network",
                                        collection="HuRI")
    # p4c.load_table_data(all_genes, data_key_column='Ensembl', table='node', table_key_column='name')
    # p4c.close_session(save_before_closing=False)
