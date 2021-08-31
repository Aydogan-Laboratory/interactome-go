import sqlite3

import owlready2 as owl
import pandas as pd
import py4cytoscape as p4c
from goatools.godag.go_tasks import get_go2parents, get_go2children
from goatools.obo_parser import GODag
from owlready2 import ThingClass

from downloads import download_files


def get_children(ontology: GODag, node: str, levels=1, all_children=None, optional_relationships=None) -> set:
    all_terms = set()
    if levels >= 1:
        if all_children is None:
            if optional_relationships is None:
                optional_relationships = {'is_a'}
            all_children = get_go2children(ontology, optional_relationships)
        if node in all_children:
            parents = all_children[node]
            for term in parents:
                all_terms = all_terms.union(
                    parents,
                    get_children(ontology, term,
                                 levels=levels - 1,
                                 all_children=all_children,
                                 optional_relationships=optional_relationships)
                )

    return all_terms


def get_parents(ontology: GODag, node: str, levels=1, all_parents=None, optional_relationships=None) -> set:
    all_terms = set()
    if levels >= 1:
        if all_parents is None:
            if optional_relationships is None:
                optional_relationships = {'is_a'}
            all_parents = get_go2parents(ontology, optional_relationships)
        if node in all_parents:
            parents = all_parents[node]
            for term in parents:
                all_terms = all_terms.union(
                    parents,
                    get_parents(ontology, term,
                                levels=levels - 1,
                                all_parents=all_parents,
                                optional_relationships=optional_relationships)

                )

    return all_terms


def obo_test():
    g = GODag("data/go-basic.obo", optional_attrs={'relationship'})
    stk_activity = g.query_term('GO:0004674', verbose=False)

    print(stk_activity.get_all_parents())
    g.draw_lineage(
        [stk_activity],
        # dpi=opts.dpi,
        # engine=opts.engine,
        # gml=opts.gml,
        # output=opts.output,
        draw_parents=True,
        draw_children=False,
    )


def carlton_cell_cycle_tags():
    onto = owl.get_ontology("file://data/go.owl")
    onto.load()
    # for c in onto.classes(): print(c.name, c.label,)
    carlton_go_tags = [
        onto.search_one(id='GO:0051301'),  # Cell  Division
        onto.search_one(id='GO:0010458'),  # Exit  from  Mitosis
        onto.search_one(id='GO:0000278'),  # Mitotic  Cell  Cycle
        onto.search_one(id='GO:0000281'),  # Mitotic  Cytokinesis
        onto.search_one(id='GO:1902412'),  # Regulation of Mitotic Cytokinesis
    ]
    all_terms = carlton_go_tags
    # cell_cycle = onto.search_one(id='GO:0007049')  # cell cycle
    # stk_activity = onto.search_one(id='GO:0004674')  # protein serine/threonine kinase activity
    # cyclin_stk_activity = onto.search_one(id='GO:0004693')  # cyclin-dependent protein serine/threonine kinase activity
    # print(stk_activity.is_a)
    # print(stk_activity.is_part)
    # exit()
    df_go = pd.DataFrame({'goId': [c.name for c in all_terms], 'name': [c.label[0] for c in all_terms]})
    with pd.ExcelWriter('cytoscape_ontology.xlsx', mode='w') as writer:
        df_go.sort_values(by='goId').to_excel(writer, sheet_name='Ontology')
    return [o.name.replace('_', ':') for o in all_terms if type(o) == ThingClass]


def my_cell_cycle_tags():
    g = GODag("data/go-basic.obo", optional_attrs={'relationship'})
    cell_cycle = g.query_term('GO:0007049', verbose=False)

    optional_relationships = {'is_a', 'part_of', 'regulates', 'negatively_regulates', 'positively_regulates'}
    cc_children = get_children(g, cell_cycle.item_id, levels=100, optional_relationships=optional_relationships)

    print(cell_cycle)
    terms = list()
    for p in cc_children:
        go_node = g.query_term(p, verbose=False)
        terms.append({'goId': go_node.item_id, 'name': go_node.name, 'ns': go_node.namespace})
        print(go_node)

    df_go = pd.DataFrame(terms)
    with pd.ExcelWriter('cytoscape_ontology.xlsx', mode='w') as writer:
        df_go.sort_values(by='goId').to_excel(writer, sheet_name='Ontology')
    return df_go['goId'].values


def interactors(db_sqlite):
    go_tags = my_cell_cycle_tags()
    print(go_tags)
    print(len(go_tags))
    db = sqlite3.connect(f"file:{db_sqlite}", uri=True)
    cur = db.cursor()

    cur.execute('DROP TABLE IF EXISTS filtered;')
    cur.execute('CREATE TABLE filtered ("Ensembl" TEXT,  "symbol" TEXT,  "alias_symbol" TEXT,  "name" TEXT);')
    csep_str = "', '"
    fensmbl_qstr = f"""
        INSERT INTO filtered
        SELECT M.Ensembl, H.symbol, H.alias_symbol, H.name
            FROM mapping AS M
            INNER JOIN hgnc H ON H.ensembl_gene_id = M.Ensembl
            WHERE M.Ensembl IN (
                SELECT DISTINCT Ensembl as FilteredGenes
                FROM annotations
                WHERE goId IN ('{csep_str.join(go_tags)}')
                --AND goEvidence IN ('EXP','IDA','IPI','IMP','IGI','IEP', 'HTP','HDA','HMP','HGI','HEP')
        );
        """
    cur.execute(fensmbl_qstr)
    db.commit()

    sql_filtered_genes = f"""
    SELECT DISTINCT I.p1 AS p1e, 
        I.p2 AS p2e, 
        F1.symbol as p1gn, 
        F2.symbol as p2gn
    FROM interactome as I
    INNER JOIN filtered F1 ON F1.Ensembl = I.p1
    INNER JOIN filtered F2 ON F2.Ensembl = I.p2;
    """
    cur.execute(sql_filtered_genes)
    df = pd.DataFrame(cur, columns=['id1', 'id2', 'GeneName1', 'GeneName2'])
    db.close()

    with pd.ExcelWriter('cytoscape_ontology.xlsx', engine="openpyxl", mode='a') as writer:
        df.to_excel(writer, sheet_name='Interactions')

    return df


def filtered_nodes(db_sqlite):
    db = sqlite3.connect(f"file:{db_sqlite}", uri=True)

    df = pd.read_sql("""
    SELECT DISTINCT I1.p1 as id, F1.symbol as name, F1.name as desc FROM interactome as I1
    INNER JOIN filtered F1 ON F1.Ensembl = I1.p1
    UNION 
    SELECT DISTINCT I2.p2 as id, F2.symbol as name, F2.name as desc  FROM interactome as I2
    INNER JOIN filtered F2 ON F2.Ensembl = I2.p2;
    """, db)
    print('CDK1' in df['name'])
    return df


if __name__ == '__main__':
    download_files()
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'

    my_cell_cycle_tags()

    print(dir(p4c))
    print(p4c.cytoscape_ping())
    print(p4c.cytoscape_version_info())
    p4c.close_session(save_before_closing=False)

    interactions = interactors(db_sqlite)
    nodes = filtered_nodes(db_sqlite)
    edges = pd.DataFrame(data={
        'source':      interactions['id1'],
        'target':      interactions['id2'],
        'interaction': ["interacts"] * len(interactions),
        'weight':      [1.0] * len(interactions)})

    print(interactions)
    print(nodes)
    print(edges)
    p4c.create_network_from_data_frames(nodes, edges,
                                        title="cell cycle network",
                                        collection="HuRI")
    p4c.set_visual_style('Marquee')
    # p4c.load_table_data(all_genes, data_key_column='Ensembl', table='node', table_key_column='name')
    # p4c.close_session(save_before_closing=False)
