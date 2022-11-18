import sqlite3

import pandas as pd
import py4cytoscape as p4c
from goatools.godag.go_tasks import get_go2parents, get_go2children
from goatools.obo_parser import GODag

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


def carlton_cell_cycle_tags():
    g = GODag("data/go.obo", optional_attrs={'relationship'})

    carlton_go_tags = [
        g.query_term('GO:0051301', verbose=False),  # Cell  Division
        g.query_term('GO:0010458', verbose=False),  # Exit  from  Mitosis
        g.query_term('GO:0000278', verbose=False),  # Mitotic  Cell  Cycle
        g.query_term('GO:0000281', verbose=False),  # Mitotic  Cytokinesis
        g.query_term('GO:1902412', verbose=False),  # Regulation of Mitotic Cytokinesis
    ]

    # search and add children of the terms to the list
    terms = list()
    optional_relationships = {'is_a', 'part_of', 'regulates', 'negatively_regulates', 'positively_regulates'}
    for cc in carlton_go_tags:
        terms.append({'goId': cc.item_id, 'name': cc.name, 'ns': cc.namespace})
        cc_children = get_children(g, cc.item_id, levels=100, optional_relationships=optional_relationships)
        for c in cc_children:
            go_node = g.query_term(c, verbose=False)
            terms.append({'goId': go_node.item_id, 'name': go_node.name, 'ns': go_node.namespace})

    df_go = pd.DataFrame(terms)
    with pd.ExcelWriter('cytoscape_ontology.xlsx', mode='w') as writer:
        df_go.sort_values(by='goId').to_excel(writer, sheet_name='CelCycle-GO')
    return df_go['goId'].values


def carlton_mitochondria_tags():
    g = GODag("data/go.obo", optional_attrs={'relationship'})

    cc_terms = [
        g.query_term('GO:0005743', verbose=False),  # Mitochondrial inner membrane
        g.query_term('GO:0005758', verbose=False),  # Mitochondrial intermembrane space
        g.query_term('GO:0005759', verbose=False),  # Mitochondrial matrix
        g.query_term('GO:0005741', verbose=False),  # Mitochondrial outer membrane
        g.query_term('GO:0005739', verbose=False),  # Mitochondrion
    ]

    # search and add children of the terms to the list
    terms = list()
    optional_relationships = {'is_a', 'part_of', 'regulates', 'negatively_regulates', 'positively_regulates'}
    for cc in cc_terms:
        terms.append({'goId': cc.item_id, 'name': cc.name, 'ns': cc.namespace})
        cc_children = get_children(g, cc.item_id, levels=100, optional_relationships=optional_relationships)
        for c in cc_children:
            go_node = g.query_term(c, verbose=False)
            terms.append({'goId': go_node.item_id, 'name': go_node.name, 'ns': go_node.namespace})

    df_go = pd.DataFrame(terms)
    with pd.ExcelWriter('cytoscape_ontology.xlsx', mode='w') as writer:
        df_go.sort_values(by='goId').to_excel(writer, sheet_name='Mitochondria-GO')
    return df_go['goId'].values


def carlton_er_tags():
    g = GODag("data/go.obo", optional_attrs={'relationship'})

    er_terms = [
        g.query_term('GO:0005783', verbose=False),  # Endoplasmic reticulum
        g.query_term('GO:0005788', verbose=False),  # Endoplasmic reticulum lumen
        g.query_term('GO:0005789', verbose=False),  # Endoplasmic reticulum membrane
        g.query_term('GO:0005793', verbose=False),  # Endoplasmic reticulum-Golgi intermediate compartment membrane
        g.query_term('GO:0030176', verbose=False),  # Integral component of endoplasmic reticulum membrane
    ]

    # search and add children of the terms to the list
    terms = list()
    optional_relationships = {'is_a', 'part_of', 'regulates', 'negatively_regulates', 'positively_regulates'}
    for cc in er_terms:
        terms.append({'goId': cc.item_id, 'name': cc.name, 'ns': cc.namespace})
        cc_children = get_children(g, cc.item_id, levels=100, optional_relationships=optional_relationships)
        for c in cc_children:
            go_node = g.query_term(c, verbose=False)
            terms.append({'goId': go_node.item_id, 'name': go_node.name, 'ns': go_node.namespace})

    df_go = pd.DataFrame(terms)
    with pd.ExcelWriter('cytoscape_ontology.xlsx', mode='w') as writer:
        df_go.sort_values(by='goId').to_excel(writer, sheet_name='EndopasmicReticulum-GO')
    return df_go['goId'].values


def my_cell_cycle_tags():
    g = GODag("data/go.obo", optional_attrs={'relationship'})
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
    db = sqlite3.connect(f"file:{db_sqlite}", uri=True)
    cur = db.cursor()
    csep_str = "', '"

    go_tags = set(carlton_cell_cycle_tags()).union(my_cell_cycle_tags())
    print(go_tags)
    print(len(go_tags))

    print("Populating table of filtered gene products.")
    cur.execute('DROP TABLE IF EXISTS filtered;')
    cur.execute('CREATE TABLE filtered ("Ensembl" TEXT, "UniProtKB" TEXT, "symbol" TEXT, "alias_symbol" TEXT, '
                '"name" TEXT, "related_to" TEXT);')
    cur.execute('''CREATE INDEX "filtered_index" ON "filtered" (
                    "Ensembl" ASC,
                    "symbol" ASC
                );''')
    db.commit()

    fensmbl_qstr = f"""
        INSERT INTO filtered
        SELECT M.Ensembl, H.uniprot_ids, H.symbol, H.alias_symbol, H.name, 'cell_cycle'
            FROM mapping AS M
            INNER JOIN hgnc H ON H.ensembl_gene_id = M.Ensembl
            WHERE H.uniprot_ids IN (
                SELECT DISTINCT a.id as FilteredGenes
                FROM annotations a
                WHERE a.db='UniProtKB' 
                AND goId IN ('{csep_str.join(go_tags)}')
                --AND evidenceCode IN ('EXP','IDA','IPI','IMP','IGI','IEP', 'HTP','HDA','HMP','HGI','HEP')
        );
        """
    cur.execute(fensmbl_qstr)
    db.commit()
    print("Table filled with cell cycle entries.")

    organelle = 'mitochondria'
    organelle_go_tags = carlton_mitochondria_tags()
    in_organelle_qstr = f"""
        INSERT INTO filtered
        SELECT M.Ensembl, H.uniprot_ids, H.symbol, H.alias_symbol, H.name, '{organelle}'
            FROM mapping AS M
            INNER JOIN hgnc H ON H.ensembl_gene_id = M.Ensembl
            WHERE H.uniprot_ids IN (
                SELECT DISTINCT a.id as FilteredGenes
                FROM annotations a
                WHERE a.db='UniProtKB' 
                AND goId IN ('{csep_str.join(organelle_go_tags)}')
                --AND evidenceCode IN ('EXP','IDA','IPI','IMP','IGI','IEP', 'HTP','HDA','HMP','HGI','HEP')
        );
    """
    print(in_organelle_qstr)
    cur.execute(in_organelle_qstr)
    db.commit()
    print("Table filled with organelle entries.")

    sql_filtered_genes = f"""
    SELECT DISTINCT 
        I.p1 AS p1e, 
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
    cur = db.cursor()
    organelle = 'mitochondria'

    print("Creating table of network nodes.")
    cur.execute('DROP TABLE IF EXISTS nodes;')
    cur.execute('CREATE TABLE nodes ("id" TEXT, "name" TEXT, "desc" TEXT, '
                '"cell_cycle" NUMERIC DEFAULT 0, '
                f'"{organelle}" NUMERIC DEFAULT 0);')
    cur.execute(f'''CREATE INDEX "node_index" ON "nodes" (
                    "id" ASC,
                    "cell_cycle" ASC,
                    "{organelle}" ASC
                );''')
    db.commit()

    sql_nodes_insert = f"""
    INSERT INTO nodes
    SELECT DISTINCT id, name, desc, 0, 0 FROM 
    (
        SELECT DISTINCT I1.p1 as id, F1.symbol as name, F1.name as desc FROM interactome as I1
        INNER JOIN filtered F1 ON F1.Ensembl = I1.p1
        UNION 
        SELECT DISTINCT I2.p2 as id, F2.symbol as name, F2.name as desc FROM interactome as I2
        INNER JOIN filtered F2 ON F2.Ensembl = I2.p2
    );
    """
    cur.execute(sql_nodes_insert)

    sql_nodes_update_1 = f"""
    UPDATE nodes SET cell_cycle=1
    WHERE id IN (
        SELECT F.Ensembl FROM filtered as F
        WHERE F.related_to='cell_cycle'
    );
    """
    cur.execute(sql_nodes_update_1)

    sql_nodes_update_2 = f"""
    UPDATE nodes SET {organelle}=1
    WHERE id IN (
        SELECT F.Ensembl FROM filtered as F
        WHERE F.related_to='{organelle}'
    );
    """
    cur.execute(sql_nodes_update_2)
    db.commit()

    df = pd.read_sql("SELECT DISTINCT * FROM nodes;", db)
    db.close()
    print('CDK1' in df['name'])
    return df


if __name__ == '__main__':
    download_files()
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'

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
                                        title="cell cycle & endoplasmic_reticulum",
                                        collection="HuRI")
    p4c.set_visual_style('Marquee')
