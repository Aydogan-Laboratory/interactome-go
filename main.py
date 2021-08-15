import os
import sys
import csv
import gzip
import time

import pandas as pd
import sqlite3
import bioservices
import logging

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
log = logging.getLogger('GO-db')
pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 100)


def import_mappings_to_sqlite(db_sqlite_file, huri_tsv_file, mapping_tsv_file):
    # create database with interactome and ID mappings in SQL
    if os.path.exists(db_sqlite_file):
        os.remove(db_sqlite_file)

    newdb = sqlite3.connect(db_sqlite_file)

    # populate interactome database
    for df in pd.read_csv(huri_tsv_file, sep='\t', encoding='utf-8',
                          chunksize=1e6, iterator=True,
                          names=['p1', 'p2']):
        df.to_sql('interactome', newdb, index=False, if_exists='append')

    # populate ID mapping database from file
    hdr = ['UniProtKB_AC', 'UniProtKB_ID', 'GeneID (EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100',
           'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI_taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
           'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed']
    for df in pd.read_csv(mapping_tsv_file, sep='\t', encoding='utf-8',
                          chunksize=1e6, iterator=True,
                          names=hdr, usecols=['UniProtKB_AC', 'Ensembl', ''],
                          dtype={'UniProtKB_AC': 'str', 'Ensembl': 'str'}):
        print(df.size)
        df = df.rename(columns={c: c.replace(' ', '') for c in df.columns})  # Remove spaces from columns

        df.to_sql('mapping', newdb, index=False, if_exists='append')
    newdb.close()


def import_raw_mappings_to_sqlite(db_sqlite_file, raw_mapping_tsv_file):
    con = sqlite3.connect(db_sqlite_file)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS rawmap;")
    cur.execute('''CREATE TABLE rawmap (
                    "UniProtKB_AC"	TEXT,
                    "ID_type"	TEXT,
                    "ID"	TEXT
                );''')
    a_file = gzip.open(raw_mapping_tsv_file, 'rt')
    rows = csv.reader(a_file, delimiter='\t')
    rows_to_insert = list()
    for r in rows:
        if 'Gene_Name' in r:
            log.debug(f"Inserting {r} into DB.")
            rows_to_insert.append(r)
    cur.executemany("INSERT INTO rawmap VALUES (?, ?, ?)", rows_to_insert)

    con.commit()
    con.close()


if __name__ == '__main__':
    huri_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/HuRI.tsv'
    map_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/idmapping_selected.tab.gz'
    raw_map_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/idmapping.dat.gz'
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'

    # load interactions from HuRI database in Ensembl format
    import_mappings_to_sqlite(db_sqlite, huri_file, map_file)
    import_raw_mappings_to_sqlite(db_sqlite, raw_map_file)
    db = sqlite3.connect(db_sqlite)
    cur = db.cursor()

    cur.execute("SELECT DISTINCT p1 as p FROM interactome UNION SELECT DISTINCT p2 as p FROM interactome;")
    all_genes = [r[0] for r in cur]
    cur.execute("SELECT DISTINCT Ensembl FROM annotations;")
    already_annotated_genes = [r[0] for r in cur]
    cur.execute("SELECT Ensembl, UniProtKB_AC FROM mapping WHERE Ensembl LIKE 'ENSG%';")
    all_maps = {r[0]: r[1] for r in cur}
    for ensembl_id in all_genes:
        if not ensembl_id in all_maps.keys():
            log.warning(f"No results found for Ensembl ID {ensembl_id}.")
            continue
        if ensembl_id in already_annotated_genes:
            log.warning(f"Skipping Ensembl ID {ensembl_id} because it's already annotated.")
            continue
        uniprotkb_id = all_maps[ensembl_id]
        log.info(f"Downloading info of gene with Ensembl ID={ensembl_id}, UniProtKB={uniprotkb_id}.")

        s = bioservices.QuickGO()
        current_page = 1
        total_pages = 1
        while current_page <= total_pages:
            annotations = s.Annotation(geneProductId=uniprotkb_id, page=current_page)
            df = pd.DataFrame(annotations['results'])
            for col in ['withFrom', 'extensions', 'targetSets']:
                if col in df.columns:
                    df.loc[:, col] = df[col].apply(lambda r: str(r) if r is not None else None)
            df.loc[:, 'Ensembl'] = [ensembl_id] * len(df)
            df.to_sql('annotations', db, index=False, if_exists='append')
            total_pages = annotations['pageInfo']['total']
            log.debug(f"Page {current_page}/{total_pages}")
            current_page += 1
            time.sleep(0.4)
    db.commit()
    db.close()
