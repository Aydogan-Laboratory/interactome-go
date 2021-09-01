import gzip
import logging
import os
import sqlite3
import sys

import pandas as pd
from Bio.UniProt.GOA import gafiterator

from downloads import download_files

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
    cur = newdb.cursor()

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
                          names=hdr, usecols=['UniProtKB_AC', 'Ensembl'],
                          dtype={'UniProtKB_AC': 'str', 'Ensembl': 'str'}):
        print(df.size)
        df = df.rename(columns={c: c.replace(' ', '') for c in df.columns})  # Remove spaces from columns

        df.to_sql('mapping', newdb, index=False, if_exists='append')

    cur.execute('''CREATE INDEX "map_index" ON "mapping" (
                    "Ensembl"	ASC,
                    "UniProtKB_AC"
                );''')
    newdb.commit()

    newdb.close()


def import_hgnc_to_sqlite(db_sqlite_file):
    hgnc_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/hgnc_complete_set.txt'
    db = sqlite3.connect(f"file:{db_sqlite_file}", uri=True)
    cur = db.cursor()
    cur.execute('DROP TABLE IF EXISTS hgnc;')
    db.commit()

    for df in pd.read_csv(hgnc_file, sep='\t', encoding='utf-8',
                          chunksize=1e6, iterator=True,
                          usecols=['hgnc_id',
                                   'symbol', 'name',
                                   'alias_symbol', 'alias_name',
                                   'locus_group', 'locus_type',
                                   'entrez_id', 'ensembl_gene_id', 'uniprot_ids']):
        print(df.size)
        df = df.rename(columns={c: c.replace(' ', '') for c in df.columns})  # Remove spaces from columns

        df.to_sql('hgnc', db, index=False, if_exists='append')
    cur.execute('''CREATE INDEX "hgnc_index" ON "hgnc" (
                    "symbol" ASC,
                    "evidenceCode" ASC,
                    "name" ASC
                );''')
    db.commit()


def import_annotations_from_gaf(db_sqlite_file):
    db = sqlite3.connect(f"file:{db_sqlite_file}", uri=True)
    cur = db.cursor()

    cur.execute("DROP TABLE IF EXISTS annotations;")
    cur.execute('''CREATE TABLE "annotations" (
                    "db"	TEXT,
                    "id"	TEXT,
                    "objectSymbol"	TEXT,
                    -- "qualifier"	TEXT,
                    "goId"	TEXT,
                    -- "reference "	TEXT,
                    "evidenceCode"	TEXT,
                    -- "withFrom"	TEXT,
                    "aspect"	TEXT,
                    "objectName"	INTEGER,
                    -- "objectSynonym"	TEXT,
                    "objectType"	TEXT,
                    -- "taxon"	TEXT,
                    "date"	TEXT,
                    "assignedBy"	TEXT,
                    "annotationExtension"	TEXT,
                    "geneProductFromId"	TEXT
                );''')
    db.commit()

    with gzip.open('data/goa_human.gaf.gz', 'rt') as handle:
        it = gafiterator(handle)
        # cur.executemany("""INSERT INTO gaf_annotations (db, id, objectSymbol, qualifier, goId,
        #                  evidenceCode, withFrom, aspect, objectName, objectSynonym,
        #                  objectType, taxon, date, assignedBy, annotationExtension, geneProductFromId) VALUES (
        #                  :DB,:DB_Object_ID,:DB_Object_Symbol,:Qualifier,:GO_ID,
        #                  :Evidence,:With,:Aspect,:DB_Object_Name,:Synonym,
        #                  :DB_Object_Type,:Taxon_ID,:Date,:Assigned_By,:Annotation_Extension,
        #                  :Gene_Product_Form_ID)""", it)
        cur.executemany("""INSERT INTO annotations (db, id, objectSymbol, goId, evidenceCode, 
                         aspect, objectName, objectType, date, assignedBy, annotationExtension, 
                         geneProductFromId) VALUES (
                         :DB,:DB_Object_ID,:DB_Object_Symbol,:GO_ID,
                         :Evidence,:Aspect,:DB_Object_Name,
                         :DB_Object_Type,:Date,:Assigned_By,:Annotation_Extension,
                         :Gene_Product_Form_ID)""", it)
    db.commit()

    cur.execute('''CREATE INDEX "gaf_index" ON "annotations" (
                    "goId" ASC,
                    "goEvidence" ASC
                );''')
    db.commit()
    db.close()


if __name__ == '__main__':
    huri_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/HI-union.tsv'
    map_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/idmapping_selected.tab.gz'
    raw_map_file = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/idmapping.dat.gz'
    db_sqlite = '/media/lab/Data/Fabio/Dev/Python-InteractomeGO/data/go-interactome.db'
    download_files()

    # load interactions from HuRI database in Ensembl format
    import_mappings_to_sqlite(db_sqlite, huri_file, map_file)
    import_hgnc_to_sqlite(db_sqlite)
    import_annotations_from_gaf(db_sqlite)
