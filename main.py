import os

import pandas as pd
import sqlite3
import bioservices
import logging

logging.basicConfig(level=logging.DEBUG)
pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', 100)


def import_mappings_to_sqlite():
    # create database with interactome and ID mappings in SQL
    os.remove('go-interactome.db')

    tsv_map_f = 'data/idmapping_selected.tab.example'
    tsv_interactome_f = 'data/HuRI.tsv'

    newdb = sqlite3.connect('go-interactome.db')

    # populate interactome database
    for df in pd.read_csv(tsv_interactome_f, sep='\t', chunksize=20, iterator=True, encoding='utf-8',
                          names=['From', 'To']):
        df.to_sql('interactome', newdb, index=False, if_exists='append')

    # populate ID mapping database from file
    hdr = ['UniProtKB_AC ', 'UniProtKB_ID ', 'GeneID (EntrezGene) ', 'RefSeq ', 'GI ', 'PDB ', 'GO ', 'UniRef100 ',
           'UniRef90 ', 'UniRef50 ', 'UniParc ', 'PIR ', 'NCBI_taxon ', 'MIM ', 'UniGene ', 'PubMed ', 'EMBL ',
           'EMBL_CDS ', 'Ensembl ', 'Ensembl_TRS ', 'Ensembl_PRO ', 'Additional PubMed']
    for df in pd.read_csv(tsv_map_f, sep='\t', chunksize=20, iterator=True, encoding='utf-8', names=hdr):
        df = df.rename(columns={c: c.replace(' ', '') for c in df.columns})  # Remove spaces from columns

        df.to_sql('mapping', newdb, index=False, if_exists='append')
    newdb.close()


if __name__ == '__main__':

    # load interactions from HuRI database in Ensembl format
    huri = pd.read_csv('data/huri.tsv', sep='\t', header=0)
    genes = set(huri.values.ravel())

    print(genes)
    print(len(genes))

    import_mappings_to_sqlite()
    db = sqlite3.connect('file:go-interactome.db?mode=rw', uri=True)
    cur = db.cursor()

    ensembl_id = "ENSG00000274474"
    cur.execute("SELECT UniProtKB_AC FROM mapping WHERE Ensembl LIKE :ensmbl_id;",
                {'ensmbl_id': b"%" + ensembl_id.encode('utf-8') + b"%"})
    rows = list(cur)
    print(rows)

    s = bioservices.QuickGO()
    uniprotkb_id = rows[0][0]
    current_page = 1
    total_pages = 1
    while current_page <= total_pages:
        annotations = s.Annotation(geneProductId=uniprotkb_id, page=current_page)
        df = pd.DataFrame(annotations['results'])
        for col in ['withFrom', 'extensions', 'targetSets']:
            df.loc[:, col] = df[col].apply(lambda r: str(r) if r is not None else None)
        df.loc[:, 'Ensembl'] = [ensembl_id] * len(df)
        print(df)
        df.to_sql('annotations', db, index=False, if_exists='append')
        total_pages = annotations['pageInfo']['total']
        print(f"Page {current_page}/{total_pages}")
        current_page += 1

    db.close()
