import os
import urllib.request

from tqdm import tqdm


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    if os.path.exists(output_path):
        return
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def download_to_data(url):
    fname = os.path.basename(url)
    download_url(url, os.path.abspath(os.path.join('data', fname)))


def download_files():
    # ontologies
    download_to_data('http://geneontology.org/ontology/go.owl')
    # ontologies from amigo
    download_to_data('http://purl.obolibrary.org/obo/cl/cl-simple.owl')
    download_to_data('http://purl.obolibrary.org/obo/wbbt.owl')
    download_to_data('http://purl.obolibrary.org/obo/uberon/basic.owl')
    download_to_data('http://purl.obolibrary.org/obo/po/imports/ncbitaxon_import.owl')
    download_to_data('http://purl.obolibrary.org/obo/eco/eco-basic.owl')
    download_to_data('http://purl.obolibrary.org/obo/chebi.owl')
    download_to_data('http://purl.obolibrary.org/obo/po.owl')
    download_to_data('http://purl.obolibrary.org/obo/go/extensions/go-gaf.owl')
    download_to_data('http://purl.obolibrary.org/obo/po/imports/ro_import.owl')
    download_to_data('http://purl.obolibrary.org/obo/go/extensions/gorel.owl')
    download_to_data('http://purl.obolibrary.org/obo/ncbitaxon/subsets/taxslim.owl')
    download_to_data('http://purl.obolibrary.org/obo/pato.owl')
    download_to_data('http://purl.obolibrary.org/obo/go/extensions/go-modules-annotations.owl')
    download_to_data('http://purl.obolibrary.org/obo/go/extensions/go-taxon-subsets.owl')

    download_to_data('http://geneontology.org/ontology/go.obo')
    download_to_data('http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt')
    download_to_data('http://www.interactome-atlas.org/data/HI-union.tsv')
    download_to_data(
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz')
    download_to_data('http://geneontology.org/gene-associations/goa_human.gaf.gz')  # annotations for homo sapiens
    download_to_data('http://current.geneontology.org/annotations/fb.gaf.gz')  # annotations for drosophila

    # annotations from amigo
    download_to_data('http://skyhook.berkeleybop.org/release/products/annotations/paint_other.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/aspgd.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/cgd.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/dictybase.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/ecocyc.gaf.gz')
    # download_to_data('http://skyhook.berkeleybop.org/release/annotations/fb.gaf.gz')
    # download_to_data('http://skyhook.berkeleybop.org/release/annotations/goa_human.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/goa_human_complex.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/goa_human_rna.gaf.gz')
    download_to_data('http://skyhook.berkeleybop.org/release/annotations/goa_uniprot_all_noiea.gaf.gz')

    # annotations from StringDB
    # download_to_data('https://stringdb-static.org/download/database.schema.v11.5.pdf')
    download_to_data('https://stringdb-static.org/download/items_schema.v11.5.sql.gz')
    download_to_data('https://stringdb-static.org/download/network_schema.v11.5.sql.gz')
    download_to_data('https://stringdb-static.org/download/evidence_schema.v11.5.sql.gz')
    download_to_data('https://stringdb.meringlab.org/download/homology_schema.v11.5.sql.gz')
