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


def download_files():
    go_ontology = 'http://geneontology.org/ontology/go.obo'
    download_url(go_ontology, os.path.abspath('data/go.obo'))

    hgnc = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
    download_url(hgnc, os.path.abspath('data/hgnc_complete_set.txt'))

    homosapiens_annotations = 'http://geneontology.org/gene-associations/goa_human.gaf.gz'
    download_url(homosapiens_annotations, os.path.abspath('data/goa_human.gaf.gz'))

    drosophila_annotations = 'http://current.geneontology.org/annotations/fb.gaf.gz'
    download_url(drosophila_annotations, os.path.abspath('data/fb.gaf.gz'))
