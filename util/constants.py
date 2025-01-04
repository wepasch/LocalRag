import os
import re

KEY_ID: str = 'id'
KEY_IDS: str = 'ids'
KEY_PMID: str = 'pmid'
KEY_NAME: str = 'name'
KEY_STATUS: str = 'status'
KEY_PAGE: str = 'page'
KEY_SRC: str = 'source'
KEY_QUERY: str = 'query'
KEY_HREF: str = 'href'
KEY_YEAR: str = 'year'
# also Pubmed keys
KEY_TITLE: str = 'Title'
KEY_PUB_DATE: str = 'PubDate'
KEY_DOI: str = 'DOI'

PATH_DATA: str = os.path.join(os.getcwd(), 'data')
PATH_LIBRARY: str = os.path.join(PATH_DATA, 'library')
PATH_LIBRARY_PDF: str = os.path.join(PATH_LIBRARY, 'pdf')
PATH_RECORD: str = os.path.join(PATH_DATA, 'record')
PATH_DOC_RECORD: str = os.path.join(PATH_RECORD, 'doc_record.json')

TEMP_URL_PMC_ID: str = 'https://pmc.ncbi.nlm.nih.gov/articles/pmid/{}/'
TEMP_URL_PM_ID: str = 'https://pubmed.ncbi.nlm.nih.gov/{}/'
URL_DOI: str = 'https://doi.org/'

RE_YEAR: re = re.compile(r'\b(?:19|20)\d{2}\b')

BREAK_TAB: str = '\n\t'
