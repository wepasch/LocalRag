import logging

from Bio import Entrez

from util.basics import add_if_available, extract_first_year
from util.config import RET_MAX_PM, USED_MAIL, ADD_INFO_KEYS
from util.constants import KEY_PMID, KEY_TITLE, KEY_PUB_DATE, KEY_YEAR

logger = logging.getLogger(__name__)


def articles_for(query: str) -> list[dict]:
    full_query: str = f'{query} AND free full text[filter]'
    Entrez.email = USED_MAIL
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=RET_MAX_PM,
                            retmode='xml',
                            term=full_query)
    results = Entrez.read(handle)
    pmids: list[str] = [str(i) for i in results['IdList']]
    return [{KEY_PMID: i} | __get_add_info(i) for i in pmids]


def __get_add_info(pmid: str) -> dict[str, str]:
    info: dict[str, str] = {}
    try:
        with Entrez.esummary(db='pubmed', id=pmid, retmode='xml') as h:
            records = Entrez.read(h)
            try:
                record: dict = records[0]
            except IndexError:
                record: dict = {}
            for key in ADD_INFO_KEYS:
                add_if_available(record, key, info)
            if KEY_PUB_DATE in info:
                year: str | None = extract_first_year(info[KEY_PUB_DATE])
                if year:
                    info[KEY_YEAR] = year

    except Exception as e:
        logger.warning(f'Can not get pubmed data for {pmid} [{e}]')
    return info
