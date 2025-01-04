import argparse
import logging
import json
import os
from urllib.parse import urljoin

from documents.document_status import DocumentStatus
from documents.pubmed import articles_for

from util.basics import setup_logging, is_pdf, delete_file
from util.constants import PATH_DOC_RECORD, KEY_STATUS, KEY_PMID, TEMP_URL_PMC_ID, KEY_QUERY, TEMP_URL_PM_ID, URL_DOI, \
    PATH_LIBRARY_PDF, BREAK_TAB, KEY_TITLE
from util.web import get_raw, download_pdf, get_ddg_filetype_urls

logger = logging.getLogger(__name__)


def main(query: str = '', download: bool = False, dry_run: bool = False):
    logger.info(f'Documents update:{BREAK_TAB}'
                f'query: {query if query else "<no query>"},{BREAK_TAB}'
                f'download: {download}{BREAK_TAB}'
                f'dry run: {dry_run}')
    if query:
        logger.info(f'Find new files for query: {query}')
        __process_query(query)
    if download:
        logger.info('Try to download new pmids.')
        __attempt_download(dry_run=dry_run)


def __process_query(query: str):
    doc_record: list[dict] = __load_doc_record()
    new_articles_data: list[dict] = __find_new_articles(query, {a[KEY_PMID] for a in doc_record})
    logger.info(f'Identified {len(new_articles_data)} new articles for query \'{query}\'.')
    for data in new_articles_data:
        data.update({KEY_STATUS: DocumentStatus.NEW.value})
        data.update({KEY_QUERY: query})
        doc_record.append(data)
    __write_doc_record(doc_record)


def __attempt_download(dry_run: bool = False) -> None:
    doc_record: list[dict] = __load_doc_record()
    new_articles_data: list[dict] = [a for a in doc_record if a[KEY_STATUS] == DocumentStatus.NEW.value]
    logger.info(f'Attempting to download {len(new_articles_data)} new articles.')
    nof_a: int = len(new_articles_data)
    ida: int
    data: dict
    for ida, data in enumerate(new_articles_data):
        print()
        logger.info(f'Attempt download {data[KEY_PMID]} [{ida + 1}/{nof_a}].')
        if __download_if_not_exist(data):
            data[KEY_STATUS] = DocumentStatus.DOWNLOADED.value
        else:
            data[KEY_STATUS] = DocumentStatus.DOWNLOAD_FAILED.value

    if dry_run:
        logger.info(f'Dry run mode. Do not write into doc record.')
    else:
        __write_doc_record(doc_record)


def __download_if_not_exist(data: dict) -> bool:
    pdf_path: str = f'{PATH_LIBRARY_PDF}/{data[KEY_PMID]}.pdf'
    if os.path.exists(pdf_path):
        logger.info(f'PDF for {data[KEY_PMID]} already exits.')
        return True
    pdf_urls: set[str] = __get_pm_pdf_urls(data)
    pdf_urls.update(__get_ddg_pdf_urls(data))
    logger.info(f'Found {len(pdf_urls)} potential pdf links.')
    if __attempt_download_of_pdfs(pdf_urls, pdf_path):
        return True
    return False

def __get_ddg_pdf_urls(data: dict) -> set[str]:
    return get_ddg_filetype_urls(data[KEY_TITLE], 'pdf')


def __attempt_download_of_pdfs(pds_urls: set[str], pdf_path: str) -> bool:
    for u in pds_urls:
        logger.info(f'Downloading pdf from {u}.')
        if download_pdf(u, pdf_path):
            if is_pdf(pdf_path):
                logger.info('...downloaded OK.')
                return True
            logger.warning('...download corrupt, delete file.')
            delete_file(pdf_path)
    return False

def __get_pm_pdf_urls(data: dict) -> set[str]:
    found_pdf_urls: set[str] = __find_pdf_urls_pmc(data)
    if not found_pdf_urls:
        found_pdf_urls = __find_pdf_url_pm(data)
    return found_pdf_urls

def __find_pdf_url_pm(data: dict) -> set[str]:
    pmid: str = data[KEY_PMID]
    article_url: str = TEMP_URL_PM_ID.format(pmid)
    html = get_raw(article_url)
    found_doi_urls: set[str] = set()
    for a in html.find_all('a', href=True):
        h = a.get('href', '<>')
        if h.startswith(URL_DOI):
            found_doi_urls.add(h)

    found_pdf_urls: set[str] = set()
    found_pdf_candidate_urls: set[str] = set()
    for d in found_doi_urls:
        found_pdf_candidate_urls.update(__find_pdf_text_urls(d))
    for c in found_pdf_candidate_urls:
        found_pdf_urls.update(__find_class_urls(c, 'download'))
    return found_pdf_urls


def __find_with_title_urls(url: str, title_part: str) -> set[str]:
    html = get_raw(url)
    print(html)
    found_download_urls: set[str] = set()
    for a in html.find_all('a', href=True):
        t: str = a.get('title', '<>')
        if title_part.lower() in t.lower():
            found_download_urls.add(a)
    return found_download_urls

def __find_class_urls(url: str, class_name: str) -> set[str]:
    html = get_raw(url)
    found_download_urls: set[str] = set()
    for a in html.find_all('a', class_=class_name):
        found_download_urls.add(a.get('href', '<>'))
    return found_download_urls


def __find_pdf_text_urls(url: str) -> set[str]:
    html = get_raw(url)
    found_pdf_urls: set[str] = set()
    for a in html.find_all('a', href=True):
        t = a.get_text('href', '<>')
        if 'pdf' in t.lower():
            found_pdf_urls.add(a.get('href', '<>'))
    return found_pdf_urls


def __find_pdf_urls(article_url: str) -> set[str]:
    html = get_raw(article_url)
    found_pdf_urls: set[str] = set()
    for a in html.find_all('a', href=True):
        h = a.get('href', '<>')
        if h.endswith('.pdf'):
            found_pdf_urls.add(urljoin(article_url, h))
    return found_pdf_urls


def __find_pdf_urls_pmc(data: dict) -> set[str]:
    pmid: str = data[KEY_PMID]
    article_url: str = TEMP_URL_PMC_ID.format(pmid)
    return __find_pdf_urls(article_url)


def __find_new_articles(query: str, existing_article_ids: set[str]) -> list[dict]:
    found_articles: list[dict] = articles_for(query)
    logger.info(f'Found {len(found_articles)} articles for query \'{query}\'.')
    return [a for a in found_articles if a[KEY_PMID] not in existing_article_ids]


def __load_doc_record() -> list[dict]:
    if os.path.isfile(PATH_DOC_RECORD):
        logger.info(f'Loading {PATH_DOC_RECORD}.')
        with open(PATH_DOC_RECORD, 'r', encoding='utf-8') as f:
            try:
                data = json.load(f)
            except Exception as e:
                logger.error(f'Error loading {PATH_DOC_RECORD}: {e}')
                return []
        if data and isinstance(data, list) and isinstance(data[0], dict):
            return data
        logger.warning(f'Error loading {PATH_DOC_RECORD}. File content does not seem to be a list of '
                       f'dictionaries: {data}')
        return []

    logger.info(f'File {PATH_DOC_RECORD} does not exist.')
    return []


def __write_doc_record(data: list[dict]) -> None:
    with open(PATH_DOC_RECORD, 'w', encoding='utf-8') as f:
        try:
            json.dump(data, f)
            logger.info(f'Writing {PATH_DOC_RECORD}.')
        except Exception as e:
            logger.error(f'Error writing {PATH_DOC_RECORD}: {e}')


def __get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', help='Query', type=str, default='')
    parser.add_argument('-d', '--download', help='Download', action='store_true')
    parser.add_argument('-t', '--dry', help='Dry Run', action='store_true')
    return parser


if __name__ == '__main__':
    setup_logging()
    arg_parser: argparse.ArgumentParser = __get_arg_parser()
    args: argparse.Namespace = arg_parser.parse_args()
    main(query=args.query, download=args.download, dry_run=args.dry)
