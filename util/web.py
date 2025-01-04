import logging
import os
import requests

from bs4 import BeautifulSoup
from duckduckgo_search import DDGS
from duckduckgo_search.exceptions import DuckDuckGoSearchException

from util.config import TIMEOUT_GLOBAL, RET_MAX_DDG
from util.constants import KEY_HREF

logger = logging.getLogger(__name__)


def get_raw(url: str):
    try:
        with requests.get(url, timeout=TIMEOUT_GLOBAL, allow_redirects=True, headers={
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
                          'Chrome/91.0.4472.124 Safari/537.36'
        }) as r:
            logger.info(f'{r.status_code} from {url}')
            return BeautifulSoup(r.text, 'html.parser')
    except IOError as e:
        logger.warning(f'Can not get html. Return empty html. Url: {url}. Error: {e}')
        return BeautifulSoup('<html></html>', 'html.parser')


def download_pdf(url: str, path: str) -> bool:
    try:
        response = requests.get(url, stream=True, timeout=TIMEOUT_GLOBAL, allow_redirects=True, headers={
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) '
                      'Chrome/91.0.4472.124 Safari/537.36'
    })
        if response.status_code == 200:
            with open(path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=1024):
                    file.write(chunk)
            return os.path.isfile(path)
        logger.warning(f'Failed to download PDF from {url} [{response.status_code}]')
    except Exception as e:
        logger.warning(f'Failed to download PDF from {url} [{e}]')
    return False


def get_ddg_filetype_urls(query: str, filetype: str) -> set[str]:
    try:
        results: list[dict[str, str]] = DDGS().text(f'{query} filetype:{filetype}',
                              safesearch='off', timelimit='y', max_results=RET_MAX_DDG)
        found_pdf_urls: set[str] = {r[KEY_HREF] for r in results if KEY_HREF in r}
        return found_pdf_urls
    except DuckDuckGoSearchException as e:
        logger.warning(f'Can not get urls for query \'{query}\',file type \'{filetype}\': {e}')