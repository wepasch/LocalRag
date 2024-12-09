import json
import os
import pymupdf
import requests
from Bio import Entrez
from bs4 import BeautifulSoup


PHYSIO_QUERY: str = 'physiotherapy internal diseases'
# add if needed
MAIL: str = ''
RECORD_PATH: str = '../../../../../record.json'
PDF_DIR: str = '../../../../../pdf'
TXT_DIR: str = '../../../../../txt'
SUFFIX_PDF: str = '.pdf'
SUFFIX_TXT: str = '.txt'
KEY_SUCCESS: str = 'successful'
KEY_FAIL: str = 'failed'


def get_pmids(query: str, email: str, max_results: int = 5) -> list[str]:
    """
    Search for query in title or abstract and return pymed ids.
    """
    Entrez.email = email
    with Entrez.esearch(
            db='pubmed',
            term=f'{query} AND full text[filter]',
            retmax=max_results,
            field='title/abstract'
    ) as handle:
        results = Entrez.read(handle)
    return results['IdList']


def get_doi(fetch_result: str) -> str | None:
    """
    Construct standard doi url from Entrez.efetch result. None if doi can not be determined.
    """
    lines: list[str] = fetch_result.split('\n')
    doi_line: str = next(filter(lambda line: 'LID - 10.' in line, lines), None)
    try:
        return doi_line.split(' ')[2]
    except IndexError:
        print(f'No doi for {fetch_result}')
        return None
    except AttributeError:
        print(f'No doi for {fetch_result} [attribute error]')


def get_doi_url(pmid: str, email: str) -> str | None:
    """
    Construct standard doi url from pymed id.
    """
    Entrez.email = email
    with Entrez.efetch(
            db='pubmed',
            id=pmid,
            rettype='medline',
            retmode='text'
    ) as handle:
        result = handle.read()
        doi: str | None = get_doi(result)
    if doi:
        return f'https://doi.org/{doi}'
    else:
        return None


def write_pdf(doi_url: str, file_name: str) -> bool:
    """
    Fetch from doi url, search for links in resulting html. The first link that contains 'pdf' is tried
    to be downloaded and written to file_name.
    """
    with requests.Session() as session:
        try:
            response = session.get(doi_url)
            response.raise_for_status()
            soup = BeautifulSoup(response.content, 'html.parser')
            link_candidates = soup.find_all('a', href=True)
            pdf_link = None
            for link in link_candidates:
                if 'pdf' in link['href'].lower():
                    pdf_link_raw = link['href']
                    pdf_link = pdf_link_raw.replace("//content", "/content")
                    break

            if pdf_link and not pdf_link.startswith('http'):
                url_comps = response.url.rsplit('/')
                prefix = '/'.join(url_comps[:3])
                pdf_link = prefix + '/' + pdf_link

            if pdf_link:
                print(f'Downloading PDF from: {pdf_link}')
                pdf_response = session.get(pdf_link)
                pdf_response.raise_for_status()

                with open(f'{PDF_DIR}/{file_name}{SUFFIX_PDF}', 'wb') as pdf_file:
                    pdf_file.write(pdf_response.content)
                print(f'PDF downloaded successfully as \'{file_name}\'.')
                return True
            else:
                print('PDF link not found on the page.')

        except requests.exceptions.HTTPError as e:
            print(f'HTTP error occurred: {e}')
        except requests.exceptions.RequestException as e:
            print(f'Error occurred: {e}')
    return False


def write_pdfs_one_by_one(pmids: list[str], mail: str, record: dict[str, list[str]]) -> None:
    """
    Tries to obtain articles of pymed ids as pdf write them as <pymed id>.pdf. Tracks success/failure in record.
    """
    nof_pmids: int = len(pmids)
    idp: int
    pmid: str
    for idp, pmid in enumerate(pmids):
        print(f'Process pmid \'{pmid}\' {idp + 1}/{nof_pmids}')
        doi_url: str = get_doi_url(pmid, mail)
        if not doi_url:
            print(f'Skip pmid {pmid}')
            continue
        was_written: bool = write_pdf(doi_url, f'{pmid}')
        key: str
        if was_written:
            key = KEY_SUCCESS
        else:
            key = KEY_FAIL
        record.get(key).append(pmid)
        write_record(record)


def get_record() -> dict:
    """
    :return: content of record file
    """
    if not os.path.isfile(RECORD_PATH):
        print(f'No record file at {RECORD_PATH}.')
        return {KEY_SUCCESS: [], KEY_FAIL: []}
    with open(RECORD_PATH, 'r') as record_file:
        record: dict = json.load(record_file)
        if KEY_SUCCESS in record and KEY_FAIL in record:
            return record
        else:
            return {KEY_SUCCESS: [], KEY_FAIL: []}


def write_record(data: dict) -> None:
    """
    :param data: is written into record
    """
    with open(RECORD_PATH, 'w') as record_file:
        json.dump(data, record_file)


def convert_pdfs() -> None:
    """
    Converts all <wd>/pdf/<name>.pdf files into <wd>/txt/<name> files by concatenating all page strings of a pdf.
    """
    full_pdfs: list[str] = next(os.walk(PDF_DIR), (None, None, []))[2]
    for thing in full_pdfs:
        doc: pymupdf.Document = pymupdf.open(f'{PDF_DIR}/{thing}')
        txt_str: str = ''
        page: pymupdf.Page
        for page in doc:
            txt_str += page.get_text("text")
        txt_str = txt_str.replace('\n', '')
        with open(f'{TXT_DIR}/{thing.replace(SUFFIX_PDF, SUFFIX_TXT)}', 'w', encoding='utf-8') as txt_file:
            txt_file.write(txt_str)


def main():
    # load record of already processed pmids
    record: dict[str, list[str]] = get_record()
    # find query related pmids
    existing_pmids: set[str] = set(record.get(KEY_SUCCESS, set())) | set(record.get(KEY_FAIL, set()))
    found_pmids: list[str] = get_pmids(PHYSIO_QUERY, MAIL, max_results=9999)
    # identify new pmids
    new_pmids: list[str] = list(set(found_pmids) - existing_pmids)
    # write pdfs from pymed ids
    write_pdfs_one_by_one(new_pmids, MAIL, record)
    # convert pdfs to text
    convert_pdfs()




if __name__ == '__main__':
    main()
