import argparse
import httpx
import json
import logging.config

from httpx import Response
from http.client import HTTPException
from langchain_core.prompts import ChatPromptTemplate
from langchain_community.document_loaders import PyPDFDirectoryLoader as PyPDFLoader
# the following import is deprecated, but change to line below breaks
from langchain_community.embeddings.ollama import OllamaEmbeddings
#from langchain_ollama.embeddings import OllamaEmbeddings
from langchain_chroma import Chroma
from langchain.schema.document import Document
from langchain_text_splitters import RecursiveCharacterTextSplitter as RCSpliter

from de.thb.util import setup_logging

CHROMA_PATH: str = 'data/chroma'
BASE_URL: str = "http://localhost:11434/api/generate"
LLM_MODEL: str = 'llama3.2:1b'
EMBEDDING_MODEL: str = 'nomic-embed-text'
PDF_LIB_DIR: str = 'data/library/pdf'
PROMPT_TEMPLATE: str = '''
    Answer the question based on following context:
    {context}
    ---
    Answer the question based on the above context only: {question} 
    '''
K_VALUE: int = 5
CHUNK_SIZE: int = 800
CHUNK_OVERLAP: int = 80
TIME_OUT: int = 1000
MAX_RETRIES: int = 5

logger = logging.getLogger(__name__)


class QueryRecord:
    __query: str
    __scored_docs: list[tuple[Document, float,]] = []
    __answer: str = None

    def __init__(self, query: str) -> None:
        self.__query = query

    @property
    def query(self) -> str:
        return self.__query

    @property
    def scored_docs(self) -> list[tuple[Document, float]]:
        """
        :return: docs with similarity score
        """
        return self.__scored_docs

    @property
    def docs(self) -> list[Document]:
        return [sd[0] for sd in self.scored_docs]

    @property
    def answer(self) -> str:
        return self.__answer

    def add_scored_docs(self, scored_docs: list[tuple[Document, float,]]) -> None:
        self.__scored_docs.extend(scored_docs)

    def get_context(self) -> str:
        """
        :return: constructed context from contents of scored docs
        """
        if len(self.__scored_docs) == 0:
            logger.error(f'No scored documents found for query: {self.__query}')
        return '\n\n--\n\n'.join(chunk.page_content for chunk in self.docs)

    def add_answer(self, answer: str) -> None:
        self.__answer = answer

    def get_full_prompt(self) -> str:
        """
        :return: prompt consisting of question and context.
        """
        prompt_template: ChatPromptTemplate = ChatPromptTemplate.from_template(PROMPT_TEMPLATE)
        return prompt_template.format(context=self.get_context(), question=self.query)

    def get_sources_list(self) -> str:
        """
        :return: lines from scored documents with page and part numbers
        """
        sources_split: list[list[str]] = [c[0].metadata['id'].split(':') for c in self.scored_docs]
        return '\n'.join([f"{f[0]} pg.{f[1]} - pt.{f[2]}" for f in sources_split])


def process_query(query_record: QueryRecord, temperature: float | None = None) -> None:
    """
    Collect fitting documents, add context from documents, send prompt with question and context to llm and add answer
    to QueryRecord
    :param query_record: with question
    :param temperature: temperature parameter for llm in [0,1]
    """
    query: str = query_record.query
    logger.info('Query: %s', query)

    db: Chroma = get_db()
    most_fitting_chunks: list[tuple[Document, float,]] = db.similarity_search_with_score(query, k=K_VALUE)
    query_record.add_scored_docs(most_fitting_chunks)
    logger.debug(f'Query \'{query}\' lead to following chunks:\n\t'
                 f'{"\n\t".join([f"{c[1]} {c[0].metadata['id']}" for c in most_fitting_chunks])} ... constructing\n\t\t'
                 f'\'{query_record.get_full_prompt()}\'')
    data: dict = {
        'model': LLM_MODEL,
        'prompt': query_record.get_full_prompt()
    }
    if temperature:
        if temperature < 0:
            logger.warning('Negative temperature detected. Set it to 0.')
            temperature = 0
        elif temperature > 1:
            logger.warning('Temperature greater 1 detected. Set it to 1.')
            temperature = 1
        logger.info(f'Using temperature of {temperature}')
        data['options'] = {"temperature": temperature}

    full_response: str = ollama_response_to_str(
        httpx.post(BASE_URL, data=json.dumps(data), headers={'Content-Type': 'application/json'}, timeout=TIME_OUT))
    query_record.add_answer(full_response)


def print_dialog(query_record: QueryRecord, verbose: bool = False) -> None:
    """
    :param query_record: is printed to console
    :param verbose: print full context too
    """
    if verbose:
        logger.info(f'Full prompt will be:\n{query_record.get_full_prompt()}\n\n\n')
    print(f'''
Question:
{query_record.query}
    
Answer:
{query_record.answer}

Sources:
{query_record.get_sources_list()}''')


def ollama_response_to_str(response: Response) -> str:
    """
    Extracts answer from ollama request as a string.
    :param response: from api request to ollama
    :return:
    """
    if response.status_code != 200:
        raise HTTPException(f'Could not retrieve data from {BASE_URL}. [Code: {response.status_code}]')
    response_lines: list[str] = [line for line in response.text.strip().split('\n') if line]
    response_dicts: list[dict] = [json.loads(line) for line in response_lines]
    return ''.join(response_dict.get('response', '<?>') for response_dict in response_dicts)


def add_to_chroma(chunks: list[Document]) -> None:
    """
    TODO research: chunks get identified as new although they were added already
    :param chunks: added to chroma db, if their ids are not in the db already
    :return:
    """
    add_chunk_ids(chunks)
    db: Chroma = get_db()
    existing_items = db.get(include=[])
    existing_ids = set(existing_items['ids'])
    new_chunks: list[Document] = list(filter(lambda c: c.metadata['id'] not in existing_ids, chunks))
    logger.info(f'Existing chunk ids: {sorted(existing_ids)}')
    logger.info(f'New chunk ids: {[c.metadata['id'] for c in sorted(new_chunks,
                                                                    key=lambda chunk: chunk.metadata['id'])]}')
    if new_chunks:
        logger.info(f'Adding {len(new_chunks)} new chunks to {len(existing_ids)} existing chunks...')
        db.add_documents(new_chunks, ids=[c.metadata['id'] for c in chunks])
    else:
        logger.info(f'No new chunks to add existing chunks.')


def add_chunk_ids(chunks: list[Document]) -> None:
    """
    Add id to each chunk's metadata. An id consists of 'relative path to original doc:page:part'
    """
    last_page: int = -1
    running_id: int = -1
    for chunk in chunks:
        if last_page != chunk.metadata['page']:
            last_page = chunk.metadata['page']
            running_id = 0
        else:
            running_id += 1
        db_id: str = f'{chunk.metadata['source']}:{chunk.metadata['page'] + 1}:{running_id + 1}'
        chunk.metadata['id'] = db_id


def load_documents(pdf_dir: str) -> list[Document]:
    """
    :return: documents at given path
    """
    logger.info(f'Load documents from {pdf_dir}.')
    loader = PyPDFLoader(pdf_dir)
    return loader.load()


def split_documents(documents: list[Document]) -> list[Document]:
    """
    :return: chunks of documents
    """
    splitter = RCSpliter(
        chunk_size=CHUNK_SIZE,
        chunk_overlap=CHUNK_OVERLAP,
        length_function=len,
        is_separator_regex=False
    )
    return splitter.split_documents(documents)


def embed_new_docs() -> None:
    """
    Embeds new docs in PDF_LIB_DIR into chroma db.
    """
    sources_documents: set = set()
    documents: list[Document] = load_documents(PDF_LIB_DIR)
    [sources_documents.add(document.metadata['source']) for document in documents]
    chunks: list[Document] = split_documents(documents)
    logger.info(f'Generated {len(chunks)} chunks from {len(documents)} pages in {len(sources_documents)} documents')
    add_to_chroma(chunks)


def get_db() -> Chroma:
    return Chroma(persist_directory=CHROMA_PATH, embedding_function=get_embedding_function_depr())


def get_embedding_function_depr():
    """
    TODO Replace import to not deprecated OllamaEmbeddings
    """
    return OllamaEmbeddings(model=EMBEDDING_MODEL)


def get_embedding_function():
    return OllamaEmbeddings(model=EMBEDDING_MODEL)


def get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", help="Query")
    parser.add_argument("-u", "--update", help="Update doc database", action='store_true')
    return parser


def main(query: str = None, update_db: bool = False) -> None:
    if update_db:
        embed_new_docs()
    if query:
        query_object: QueryRecord = QueryRecord(query)
        process_query(query_object)
        print_dialog(query_object, verbose=True)


if __name__ == '__main__':
    arg_parser: argparse.ArgumentParser = get_arg_parser()
    args: argparse.Namespace = arg_parser.parse_args()
    setup_logging()
    main(args.query if args.query else '', args.update)

