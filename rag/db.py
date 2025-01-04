import argparse
import logging
import os.path

from langchain_community.document_loaders import PyPDFDirectoryLoader as PyPDFLoader
from langchain_ollama import OllamaEmbeddings
from langchain_chroma import Chroma
from langchain.schema.document import Document
from langchain_text_splitters import RecursiveCharacterTextSplitter as RCSpliter

from util.config import CHUNK_SIZE, CHUNK_OVERLAP, EMBEDDING_MODEL, AUTH_CHROMA, UPDATE_SIZE
from util.basics import setup_logging, verify, __split_in_lists_of
from util.constants import KEY_SRC, KEY_IDS, KEY_ID, KEY_PAGE, PATH_LIBRARY_PDF

logger = logging.getLogger(__name__)

CHROMA_PATH: str = 'data/chroma'

def __embed_new_docs(verbose: bool = False) -> None:
    """
    Embeds new docs in PDF_LIB_DIR into chroma db.
    """
    documents: list[Document] = __load_db_documents(PATH_LIBRARY_PDF)
    chunks: list[Document] = __split_documents(documents)
    logger.info(f'Generated {len(chunks)} chunks from {len(documents)} pages in {len(documents)} documents')
    __add_to_chroma(chunks, verbose=verbose)


def __split_documents(documents: list[Document]) -> list[Document]:
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


def __add_to_chroma(chunks: list[Document], verbose: bool = False) -> None:
    """
    :param chunks: added to chroma db, if their ids are not in the db already
    :param verbose: log existing chunk ids
    :return:
    """
    __add_chunk_ids(chunks)
    db: Chroma = __get_db()
    existing_ids: list[str] = __get_all_ids(db)
    new_chunks: list[Document] = list(filter(lambda c: c.metadata[KEY_ID] not in existing_ids, chunks))
    if verbose:
        logger.info(f'Existing chunk ids: {sorted(existing_ids)}')
    logger.info(f'New chunk ids: {[c.metadata[KEY_ID] for c in sorted(new_chunks,
                                                                    key=lambda chunk: chunk.metadata[KEY_ID])]}')
    if new_chunks:
        logger.info(f'Adding {len(new_chunks)} new chunks to {len(existing_ids)} existing chunks...')
        __db_add_in_batches(new_chunks, db)
    else:
        logger.info(f'No new chunks to add existing {len(existing_ids)} chunks.')


def __db_add_in_batches(chunks: list[Document], db: Chroma) -> None:
    batches: list[list[Document]] = __split_in_lists_of(chunks, UPDATE_SIZE)
    idb: int
    batch: list[Document]
    for idb, batch in enumerate(batches):
        logger.info(f'Add chunks {idb * UPDATE_SIZE}-{idb * UPDATE_SIZE + len(batch)}.')
        db.add_documents(batch, ids=[b.metadata['id'] for b in batch])


def __add_chunk_ids(chunks: list[Document]) -> None:
    """
    Add id to each chunk's metadata. An id consists of 'relative path to original doc:page:part'
    """
    last_page: int = -1
    running_id: int = -1
    for chunk in chunks:
        if last_page != chunk.metadata[KEY_PAGE]:
            last_page = chunk.metadata[KEY_PAGE]
            running_id = 0
        else:
            running_id += 1
        db_id: str = f'{os.path.basename(chunk.metadata[KEY_SRC])}:{chunk.metadata[KEY_PAGE] + 1}:{running_id + 1}'
        chunk.metadata[KEY_ID] = db_id


def __load_db_documents(pdf_dir: str) -> list[Document]:
    """
    :return: documents at given path
    """
    logger.info(f'Load documents from {pdf_dir}.')
    loader = PyPDFLoader(pdf_dir)
    return loader.load()

def get_most_fitting_chunks(query: str, k_value: int) -> list[tuple[Document, float]]:
    db: Chroma = __get_db()
    return db.similarity_search_with_score(query, k=k_value)


def __get_db() -> Chroma:
    return Chroma(persist_directory=CHROMA_PATH, embedding_function=__get_embedding_function())


def __get_embedding_function():
    return OllamaEmbeddings(model=EMBEDDING_MODEL)


def __get_all_ids(db: Chroma) -> list[str]:
    existing_items = db.get(include=[])
    return existing_items[KEY_IDS]


def __get_db_stats() -> dict:
    ids: list = __get_all_ids(__get_db())
    return {
        'size': len(ids),
        'ids': ids
    }

def __delete_all(auth: str) -> bool:
    if verify(auth, AUTH_CHROMA):
        db: Chroma = __get_db()
        db.delete(ids=__get_all_ids(db))
        logger.info('Deleted chroma db contents.')
        return True
    else:
        logger.error(f'Auth error for chroma delete.')


def main(update: bool = False, verbose: bool = False, delete: str = '') -> None:
    if verbose:
        db_info: dict = __get_db_stats()
        logger.info(f'DB stats: {db_info}')
    if delete:
        if __delete_all(delete):
            exit(0)
        else:
            exit(1)
    if update:
        __embed_new_docs(verbose=verbose)


def __get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--update', help='Update', action='store_true')
    parser.add_argument('-v', '--verbose', help='Verbose', action='store_true')
    parser.add_argument('-x', '--delete', help='Delete Chroma db', type=str, default='')
    return parser


if __name__ == '__main__':
    setup_logging()
    arg_parser: argparse.ArgumentParser = __get_arg_parser()
    args: argparse.Namespace = arg_parser.parse_args()
    main(update=args.update, verbose=args.verbose, delete=args.delete)
