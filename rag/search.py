import argparse
import httpx
import json
import logging.config

from httpx import Response
from http.client import HTTPException
from langchain.schema.document import Document

from util.config import LLM_MODEL, K_VALUE, LLM_URI, TIME_OUT
from rag.query_record import QueryRecord
from rag.db import get_most_fitting_chunks
from util.basics import setup_logging

logger = logging.getLogger(__name__)

def process_query(query_record: QueryRecord, temperature: float | None = None) -> None:
    """
    Collect fitting documents, add context from documents, send prompt with question and context to llm and add answer
    to QueryRecord
    :param query_record: with question
    :param temperature: temperature parameter for llm in [0,1]
    """
    query: str = query_record.query
    logger.info('Query: %s', query)
    most_fitting_chunks: list[tuple[Document, float,]] = get_most_fitting_chunks(query, K_VALUE)
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
        httpx.post(LLM_URI, data=json.dumps(data), headers={'Content-Type': 'application/json'}, timeout=TIME_OUT))
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
        raise HTTPException(f'Could not retrieve data from {LLM_URI}. [Code: {response.status_code}]')
    response_lines: list[str] = [line for line in response.text.strip().split('\n') if line]
    response_dicts: list[dict] = [json.loads(line) for line in response_lines]
    return ''.join(response_dict.get('response', '<?>') for response_dict in response_dicts)


def get_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", help="Query", type=str, required=True)
    return parser


def main(query: str = None) -> None:
    query_object: QueryRecord = QueryRecord(query)
    process_query(query_object)
    print_dialog(query_object, verbose=True)


if __name__ == '__main__':
    arg_parser: argparse.ArgumentParser = get_arg_parser()
    args: argparse.Namespace = arg_parser.parse_args()
    setup_logging()
    main(args.query)

