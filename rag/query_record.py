import logging

from langchain_core.prompts import ChatPromptTemplate
from langchain.schema.document import Document

from util.config import PROMPT_TEMPLATE
from rag.db import KEY_ID

logger = logging.getLogger(__name__)


class QueryRecord:
    __query: str
    __scored_docs: list[tuple[Document, float,]]
    __answer: str = None

    def __init__(self, query: str) -> None:
        self.__query = query
        self.__scored_docs = []

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

    def get_sources_list(self) -> list[str]:
        """
        :return: lines from scored documents with page and part numbers
        """
        sources_split: list[list[str]] = [c[0].metadata[KEY_ID].split(':') for c in self.scored_docs]
        return [f"{f[0]} pg.{f[1]} - pt.{f[2]}" for f in sources_split]
