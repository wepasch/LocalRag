from __future__ import annotations
from enum import Enum

class DocumentStatus(Enum):
    __value: int
    __label: str

    def __init__(self, value: int, label: str):
        self.__value = value
        self.__label = label

    NEW = 0, 'new'
    DOWNLOADED = 1, 'downloaded'
    DOWNLOAD_FAILED = -1, 'download failed'
    TO_DELETE = 2, 'to delete'

    @property
    def value(self) -> int:
        return self.__value

    @property
    def label(self) -> str:
        return self.__label

    @classmethod
    def get_enum_label(cls) -> str:
        return 'document type'

    @classmethod
    def from_value(cls, i: int) -> DocumentStatus | None:
        return next(filter(lambda x: x.value == i, DocumentStatus.__members__.values()), None)
