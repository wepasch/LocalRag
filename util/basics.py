import os
import yaml
import logging.config

from pathlib import Path
from pypdf import PdfReader
from pypdf.errors import PdfReadError
from hashlib import sha256
from typing import List, TypeVar

from util.constants import RE_YEAR

LOG_CONFIG_PATH: str = 'resources/util/logging_config.yaml'
T = TypeVar('T')


def setup_logging() -> None:
    with open(LOG_CONFIG_PATH, 'r') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)


def get_resource_dir() -> Path:
    return Path(get_root_dir().joinpath('resources'))


def get_root_dir() -> Path:
    return Path(__file__).parent.parent.parent.parent


def is_pdf(path: str) -> bool:
    try:
        PdfReader(path)
        return True
    except PdfReadError:
        return False


def delete_file(path: str) -> None:
    os.remove(path)


def add_if_available(src: dict, key: str, dst: dict) -> None:
    if key in src:
        dst[key] = src[key]


def extract_first_year(s: str) -> str | None:
    year: list[str] = RE_YEAR.findall(s)
    if year:
        return year[0]
    return None


def __split_in_lists_of(lst: List[T], max_size: int) -> List[List[T]]:
    return [lst[i:i + max_size] for i in range(0, len(lst), max_size)]


def verify(auth: str, _hash: str) -> bool:
    auth_hash = sha256(auth.encode('utf-8')).hexdigest()
    return _hash == auth_hash
