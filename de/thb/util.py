import logging.config
from pathlib import Path

import yaml

LOG_CONFIG_PATH: str = 'de/thb/logging_config.yaml'


def setup_logging() -> None:
    with open(get_resource(LOG_CONFIG_PATH), 'r') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)


def get_resource_dir() -> Path:
    return Path(get_root_dir().joinpath('resources'))


def get_root_dir() -> Path:
    return Path(__file__).parent.parent.parent


def get_resource(rel_path: str) -> Path:
    return get_resource_dir().joinpath(rel_path)
