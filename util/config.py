# Chat
from util.constants import KEY_TITLE, KEY_PUB_DATE, KEY_DOI

CHAT_HOST: str = '127.0.0.1'
CHAT_PORT: int = 5000

# LLM
LLM_MODEL: str = 'llama3.2:1b'
K_VALUE: int = 5
TEMPERATURE: float = 0.5
LLM_URI: str = 'http://127.0.0.1:11434/api/generate'
PROMPT_TEMPLATE: str = '''
    Answer the question based on following context:
    {context}
    ---
    Answer the question based on the above context only: {question} 
    '''
TIME_OUT: int = 1000
MAX_RETRIES: int = 5

# Embedding
EMBEDDING_MODEL: str = 'nomic-embed-text'
CHUNK_SIZE: int = 800
CHUNK_OVERLAP: int = 80
UPDATE_SIZE: int = 250
# enables local reset of chroma db
AUTH_CHROMA: str = '2d1a968ef85ee447222271b6b0624da8174092da3b74b21a5ede09d557cc30cb'

# Pubmed
RET_MAX_PM: int = 25
ADD_INFO_KEYS: list[str] = [KEY_TITLE, KEY_PUB_DATE, KEY_DOI]
USED_MAIL: str = ''

# DuckDuckGo
RET_MAX_DDG: int = 3

# HTTP
TIMEOUT_GLOBAL: int = 10
