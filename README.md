# Local Rag
### Prerequisites
- running ollama docker container with llm and embedding (e.g. llama3.2:1b, nomic embed text)<br>
  **Steps**
  - initialize container from terminal with running *docker desktop* and command:<br>
  _docker run -d -v ollama:/root/.ollama -p 11434:11434 --name ollama ollama/ollama_
  - run the following commands in *docker desktop* under the *Exec* tab of the **ollama** container :<br>
  >ollama pull llama3.2:1b<br>

  >ollama pull nomic-embed-text

If other embedding or llm is used, edit entries EMBEDDING_MODEL or LLM_MODEL in util/config.py

- activated environment with python interpreter (requirements.txt coming up)
- (set CHAT_HOST in util/config.py to host ip, default is 127.0.0.1)
- (add pdfs custom into ./data/library/pdf)

If executed modules can not be found, add project path to PYTHONPATH.
___
## Update documents
> python rag/update.py -q &lt;query&gt;

Enter search term after -q. Pubmed will be searched for free, fulltext articles. 
The articles' data is recorded in *data/record/doc_record.json*.

> python rag/update.py -q &lt;query&gt;

Attempts to download recorded articles with new property status==NEW. Tries to get download url 
from pubmed and then duckduckgo. If successfull downloads are saved into *./data/library/pdf*.
___
## Update chroma db
> python rag/db.py -u

Chunks and adds documents in *data/library/pdf* to chroma db. Run only with flag -v to list db size and ids. 
___
## Run app
>python rag/chat.py

App is running, when following line is printed *...INFO: Start with 127.0.0.1:5000...* 
and can be accessed under the given url and port. 