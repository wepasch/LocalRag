# rag
## cli rag
### Prerequisites
- running ollama docker container with language and embedding (e.g. llama3.2:1b, nomic embed text)<br>
  **Steps**
  - initialize container from terminal with running *docker desktop* and command:<br>
  _docker run -d -v ollama:/root/.ollama -p 11434:11434 --name ollama ollama/ollama_
  - run the following commands in *docker desktop* under the *Exec* tab of the **ollama** container :<br>
  >ollama pull llama3.2:1b_<br>

  >ollama pull nomic-embed-text

### Usage
- de.thb.tm.rag.rag_cli.py can be run with the *-u* flag to update chroma db with files from PDF_LIB_DIR. Files in 
that directory will be embedded by the embedding model into the chroma db. If there is no existing chroma db, it should 
be created at CHROMA_PATH in the first run.
- de.thb.tm.rag.rag_cli.py can be run with the *-q* option . *-q* is followed by the question for the rag system. 
The resulting dialog is printed to console. Example:<br>
> python de/thb/tm/rag/rag_cli.py -q "Can the mckenzie method reduce pain and if it does which kind of pain?"


## rag frontend
### Usage
Run from project directory:<br>
_streamlit run de/thb/tm/rag/rag_frontend.py --server.port 12345_
