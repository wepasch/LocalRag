import eventlet

from util.config import CHAT_HOST, CHAT_PORT

eventlet.monkey_patch()

import logging
import re
import os

from random import choice, randint
from langchain_core.documents import Document
from flask import Flask, render_template, request, send_from_directory
from flask_socketio import SocketIO, emit

from rag.search import QueryRecord, process_query
from util.basics import setup_logging
from util.constants import KEY_ID, KEY_SRC, PATH_LIBRARY

CHAT_URI: str = f'http://{CHAT_HOST}:{CHAT_PORT}'

PATH_PROJECT: str = f'{os.getcwd()}{os.sep}'
PATH_STATIC: str = f'{PATH_PROJECT}static{os.sep}'
PATH_IMAGES: str = f'{PATH_STATIC}images{os.sep}'
PATH_USR_IMS: str = f'{PATH_IMAGES}users{os.sep}'
REGEX_NO_HTML = re.compile('<.*?>|&([a-z0-9]+|#[0-9]{1,6}|#x[0-9a-f]{1,6});')

logger = logging.getLogger(__name__)
app = Flask(__name__,
            static_folder=f'{os.getcwd()}{os.sep}static',
            template_folder=f'{os.getcwd()}{os.sep}templates')
socketio = SocketIO(app, async_mode='eventlet')


users: dict[str, dict] = {}

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/library/<path:filename>')
def serve_pdf(filename: str):
    return send_from_directory(PATH_LIBRARY, filename)

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(PATH_IMAGES, 'favicon.ico')

@socketio.on('connect')
def handle_connect() -> None:
    rand_user_image_path: str = PATH_USR_IMS + choice(os.listdir(PATH_USR_IMS))
    rand_user_image_path = rand_user_image_path.replace(PATH_PROJECT, '')
    new_user: dict[str, str] = {
        'username': f'User_{randint(1001,9999)}',
        'avatar_path': rand_user_image_path
    }
    users[request.sid] = new_user
    emit('user_joined', new_user, broadcast=True)
    emit('set_username', {'username': new_user['username']})


@socketio.on('disconnect')
def handle_disconnect() -> None:
    leaving_user: dict[str, str] = users.pop(request.sid, None)
    if leaving_user is None:
        return
    emit('user_left', leaving_user['username'], broadcast=True)


@socketio.on('send_message')
def handle_message(message: dict, is_ask: bool = False) -> None:
    sending_user: dict[str, str] = users.get(request.sid, None)
    if sending_user is None:
        return

    return_message: str = message.get('message', 'ERROR')
    if is_ask:
        return_message = f'<u>{return_message}</u>'
    message_data: dict[str, str] = {
        'username': sending_user['username'],
        'avatar_path': sending_user['avatar_path'],
        'message': return_message
    }
    emit('new_message', message_data, broadcast=True)

@socketio.on('update_username')
def handle_update_username(message_data) -> None:
    user: dict[str, str] = users.get(request.sid, None)
    if user is None:
        return
    old_username: str = user['username']
    new_username: str = __remove_html_tags(message_data['username'])

    user['username'] = new_username
    emit('username_updated', {'old_username': old_username,
                              'new_username': new_username}, broadcast=True, include_self=False)

    emit('username_updated', {'old_username': old_username,
                              'new_username': new_username}, room=request.sid)

@socketio.on('ask_message')
def handle_ask_message(message):
    emit('bot_working', broadcast=True)
    handle_message(message, is_ask=True)
    emit('bot_answer', {}, broadcast=True)

    user_message = message.get('message')
    bot_response: dict = __get_bot_answer(user_message)

    emit('bot_answer', bot_response, broadcast=True)
    emit('bot_done', broadcast=True)


def __get_bot_answer(question: str) -> dict:
    try:
        query_record: QueryRecord = QueryRecord(question)
        process_query(query_record)
        answer: str = query_record.answer
        links_data: list[dict[str, str]] = __extract_links(query_record.docs)
        return {
            'answer': answer,
            'links': links_data
        }
    except Exception as e:
        logger.error(f'Error processing question: {e}', exc_info=True)
        return {
            'answer': 'I\'m sorry, I couldn\'t process that question – mby I need help...',
            'links': []
        }

def __extract_links(docs: list[Document]) -> list[dict[str, str]]:
    links_data: list[dict] = []
    existing_names: set[str] = set()
    doc: Document
    for d in docs:
        id_arr: list[str] = d.metadata[KEY_ID].split(':')
        file_name: str = id_arr[0]
        page: str = id_arr[1]
        name: str = file_name
        if len(name) > 50:
            name = f'{name[:40]}...'
        name += f' (pg {page})'
        link: str = f"{CHAT_URI}/library/pdf/{file_name}#page={page}"
        if name not in existing_names:
            links_data.append({
                'name': name,
                'link': link
            })
            existing_names.add(name)
    return links_data

def __remove_html_tags(text: str) -> str:
    text = re.sub(REGEX_NO_HTML, '', text)
    return re.sub(re.compile('[<>]'), '', text)


def main(mode: str = '?'):
    if mode == 'dev':
        logger.info(f'Start dev mode with 127.0.0.1:{CHAT_PORT}.')
        socketio.run(app, port=CHAT_PORT)
    else:
        logger.info(f'Start with {CHAT_HOST}:{CHAT_PORT}.')
        socketio.run(app, host=CHAT_HOST, port=CHAT_PORT)


if __name__ == '__main__':
    setup_logging()
    main()
