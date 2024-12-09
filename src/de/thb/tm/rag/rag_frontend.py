import enum

import streamlit as st


KEY_CONTENT: str = 'content'
KEY_USER: str = 'user'
KEY_ROLE: str = 'role'


class Role(enum.Enum):
    __label: str

    def __init__(self, label: str) -> None:
        self.__label = label

    USER = KEY_USER
    ROLE = KEY_ROLE

    @property
    def label(self) -> str:
        return self.__label



def get_query():
    if prompt := st.chat_input('Question?'):
        with st.chat_message('user'):
            st.markdown(prompt)
        st.session_state.messages.append({'role': 'user', 'content': prompt})


def display_messages():
    for messages in st.session_state.messages:
        with st.chat_message(messages['role']):
            st.markdown(messages['content'])


def process_last_question():
    if len(st.session_state.messages) and st.session_state.messages[-1]['role'] == 'user':
        answer: str = f'You said: {st.session_state.messages[-1]['content']}'
        st.session_state.messages.append({'role': 'assistant', 'content': answer})


if len(st.session_state) == 0:
    st.session_state.messages = []
display_messages()
get_query()
process_last_question()
