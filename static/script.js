const socket = io();

const currentUsername = document.getElementById('current-username');
const usernameInput = document.getElementById('username-input');
const updateUsernameButton = document.getElementById('button-update-username');
const chatMessages = document.getElementById('chat-messages');
const messageInput = document.getElementById('message-input');
const sendButton = document.getElementById('send-chat-button');
const askButton = document.getElementById('send-ask-button');

function addMessage(content, isSystem = false) {
    const messageElement = document.createElement('div');
    messageElement.classList.add('message');
    if (isSystem) {
        messageElement.classList.add('system-message');
    }
    messageElement.innerHTML = content;
    chatMessages.appendChild(messageElement);
    chatMessages.scrollTop = chatMessages.scrollHeight;
}

function sendMessage() {
    const message = messageInput.value.trim();
    if (message) {
        socket.emit('send_message', {message});
        messageInput.value = '';
    }
}

function submitNewName() {
    const newUsername = usernameInput.value.trim();
    if (newUsername) {
        socket.emit('update_username', {username: newUsername});
        usernameInput.value = '';
    }
}

function askBot() {
    const message = messageInput.value.trim();
    if (message) {
        socket.emit('ask_message', {message});
        messageInput.value = '';
    }
}

socket.on('set_username', (data) => {
    currentUsername.textContent = `Your username: ${data.username}`;
    addMessage(`<em>Welcome! Your username is <strong>${data.username}</strong>.</em>`, true);
});

socket.on('user_joined', (data) => {
    addMessage(`<em>${data.username} joined the chat.</em>`, true);
});

socket.on('user_left', (username) => {
    addMessage(`<em>${username} left the chat.</em>`, true);
});

socket.on('new_message', (data) => {
    const {username, avatar_path, message} = data;
    addMessage(`
        <div class="chat-message">
            <img src="${avatar_path}" alt="${username}" class="avatar">
            <strong>${username}</strong>: ${message}
        </div>
    `);
});

sendButton.addEventListener('click', () => {
    sendMessage()
});

messageInput.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        sendMessage()
    }
});

updateUsernameButton.addEventListener('click', () => {
    submitNewName()
});

usernameInput.addEventListener('keypress', (event) => {
    if (event.key === 'Enter') {
        submitNewName()
    }
});

askButton.addEventListener('click', () => {
    askBot()
});

socket.on('username_updated', (data) => {
    const {old_username, new_username} = data;
    if (old_username !== new_username) {
        if (currentUsername.textContent.includes(old_username)) {
            currentUsername.textContent = `Your username: ${new_username}`;
        } else {
            addMessage(`<em>${old_username} changed their name to ${new_username}.</em>`, true);
        }
    }
});

socket.on('bot_answer', (data) => {
    const {answer, links} = data
    if (Object.keys(data).length === 0) {
        addMessage(`
                <div class="bot-answer">
                    <em>...accepting no questions while searching...</em>
                </div>
        `);
    } else {
        const botMessages = document.querySelectorAll('.bot-answer');
        if (botMessages.length > 0) {
            const lastBotMessage = botMessages[botMessages.length - 1];
            if (links.length > 0) {
                lastBotMessage.innerHTML = `
            <p>${answer}</p>
            <strong>References:</strong>
            <ul class="references-list">
                ${links.map(link => `
                    <li><a href="${link.link}" target="_blank">${link.name}</a></li>
                `).join('')}
            </ul>
        `;
            } else {
                lastBotMessage.innerHTML = `
            <p>${answer}</p>
            `;
            }
        }
    }
});

socket.on('bot_working', () => {
    askButton.disabled = true
    askButton.style.opacity = 0.5
});

socket.on('bot_done', () => {
   askButton.disabled = false
    askButton.style.opacity = 1
});