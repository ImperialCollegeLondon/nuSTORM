FROM rootproject/root:latest

# Set up user/groups
ARG USERNAME=nustorm
ARG UID=1000
ARG GID=1000

RUN groupadd -g $GID -o $USERNAME \
    && useradd -m -u $UID -g $GID -o -s /bin/bash $USERNAME \
    && adduser $USERNAME sudo \
    && printf "$USERNAME ALL= NOPASSWD: ALL\\n" >> /etc/sudoers

# Run in dedicated directory
WORKDIR /nustorm

RUN source ./startup.bash

