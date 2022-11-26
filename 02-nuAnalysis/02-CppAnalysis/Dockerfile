FROM rootproject/root:latest

LABEL maintainer="Stephen Kyberd <stephen.kyberd@oxbotica.com>"

# =============================================================

# Set up user/groups
ARG USERNAME=nustorm
ARG UID=1000
ARG GID=1000

RUN groupadd -g $GID -o $USERNAME \
    && useradd -m -u $UID -g $GID -o -s /bin/bash $USERNAME \
    && adduser $USERNAME sudo \
    && printf "$USERNAME ALL= NOPASSWD: ALL\\n" >> /etc/sudoers

# =============================================================

# Install minimally required system packages
RUN apt-get -qq update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    ninja-build \
    && rm -rf /var/lib/apt/lists/*


# Run in dedicated directory
WORKDIR /nustorm

# Copy the source code in to the container
COPY --chown=$USERNAME . /nustorm

# Create a build folder
WORKDIR /nustorm/build

# Run CMake
RUN cmake ../src --preset default

# Run Ninja to build the code
RUN ninja


ENTRYPOINT ["/bin/bash"]
