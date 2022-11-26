#!/bin/bash
# 
# Run nuSTORM tests inside Docker container
# USAGE: ./run_tests_in_docker.sh
#

# Fail on any errors, warn about unset variables. Generally make bash more usable
set -o errexit
set -o errtrace
set -o nounset
set -o pipefail

# ========================================================================================

echo "Running tests in Docker"

# Check docker daemon is running, if not give a useful error!
if (! docker stats --no-stream ); then
	echo "Docker daemon is not running or there is an error. Exiting, please start Docker!"
	exit 1
fi

# ========================================================================================

# Get path to script dir
export SOURCE="${BASH_SOURCE[0]}"
while [ -h "${SOURCE}" ]; do # resolve ${SOURCE} until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "${SOURCE}" )" && pwd )"
  SOURCE="$(readlink "${SOURCE}")"
  [[ ${SOURCE} != /* ]] && SOURCE="${DIR}/${SOURCE}" # if ${SOURCE} was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
declare DIR
DIR="$( cd -P "$( dirname "${SOURCE}" )" && pwd )"

cd "$DIR" || exit 0

# ========================================================================================

DOCKERFILE="${DIR}/Dockerfile"
# Tag the image with the current git description
DOCKER_TAG=$(git rev-parse HEAD)

echo "Building Docker Image"

# Build the nuStorm container image
docker build \
    -t "nustorm:${DOCKER_TAG}" \
    -t "nustorm:latest" \
    -f "${DOCKERFILE}" "${DIR}"

# Run the container and the tests inside
echo "Running tests in docker"
docker run --rm \
    "nustorm:${DOCKER_TAG}" \
    /nustorm/02-Tests/intregration-test.bash
