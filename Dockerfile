# The default image is available on DockerHub:
#   $ docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner naokihori/snss:latest
# You can customise it by modifying this file and build a Docker image by yourself.
#   $ docker build -t naokihori/snss:latest .
#   $ docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner naokihori/snss:latest
# NOTE: Change arg of cpuset-cpus "0-1" depending on your CPU resources

FROM ubuntu:latest
RUN apt-get -y update
RUN apt-get -y install make libopenmpi-dev libfftw3-dev

ARG RUNNER=runner
RUN adduser --disabled-password ${RUNNER}
RUN echo "${RUNNER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN chown -R ${RUNNER}:${RUNNER} /home/${RUNNER}
WORKDIR /home/${RUNNER}
USER ${RUNNER}
