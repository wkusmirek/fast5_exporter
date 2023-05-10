FROM        ubuntu:20.04
MAINTAINER  Kuśmirek Wiktor <kusmirekwiktor@gmail.com>

ARG TARGETARCH
ARG BIN_DIR=.

RUN apt-get update
RUN apt-get install -y python3 pip

RUN pip3 install ont-fast5-api

COPY ${BIN_DIR}/fast5_exporter /fast5_exporter
COPY ${BIN_DIR}/python /python

EXPOSE 9307
ENTRYPOINT [ "/fast5_exporter" ]
