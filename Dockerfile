# Dockerfile with necessary packages to run popcorn as executable

FROM ubuntu:14.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update &&
    apt-get -y install git python-dev python-pip python-setuptools

RUN git clone https://github.com/brielin/popcorn /opt/popcorn &&
    cd /opt/popcorn &&
    python setup.py install

ENTRYPOINT ["popcorn"]
