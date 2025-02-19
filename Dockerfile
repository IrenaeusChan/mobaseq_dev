FROM debian:bullseye-slim
MAINTAINER "Irenaeus Chan <chani@wustl.edu>"

# Volumes
VOLUME /build

ARG TAG=v0.0.0
ARG PYTHON_VERSION=3.11.1

# Install build dependencies
RUN apt-get update -qq \
    && apt-get -y install apt-transport-https ca-certificates \
    && apt-get update -qq \
    && apt-get -y install \
    build-essential \
    vim \
    git \
    less \
    procps \
    libnss-sss \
    libcurl4-openssl-dev \
    curl \
    wget \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libffi-dev \
    --no-install-recommends

# Install R and dependencies
RUN apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && R --version \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get clean all

# Download and install Python
RUN wget https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz \
    && tar -xzf Python-$PYTHON_VERSION.tgz \
    && cd Python-$PYTHON_VERSION \
    && ./configure --prefix=/opt/python-$PYTHON_VERSION \
    && make \
    && make install

# Update PATH
ENV PATH="/opt/python-$PYTHON_VERSION:$PATH"

# Install pip and the package
RUN /opt/python-$PYTHON_VERSION/bin/pip3 install --upgrade pip
RUN /opt/python-$PYTHON_VERSION/bin/pip3 install --upgrade git+https://github.com/IrenaeusChan/python_dev_template.git
RUN /opt/python-$PYTHON_VERSION/bin/pip3 install --upgrade git+https://github.com/IrenaeusChan/mobaseq_dev.git
