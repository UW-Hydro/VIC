FROM ubuntu:trusty

RUN apt-get update -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        build-essential \
        ssh \
        netcdf-bin \
        libnetcdf-dev \
        libopenmpi-dev \
        openmpi-bin \
        git \
        make

CMD git clone https://github.com/jhamman/VIC.git && \
    cd VIC/vic/drivers/classic && \
    make && \
    cd ../image && \
    make
