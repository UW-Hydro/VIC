FROM ubuntu:trusty

# Get the basic VIC dependencies
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

# Dedicated work directory for output
ENV WORKDIR $HOME/workdir
RUN mkdir -p $WORKDIR

# Put the UW-Hydro version of VIC in the container
RUN git clone https://github.com/UW-Hydro/VIC.git

# Command to run when this image is "run", just build the classic and image drivers
CMD git checkout develop && \
    git pull origin develop && \
    cd VIC/vic/drivers/classic && \
    make && \
    ./vic_classic.exe -o && \
    cd ../image && \
    make && \
    ./vic_image.exe -o
