Using VIC with Docker
========

[Docker](http://www.docker.com) is an open-source application container engine. Docker containers allow us to distribute a light-weight, virtual machine that includes the required dependencies to run VIC. VIC developers maintain a single [Ubuntu](http://www.ubuntu.com/) docker image with the minimum requirements to run the VIC **classic** and **image** drivers.

# Getting Started with Docker
Docker images may be run on most computers running Linux, Mac OS X, or Windows operating systems. To run docker images, **Docker Engine** or "Docker" must be installed on the host computer. The Docker [documentation](https://docs.docker.com/) provides specific instructions for installing Docker for each support operating system.

# Getting the VIC Docker Image
We provide a single stable light-weight Docker image that includes the minimum requirements to run VIC (e.g. GCC, netCDF, MPI). This image is distributed via [Docker Hub](https://hub.docker.com/) and is called `vic_base_image`.

```Bash
# Pull the vic_base_image
docker pull uwhydro/vic:vic_base_image
# List available docker images, vic_base_image should be listed
docker images
# Run the vic_base_image, this will just build the VIC drivers and exit
docker run -i uwhydro/vic:vic_base_image
# Open an interactive bash command line in the vic_base_image
docker run -it uwhydro/vic:vic_base_image /bin/bash
```

# Running VIC Inside a Docker Image
Once inside an interactive session inside the `vic_base_image`, we already have the VIC Git repository and all the dependencies to build and run VIC.

```Bash
# change directories to the VIC Classic Driver
cd VIC/vic/drivers/classic
# Build the Classic Driver
make
# Run VIC
./vic_classic.exe -g global_param_file.txt
```

# Testing with Docker Images
We are currently building a test platform that will be run using Docker.  Stay tuned for more information on this.

References:
- [Github Issue #50](https://github.com/UW-Hydro/VIC/issues/50)
- [Github Issue #79](https://github.com/UW-Hydro/VIC/issues/50)
- [Github Issue #350](https://github.com/UW-Hydro/VIC/issues/350)
- [wercker](http://wercker.com/)

# Other Useful Docker Tips

### Adding other software to the `vic_base_image`
Docker images can be treated just like a normal unix machine, so software can be added using package managers like [`apt-get`](https://help.ubuntu.com/community/AptGet/Howto) or [`anaconda`](https://anaconda.org).  For example, to add the `vim` text editor:

```Bash
apt-get install vim
```

### Mounting external directories to the `vic_base_image`
The `vic_base_image` doesn't include any data to run VIC, so most users will need to either download data into their image or, more commonly, mount a local volume to the Docker image. Docker makes this easy, this command mounts a local directory (e.g. `~/workdir/VIC_data/`) to the `vic_base_image`.

```Bash
# Open an interactive bash command line in the vic_base_image
# mounting a local directory
docker run -it -v ~/workdir/VIC_data:~/workdir/VIC_data uwhydro/vic:vic_base_image /bin/bash
```
