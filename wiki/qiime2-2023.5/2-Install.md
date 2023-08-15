## Dependencies

Tourmaline requires the following software:

* Conda
* QIIME 2 version 2023.5
* QIIME 2 plugins: deicode, empress
* Snakemake version 7.30.1
* Python packages: biopython, tabulate
* R packages: msa, odseq
* Multiple sequence alignment tools: clustalo, muscle v5
* Other command-line tools: pandoc

## Installation options

### Option 1: Native installation

The native installation builds on the Conda installation of QIIME 2. 

#### Conda

First, if you don't have Conda installed on your machine, install [Miniconda](https://conda.io/miniconda.html) for your operating system (Python 3.8+ version).

#### Snakemake

Second, install snakemake into its own environment. Ensure it is version < 7.30.1, as v7.30.2 only works with Python 3.9+.

```
conda install -c conda-forge -c bioconda snakemake=7.30.1
```

#### QIIME 2

Third, install QIIME 2 in a Conda environment, if you haven't already. See the instructions at [qiime2.org](https://docs.qiime2.org/2023.5/install/native/). For example, on macOS these commands will install QIIME 2 inside a Conda environment called `qiime2-2023.5` (for Linux, change "osx" to "linux"):

```
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-osx-conda.yml
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-osx-conda.yml
```

#### Snakemake and other dependencies

Third, activate your QIIME 2 environment and install dependencies:

```
conda activate qiime2-2023.5
conda install -c conda-forge -c bioconda biopython muscle clustalo tabulate
conda install -c conda-forge deicode
pip install empress
qiime dev refresh-cache
conda install -c bioconda bioconductor-msa bioconductor-odseq
```

### Option 2: Docker container

An alternative to the native installation is a Docker container. Here are the steps to follow:

#### Docker Desktop

Make sure Docker is installed on your system. If you are using a laptop or desktop machine, you can [install](https://docs.docker.com/get-docker/) Docker Desktop (Mac or Windows) or Docker (Linux). If you are on a compute cluster, you may need to contact your system administrator. This command will list your Docker images and make sure the Docker daemon is running:

```
docker images
```

Make sure Docker has enough memory. On Docker for Mac, the default memory is 2 GB. Go to Preferences -> Resources -> Advanced -> Memory and increase the maximum memory to 8 GB or more if possible.

#### Docker image

Download the [Tourmaline Docker image](https://hub.docker.com/repository/docker/aomlomics/tourmaline) from DockerHub:

```
docker pull aomlomics/tourmaline
```

List your Docker images again to make sure the `tourmaline` image is there:

```
docker images
```

#### Docker container

Now create and run a container:

```
docker run -v $HOME:/data -it aomlomics/tourmaline
```

If installing on a Mac with an Apple M1 chip, run the Docker image with the `--platform linux/amd64` command. It will take a few minutes for the image to load the first time it is run.

```
docker run --platform linux/amd64 -v $HOME:/data -it aomlomics/tourmaline
```

#### External files

The `-v` (volume) flag above allows you to mount a local file system volume (in this case your home directory) to read/write from your container. Note that symbolic links in a mounted volume will not work.

Use mounted volumes to:

* copy metadata and manifest files to your container;
* create symbolic links from your container to your FASTQ files and reference database;
* copy your whole Tourmaline directory out of the container when the run is completed (alternatively, you can clone the Tourmaline directory inside the mounted volume).

To access files that aren't on a mounted volume, from outside the container (running or not), use these commands:

```
docker cp container:source_path destination_path # container to file system
docker cp source_path container:destination_path # file system to container
```

#### Common errors

If you get an error "Plugin error from feature-classifier: Command ... died with <Signals.SIGKILL: 9>.", this means your Docker container ran out of memory. Go to Preferences -> Resources -> Advanced -> Memory and increase the maximum memory to 8 GB or more if possible.

#### Restart a container

You can stop a container by typing `exit`. The container will still be running.

See all your containers with this command:

```
docker ps -a
```

Start a container using its container ID or name with this command:

```
docker start -ia CONTAINER
```

#### Remove containers and images

Remove a Docker container:

```
docker container rm CONTAINER
```

Remove a Docker image:

```
docker rmi IMAGE
```

Remove all stopped containers, all networks not used by at least one container, all images without at least one container associated to them, and all build cache:

```
docker system prune -a
```
