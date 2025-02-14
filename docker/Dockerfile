# Use the official QIIME 2 image as a parent image
FROM quay.io/qiime2/amplicon:2024.10

# Label information
LABEL maintainer="Katherine Silliman"
LABEL description="Docker image to build the Tourmaline2 Snakemake workflow for QIIME 2 v.2024.10"

# Set up bash environment: aliases, colors, history
RUN echo "alias cd..='cd ..'" >> ~/.bashrc
RUN echo "alias l='ls --color'" >> ~/.bashrc
RUN echo "alias la='ls -a --color'" >> ~/.bashrc
RUN echo "alias lal='ls -alh --color'" >> ~/.bashrc
RUN echo "alias ll='ls -lh --color'" >> ~/.bashrc
RUN echo "alias rm='rm -i'" >> ~/.bashrc
RUN echo "alias taill='ls -lrt | tail'" >> ~/.bashrc
RUN echo "bind '\"\e[A\": history-search-backward'" >> ~/.bashrc
RUN echo "bind '\"\e[B\": history-search-forward'" >> ~/.bashrc

# Install dependencies using apt-get
RUN apt-get update -y && \
 apt-get install -y \
 build-essential \
 apt-utils \
 nano \
 unzip \
 pandoc && \
 rm -rf /var/lib/apt/lists/*
 

# Add conda installation dir to $PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/bin:${PATH}

# Install snakemake environment
RUN /bin/bash -c "conda update -n base -c conda-forge conda"
RUN /bin/bash -c "conda create -c conda-forge -c bioconda -n snakemake-tour2 snakemake biopython yq parallel"
RUN echo "source activate snakemake-tour2" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}

# Install bowtie2-blca environment
RUN /bin/bash -c "conda create -c conda-forge -c bioconda -n bt2-blca biopython muscle=3.8 bowtie2"
RUN echo "source activate snakemake-tour2" > ~/.bashrc
ENV PATH /opt/conda/envs/snakemake/bin:${PATH}