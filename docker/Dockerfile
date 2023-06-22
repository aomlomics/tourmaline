# Use the official QIIME 2 image as a parent image
FROM quay.io/qiime2/core:2023.2

# Label information
LABEL maintainer="Luke Thompson"
LABEL description="Docker image to build the Tourmaline Snakemake workflow for QIIME 2 v.2023.2"

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
 
RUN pip install empress 

# Add conda installation dir to $PATH (instead of doing 'conda activate')
ENV PATH="/opt/conda/envs/qiime2-2023.2/bin:$PATH"

# This is necessary to install snakemake using conda
SHELL ["conda", "run", "-n", "qiime2-2023.2", "/bin/bash", "-c"]
 
# Install snakemake and other dependencies using conda
RUN conda update -n base -c defaults conda
RUN conda install -c conda-forge -c bioconda snakemake biopython muscle clustalo tabulate
RUN conda install -c conda-forge deicode
RUN qiime dev refresh-cache
RUN conda install -c bioconda bioconductor-msa bioconductor-odseq
