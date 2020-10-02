# Use the official QIIME 2 image as a parent image
FROM qiime2/core:2020.8

# Label information
LABEL maintainer="Luke Thompson"
LABEL description="Docker image to build the Tourmaline Snakemake workflow for QIIME 2"

# Set the working directory for below commands - necessary? Causes container to open in this directory
#WORKDIR /usr/src/app

# Install dependencies
RUN apt-get update && \
 apt-get install -y \
	snakemake \
	python3-biopython \
	python3-tabulate \
	pandoc && \
 rm -rf /var/lib/apt/lists/* && \
 R -e "install.packages(c('BiocManager'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
 R -e "BiocManager::install(c('msa', 'odseq'))" && \
 pip install git+https://github.com/biocore/empress.git && \
 pip install tabview && \
 qiime dev refresh-cache
