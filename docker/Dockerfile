# Set the base image 
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Leszek Pryszcz <lpryszcz@crg.es>

# miniconda
RUN apt update && apt install -y wget python3-pip git

# dependencies - ont-tombo error due to 
RUN pip install "numpy<1.20" ont_fast5_api matplotlib mappy numba pysam pandas seaborn sklearn

# install the tool
RUN wget -q https://github.com/novoalab/nanoRMS/archive/refs/heads/master.tar.gz && \
    tar xpfz master.tar.gz && rm master.tar.gz && mv nanoRMS-master /opt/nanoRMS && \
    rm -r /opt/nanoRMS/predict_rna_mod /opt/nanoRMS/visualization_per_read

# numpy<1.20 has to be installed before! note, it'll fail with pip > 20.0.2, 19 also fails
RUN pip install git+https://github.com/nanoporetech/tombo.git # this works but with errors for wheel

# clean-up
RUN apt purge -y wget python-pip git && apt autoremove -y && apt autoclean -y && rm -rf /var/lib/apt/lists/*