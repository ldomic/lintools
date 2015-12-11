FROM phusion/baseimage:0.9.17
MAINTAINER “Laura Domicevica <ldomicevica@gmail.com>”

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

ENV CONDA_ENV_PATH /opt/miniconda
ENV MY_CONDA_ENV “my-rdkit-env”
ENV CONDA_ACTIVATE “source $CONDA_ENV_PATH/bin/activate $MY_CONDA_ENV”

RUN apt-get update

RUN apt-get install -y gfortran build-essential make gcc build-essential cmake libboost-all-dev
RUN apt-get install -y libatlas-dev liblapack-dev libgeos-dev
RUN apt-get install -y git-core curl wget nano unzip
RUN apt-get install -y python-matplotlib 

RUN wget --quiet \
    https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    bash Miniconda-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda-latest-Linux-x86_64.sh && \
    chmod -R a+rx /opt/miniconda

ENV PATH $CONDA_ENV_PATH/bin:$PATH
WORKDIR /tmp

RUN conda create -c https://conda.anaconda.org/rdkit -n my-rdkit-env rdkit

RUN conda install -y -n my-rdkit-env pip matplotlib jupyter ipython notebook ipywidgets shapely scipy numpy

RUN conda clean -y -t

#RUN $CONDA_ACTIVATE && pip install MDAnalysis

# Script: Activate virtualenv and launch notebook
ENV STARTSCRIPT /opt/start
RUN source activate my-rdkit-env && pip install MDAnalysis


WORKDIR /data


