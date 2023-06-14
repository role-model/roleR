#!/bin/bash

apt update; apt install libgsl-dev libpng-dev libxml2-dev libglpk-dev libfontconfig1-dev libharfbuzz-dev  libfribidi-dev libtiff-dev -y
# libpng-dev is for 'remotes'
# libgsl-dev is for 'msprime' (I think)
# libxml2-dev is for 'taxize' <- RoLE workshop Part I
# libglpk-dev is for hillR
# libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff-dev for tidyverse

sudo -u rstudio bash -i -c '\
    cd ~
    pwd
    if uname -a | grep aarch64 >> /dev/null;
    then
        echo "Detected arm"
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh
        bash Miniconda3-latest-Linux-aarch64.sh -b
    elif uname -a | grep x86_64 >> /dev/null;
    then
        echo "Detected x86_64"
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b
    else
        echo "Failed to detect host architecture. Must be arm (Mac) or x86_64 (Linux/Windows)"
    fi
    ./miniconda3/bin/conda init
'

sudo -u rstudio bash -i -c '\
    cd ~
    which python
    conda install -c bioconda python-newick -y
    conda install -c conda-forge pandas matplotlib scikit-learn toytree ipyparallel -y
    rm -rf msprime
    git clone --recurse-submodules https://github.com/tskit-dev/msprime.git
    cd msprime; pip install .
    mkdir src
    cd src
    git clone https://github.com/iBioGen/iBioGen.git
    cd iBioGen
    pip install -e .
'

echo "install.packages(\"remotes\")
install.packages(\"reticulate\") # <- Necessary for installing/configuring conda base
install.packages(c(\"tidyverse\", \"ape\", \"taxize\", \"hillR\", \"spocc\", \"rotl\", \"rentrez\", \"vcfR\")) # <- Workshop Part I reqs

if (!require(\"BiocManager\", quietly = TRUE)) # <- msa is special and comes from BiocManager
  install.packages(\"BiocManager\")
BiocManager::install(\"msa\")

install.packages(\"tidymodels\") <- For part II inference
install.packages(\"shinyWidgets\") <- For RoLE-Shiny
library(reticulate)
library(remotes)
reticulate::use_condaenv('base') # <- Tell reticulate to use the externally installed conda env
# NB: use reticulate::conda_list() to check the env name if this is not correct, but it should be by default
remotes::install_github('role-model/roleR', force=TRUE)
remotes::install_github(\"role-model/roleShiny\")" > tmp.txt

sudo -u rstudio bash -i -c 'Rscript tmp.txt'

echo "library(reticulate)
reticulate::use_condaenv('base')" > /home/rstudio/.Rprofile
chown rstudio:rstudio /home/rstudio/.Rprofile
