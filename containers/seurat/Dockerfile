FROM rocker/r-ver:4.5.2

RUN apt-get update && apt-get install -y \
        procps libcurl4-openssl-dev libssl-dev zlib1g-dev libglpk-dev libxml2-dev libpng-dev \
        libicu-dev libmagick++-dev gsfonts pandoc python3 git make \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages("remotes")'

RUN R -e 'options(warn = 2, timeout = 1800, repos = c(getOption("repos"), remotes::bioc_install_repos())); \
          remotes::install_version("SeuratObject", "5.3.0", upgrade = "never"); \
          remotes::install_version("Seurat", "5.4.0", upgrade = "never"); \
          remotes::install_version("harmony", "1.2.4", upgrade = "never"); \
          remotes::install_version("R.utils", "2.13.0", upgrade = "never"); \
          remotes::install_version("Banksy", "1.6.0", upgrade = "never"); \
          remotes::install_github("satijalab/seurat-wrappers@9586c29", upgrade = "never")'

RUN R -e 'library(Seurat); library(SeuratWrappers); library(Banksy); library(harmony)'
