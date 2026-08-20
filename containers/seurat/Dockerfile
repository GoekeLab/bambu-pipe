FROM rocker/r-ver:4.5.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages("pak")'

RUN R -e 'options(timeout = 1800); \
          pak::pkg_install(c("Seurat", "harmony", "R.utils", "bioc::Banksy", \
                             "github::satijalab/seurat-wrappers@9586c29"))'

RUN R -e 'library(Seurat); library(SeuratWrappers); library(Banksy); library(harmony)'