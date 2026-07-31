FROM rocker/r-ver:4.5.2

# procps supplies the ps that Nextflow uses to collect task metrics
# git is needed to clone Matrix from Codeberg
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps git \
    && rm -rf /var/lib/apt/lists/*

# pak resolves system requirements automatically, so no hand-maintained apt list
RUN R -e 'install.packages("pak")'

RUN R -e 'options(timeout = 1800); \
          pak::pkg_install(c("Seurat", "harmony", "R.utils", \
                             "github::satijalab/seurat-wrappers"))'

RUN R -e 'library(Seurat); library(SeuratWrappers); library(Banksy); library(harmony)'