FROM rocker/r-ver:4.5.2

# procps supplies the ps that Nextflow uses to collect task metrics
# git is needed to clone Matrix from Codeberg
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps git \
    && rm -rf /var/lib/apt/lists/*

# Refer to BANKSY issue: https://github.com/prabhakarlab/Banksy/issues/76 (removed when fixed)
RUN cd /tmp \
    && git clone --depth 1 https://codeberg.org/jaganmn/Matrix.git \
    && R CMD INSTALL --no-docs Matrix \
    && rm -rf /tmp/Matrix

# pak resolves system requirements automatically, so no hand-maintained apt list
RUN R -e 'install.packages("pak")'

RUN R -e 'options(timeout = 1800); \
          pak::pkg_install(c("Seurat", "harmony", "R.utils", \
                             "github::satijalab/seurat-wrappers"))'

RUN R -e 'library(Seurat); library(SeuratWrappers); library(Banksy); library(harmony); \
          stopifnot(packageVersion("Matrix") >= "1.8-0")'