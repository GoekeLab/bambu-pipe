# Seurat 5.4.0 and harmony are already installed in the base image
FROM satijalab/seurat:5.4.0

# procps supplies the ps that Nextflow uses to collect task metrics
# libmagick++-dev is needed by magick, an import of SpatialExperiment (a Banksy dependency)
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Seurat Wrappers, plus Banksy for spatially aware clustering (Banksy is only a
# suggested dependency of Seurat Wrappers, so it is not installed alongside it)
RUN R -e 'options(timeout = 600); remotes::install_github("satijalab/seurat-wrappers", upgrade = "always")' \
    && R -e 'library(SeuratWrappers); library(Banksy); library(harmony)'
