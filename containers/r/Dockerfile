FROM mambaorg/micromamba:1.5.10-noble

USER root
RUN apt-get update && apt-get install -y libcurl4-openssl-dev procps git \
    && rm -rf /var/lib/apt/lists/*

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Install bioconductor version of bambu for its dependencies
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    r-base=4.5.3 \
    r-devtools=2.5.2 \
    r-biocmanager=1.30.27 \
    bioconductor-bambu=3.12.1 \
    r-seurat=5.4.0 \
    gxx=15.2.0 \
    && micromamba clean --all --yes

# Clone bambu single cell feature branch 
RUN git clone --branch exportNewFunctions \
    https://github.com/GoekeLab/bambu.git /opt/bambu \
    && R -e 'devtools::install("/opt/bambu")'

ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"