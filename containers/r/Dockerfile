FROM rocker/r-base:4.4.1

# install system dependencies
RUN apt-get update && apt-get install -y \
	  libcurl4-openssl-dev \
	  procps && rm -rf /var/lib/apt/lists/*

# install Seurat Object (v5.3.0), Seurat (v5.4.0), and Bambu
RUN R -e "install.packages(c('pak', 'BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "pak::pkg_install(c('SeuratObject@5.3.0', 'Seurat@5.4.0', 'GoekeLab/bambu@devel_pre_v4'))"