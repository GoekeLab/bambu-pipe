FROM debian:trixie

# install R 4.4.1 and all CLI tools 
ENV R_VERSION=4.4.1 \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update \ 
&& apt-get install -y wget make g++ git-all zlib1g zlib1g-dev r-base python-is-python3 python3-pip python3-venv gnupg2 libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 

RUN apt-get update -qq && apt-get -y install --no-install-recommends \
    ca-certificates \
    build-essential \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    git-all \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    libopenblas-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libsodium-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

# make a directory to store all required softwares
RUN mkdir mnt/software

# install pigz (v2.8-1)
RUN apt-get update && apt-get install pigz=2.8-1

# install chopper (v0.12.0b)
RUN cd mnt/software && mkdir chopper && cd chopper && wget https://github.com/wdecoster/chopper/releases/download/v0.12.0b/chopper-linux-musl \
&& mv chopper-linux-musl chopper && chmod +x chopper

# install cutadapt (v5.2)
RUN python3 -m venv mnt/software/cutadapt-venv && mnt/software/cutadapt-venv/bin/pip install --no-cache-dir cutadapt==5.2
    
# install flexiplex (v1.02.5)
RUN cd mnt/software && wget https://github.com/DavidsonGroup/flexiplex/archive/refs/tags/v1.02.5.tar.gz \
&& tar -xvf v1.02.5.tar.gz && rm v1.02.5.tar.gz && cd flexiplex-1.02.5 && make
 
# install flexiplex-filter (v1.02.5)
RUN cd mnt/software/flexiplex-1.02.5/scripts && python3 -m venv /mnt/software/flexiplex-filter-venv && /mnt/software/flexiplex-filter-venv/bin/pip install --no-cache-dir .

# install cellranger (v10.0.0) **IMPORTANT: Update the Download Key when building the DockerFile
RUN cd mnt/software && curl -o cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1770734790&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=XJngvA75zSPumb8-fcpq4j4SD31AHtB5Fa0Jjt3phcdFFuhu7UOtz2G6hgUgeOcCRCVhWrZg8Mz-jMYWFyy~XHcB2S7diD3fyedIqhq6ewawhOxPoGOMOGCYrA8rQZrUJHnQD9TE-Bh39g7-7QJfdWG6CXaDI8p6SMb3L4X3dQ-cOl4bRVAS9IVjgRXLpgFdFQ7oS4bHNXlwkpv9S5gPOulHYWd4CDL0bg0DPaoqCqpuMoZabRmbjvcrTeS~Kj-RqFGo4piAbNQpWhbZb3dbyXkOpTKo-jbvtzUsTskRZRn9MZTb8~qO50RTTvHl2~J-tYlQmLCkVMUQr4LEuXEWig__" && tar -xzvf cellranger-10.0.0.tar.gz && rm cellranger-10.0.0.tar.gz

# install minimap2 (v2.30)
RUN cd mnt/software && wget https://github.com/lh3/minimap2/archive/refs/tags/v2.30.tar.gz && tar -xvf v2.30.tar.gz && rm -rf v2.30.tar.gz && cd minimap2-2.30 && make

# install k8 javascript shell (v0.2.4)
RUN cd mnt/software && curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && cp k8-0.2.4/k8-`uname -s` k8

# install samtools (v1.18)
RUN cd mnt/software && wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
&& tar -xvf htslib-1.18.tar.bz2 && rm htslib-1.18.tar.bz2 && cd htslib-1.18 && make && make install \ 
&& wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
&& tar -xvf samtools-1.18.tar.bz2 && rm samtools-1.18.tar.bz2 && cd samtools-1.18 && make && make install

# install jaffa (v2.3)
RUN cd mnt/software && wget https://github.com/Oshlack/JAFFA/releases/download/version-2.3/JAFFA-version-2.3.tar.gz \
 && tar -xvf JAFFA-version-2.3.tar.gz && rm JAFFA-version-2.3.tar.gz && cd JAFFA-version-2.3 $$ ./install_linux64.sh

# install Seurat (v5.4.0)
RUN R -e "install.packages('remotes', repos = 'https://cloud.r-project.org')" 
RUN R -e "remotes::install_version('SeuratObject', version = '5.3.0', repos = 'https://cloud.r-project.org')"
RUN R -e "remotes::install_version('Seurat', version = '5.4.0', repos = 'https://cloud.r-project.org')"

# install bambu
# RUN R -e "install.packages('R.utils', repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages(c('devtools', 'BiocManager'), repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('bambu')"  
RUN cd mnt/software && git clone -b devel_pre_v4 --single-branch https://github.com/GoekeLab/bambu.git
RUN cd mnt/software && R -e "library('devtools'); devtools::load_all('bambu')"

# environment variables 
ENV PATH=$PATH:/mnt/software:/mnt/software/chopper:/mnt/software/cutadapt-venv/bin:/mnt/software/flexiplex-1.02.5:/mnt/software/flexiplex-filter-venv/bin:/mnt/software/minimap2-2.30:/mnt/software/minimap2-2.30/misc:/mnt/software/samtools-1.18:/mnt/software/htslib-1.18:/mnt/software/bambu