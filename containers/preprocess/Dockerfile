FROM mambaorg/micromamba:git-c0f93d2-amazon2023

RUN micromamba install -y -n base \
	-c conda-forge \
	-c bioconda \
	chopper=0.12.0 \
	flexiplex=1.02.5 \
	cutadapt=5.2 \
	pigz=2.8 \
&& micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH