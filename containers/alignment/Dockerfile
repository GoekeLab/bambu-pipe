FROM mambaorg/micromamba:git-c0f93d2-amazon2023

RUN micromamba install -y -n base \
	-c conda-forge \
	-c bioconda \
	minimap2=2.30 \
	samtools=1.23 \
&& micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH