FROM continuumio/miniconda3

RUN conda install -c bioconda -c conda-forge \
    samtools gatk freebayes bwa-mem2 trim-galore bcftools python=3.10 click

COPY . /plasmo
WORKDIR /plasmo

ENTRYPOINT ["python", "-m", "plasmo_genomics.cli"]
