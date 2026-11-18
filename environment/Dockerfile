FROM ubuntu:24.04@sha256:77eb101c607d6e3d4f5ad0c8a02c2136c689b98aba1f359e6502ba517e5c57ec

ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8 LC_ALL=C.UTF-8

# Bootstrap TLS trust, then resolve analysis dependencies from the dated snapshot.
RUN apt-get -o Acquire::Retries=3 update -qq \
    && apt-get -o Acquire::Retries=3 install -y --no-install-recommends ca-certificates
COPY ubuntu-snapshot.sources /etc/apt/sources.list.d/ubuntu.sources
RUN apt-get -o Acquire::Retries=3 update -qq \
    && apt-get -o Acquire::Retries=3 install -y --no-install-recommends \
       r-base-core=4.3.3-2build2 r-bioc-deseq2=1.42.0+dfsg-1 r-bioc-limma=3.58.1+dfsg-1 \
    && Rscript --vanilla -e 'stopifnot(getRversion() == "4.3.3", packageVersion("DESeq2") == "1.42.0", packageVersion("limma") == "3.58.1")' \
    && rm -rf /var/lib/apt/lists/*
