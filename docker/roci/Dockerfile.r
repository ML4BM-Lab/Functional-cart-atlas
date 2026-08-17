# syntax=docker/dockerfile:1

# Rocinante R environment
# Original host: R 4.5.1, Bioconductor 3.21
FROM rocker/r-ver:4.5.1 AS r-builder

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_CACHE_SYMLINKS=FALSE \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    RENV_CONFIG_PAK_ENABLED=FALSE \
    TZ=Etc/UTC

# Build dependencies for the 312-package Rocinante runtime closure. The rocker
# base already supplies the R compiler toolchain, BLAS/LAPACK and common libs.
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        pkg-config \
        pandoc \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff-dev \
        libjpeg-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libcairo2-dev \
        libpango1.0-dev \
        libx11-dev \
        libxt-dev \
        libgraphviz-dev \
        graphviz \
        libgsl-dev \
        libnlopt-dev \
        libsqlite3-dev \
        libicu-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/functional-cart-atlas

COPY environments/roci/Docker/renv.runtime.lock /opt/functional-cart-atlas/renv.lock

# Install into a site library that remains visible when the final container
# starts in /workspace; a project-local renv library would require activation.
RUN R -q -e "install.packages('https://cloud.r-project.org/src/contrib/Archive/renv/renv_1.2.3.tar.gz', repos = NULL, type = 'source')" \
    && R -q -e "renv::restore(lockfile = '/opt/functional-cart-atlas/renv.lock', library = '/usr/local/lib/R/site-library', prompt = FALSE)"


# reticulate embeds Python and therefore needs libpython, which Ubuntu only
# ships in a development package. Isolate that package here and copy just the
# 8.7 MB shared library into the runtime image, without headers or toolchains.
FROM rocker/r-ver:4.5.1 AS python-shared-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        libpython3.12-dev \
    && rm -rf /var/lib/apt/lists/*


FROM rocker/r-ver:4.5.1

LABEL org.opencontainers.image.title="Functional CAR-T Atlas - R (roci)" \
      org.opencontainers.image.description="R runtime for Functional CAR-T Atlas scripts executed on Rocinante" \
      org.opencontainers.image.version="0.2.0"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    RETICULATE_PYTHON=/opt/reticulate-venv/bin/python \
    TZ=Etc/UTC

# Runtime-only additions. Most graphical and numeric libraries already belong
# to rocker/r-ver; Graphviz, GSL and NLopt back compiled packages copied below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-venv \
        pandoc \
        graphviz \
        libgsl27 \
        libharfbuzz-icu0 \
        libnlopt0 \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/reticulate-venv

COPY environments/roci/Docker/requirements.r-reticulate.txt /tmp/requirements.txt

RUN /opt/reticulate-venv/bin/python -m pip install --no-cache-dir --no-deps \
        -r /tmp/requirements.txt \
    && rm -f /tmp/requirements.txt

COPY --from=python-shared-builder \
    /usr/lib/x86_64-linux-gnu/libpython3.12.so.1.0 \
    /usr/lib/x86_64-linux-gnu/libpython3.12.so

COPY --from=r-builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

WORKDIR /workspace

CMD ["R"]
