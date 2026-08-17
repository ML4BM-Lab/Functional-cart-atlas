# syntax=docker/dockerfile:1

# Margaret R environment
# Original host: R 4.1.3, Bioconductor 3.14
FROM rocker/r-ver:4.1.3 AS r-builder

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_CACHE_SYMLINKS=FALSE \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    RENV_CONFIG_PAK_ENABLED=FALSE \
    TZ=Etc/UTC

# Compile the audited 309-package mandatory runtime closure in a disposable
# stage. The
# final image receives only the installed R site library, not this toolchain.
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        curl \
        ca-certificates \
        cmake \
        pkg-config \
        pandoc \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        libssh2-1-dev \
        libsodium-dev \
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
        libxext-dev \
        libxrender-dev \
        libxt-dev \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libhdf5-dev \
        libglpk-dev \
        libgsl-dev \
        libgraphviz-dev \
        graphviz \
        libgmp3-dev \
        libmpfr-dev \
        libnlopt-dev \
        libopenblas-dev \
        liblapack-dev \
        libsqlite3-dev \
        libicu-dev \
        libfftw3-dev \
        libmagick++-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libudunits2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# GitHub DESCRIPTION files can name optional repositories in Remotes. The
# audited lock already contains every mandatory dependency, so do not let renv
# expand that exact closure with unlocked optional packages during restore.
ENV RENV_CONFIG_INSTALL_REMOTES=FALSE

WORKDIR /opt/functional-cart-atlas

COPY environments/marg/Docker/renv.runtime.lock /opt/functional-cart-atlas/renv.lock

# Install into the visible site library; a project-local renv library would no
# longer be active after the final container changes to /workspace.
RUN R -q -e "install.packages('https://cloud.r-project.org/src/contrib/Archive/renv/renv_1.2.3.tar.gz', repos = NULL, type = 'source')" \
    && R -q -e "stopifnot(identical(renv::config\$install.remotes(), FALSE))" \
    && R -q -e "renv::restore(lockfile = '/opt/functional-cart-atlas/renv.lock', library = '/usr/local/lib/R/site-library', prompt = FALSE)"


# reticulate needs an unversioned libpython symlink to embed Python. Obtain it
# in isolation so Python headers and development packages do not enter runtime.
FROM rocker/r-ver:4.1.3 AS python-shared-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        libpython3.8-dev \
    && rm -rf /var/lib/apt/lists/*


FROM rocker/r-ver:4.1.3 AS r-runtime

LABEL org.opencontainers.image.title="Functional CAR-T Atlas - R (marg)" \
      org.opencontainers.image.description="R runtime for Functional CAR-T Atlas scripts executed on Margaret" \
      org.opencontainers.image.version="0.2.0"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    RETICULATE_PYTHON=/opt/reticulate-venv/bin/python \
    TZ=Etc/UTC

# Runtime libraries only. A later shared-object audit verifies this list before
# the candidate can receive the final tag.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-distutils \
        python3-venv \
        pandoc \
        graphviz \
        libcurl4 \
        libssl1.1 \
        libxml2 \
        libssh2-1 \
        libsodium23 \
        libharfbuzz-icu0 \
        libhdf5-103 \
        libglpk40 \
        libgsl23 \
        libgmp10 \
        libmpfr6 \
        libnlopt0 \
        libopenblas0 \
        liblapack3 \
        libsqlite3-0 \
        libicu66 \
        libfftw3-double3 \
        libmagick++-6.q16-8 \
        libgdal26 \
        libgeos-c1v5 \
        libproj15 \
        libudunits2-0 \
    && python3 -m venv /opt/reticulate-venv \
    && apt-get purge -y --auto-remove \
        gcc \
        g++ \
        gfortran \
        make \
        cpp \
        gcc-9 \
        g++-9 \
        gfortran-9 \
        binutils \
    && rm -rf /var/lib/apt/lists/*

COPY environments/marg/Docker/requirements.r-reticulate.txt /tmp/requirements.txt

RUN /opt/reticulate-venv/bin/python -m pip install --no-cache-dir --no-deps \
        -r /tmp/requirements.txt \
    && rm -f /tmp/requirements.txt

COPY --from=python-shared-builder \
    /usr/lib/x86_64-linux-gnu/libpython3.8.so.1.0 \
    /usr/lib/libpython3.8.so

COPY --from=r-builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

WORKDIR /workspace

CMD ["R"]


# Flatten the purged runtime into one final layer. Without this stage, Docker
# keeps the deleted compiler bytes in Rocker's lower base layer even though
# they are no longer visible inside the container.
FROM scratch

LABEL org.opencontainers.image.authors="Carl Boettiger <cboettig@ropensci.org>" \
      org.opencontainers.image.base.name="docker.io/library/ubuntu:focal" \
      org.opencontainers.image.licenses="GPL-2.0-or-later" \
      org.opencontainers.image.source="https://github.com/rocker-org/rocker-versioned2" \
      org.opencontainers.image.vendor="Rocker Project" \
      org.opencontainers.image.title="Functional CAR-T Atlas - R (marg)" \
      org.opencontainers.image.description="R runtime for Functional CAR-T Atlas scripts executed on Margaret" \
      org.opencontainers.image.version="0.2.0"

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    R_VERSION=4.1.3 \
    R_HOME=/usr/local/lib/R \
    TZ=Etc/UTC \
    CRAN=https://p3m.dev/cran/__linux__/focal/2022-04-21 \
    LANG=en_US.UTF-8 \
    DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_AUTOLOADER_ENABLED=FALSE \
    RETICULATE_PYTHON=/opt/reticulate-venv/bin/python

COPY --from=r-runtime / /

WORKDIR /workspace

CMD ["R"]
