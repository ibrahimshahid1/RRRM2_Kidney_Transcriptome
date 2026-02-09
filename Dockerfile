# RRRM-2 Kidney Transcriptome Pipeline
# Multi-stage Dockerfile for reproducible cloud execution

FROM rocker/r-ver:4.3.2 AS r-base

# Install R system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    pandoc \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')" && \
    R -e "BiocManager::install(c('DESeq2', 'limma', 'sva', 'edgeR', 'SingleCellExperiment', 'Biobase', 'zellkonverter', 'AnnotationDbi', 'org.Mm.eg.db'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('data.table', 'Matrix', 'ggplot2', 'dplyr'), repos='https://cloud.r-project.org/')"

# Install MuSiC from GitHub
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org/')" && \
    R -e "devtools::install_github('xuranw/MuSiC')"

# Python stage
FROM python:3.11-slim AS python-base

# Install Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Final stage - combine R and Python
FROM rocker/r-ver:4.3.2

# Copy R packages from r-base stage
COPY --from=r-base /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Install Python 3.11
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-venv \
    python3-pip \
    libhdf5-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy Python packages
COPY --from=python-base /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages

# Set up working directory
WORKDIR /app

# Copy repository
COPY . /app

# Install Python package in editable mode
RUN pip3 install -e .

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Default entrypoint
ENTRYPOINT ["python3", "src/run_all_phases.py"]
CMD ["--help"]
