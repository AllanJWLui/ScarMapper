# Multi-stage build for ScarMapper
# Stage 1: Builder - compile Python dependencies with C extensions
FROM python:3.11-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libc6-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libmagic-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
WORKDIR /build
COPY pyproject.toml setup.py MANIFEST.in ./
COPY scarmapper/ ./scarmapper/

# Install Python package to /install prefix
RUN pip install --no-cache-dir --prefix=/install .

# Stage 2: Runtime - minimal image with runtime dependencies only
FROM python:3.11-slim

# Install runtime libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libmagic1 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    && rm -rf /var/lib/apt/lists/*

# Copy installed Python packages from builder
COPY --from=builder /install /usr/local

# Copy bundled PEAR binary (statically-linked Linux x86_64 ELF)
COPY Pear/bin/pear /usr/local/bin/pear
RUN chmod +x /usr/local/bin/pear

# Set environment variables for read-only Singularity environments
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV HOME=/tmp

# Create mount points for data and output
RUN mkdir -p /data /output
VOLUME ["/data", "/output"]

WORKDIR /data

# Set entrypoint to scarmapper CLI
ENTRYPOINT ["scarmapper"]
