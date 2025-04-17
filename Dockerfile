# Dockerfile for ews-nf
FROM bioconductor/r-ver:3.20

# Labels following the Open Containers Initiative (OCI) recommendations
# For more information, see https://specs.opencontainers.org/image-spec/annotations/?v=v1.0.1
LABEL org.opencontainers.image.title="ews-nf"
LABEL org.opencontainers.image.description="Docker image for the ews-nf workflow"
LABEL org.opencontainers.image.source="https://github.com/AlexsLemonade/ews-nf/docker"

# Disable the renv cache to install packages directly into the R library
ENV RENV_CONFIG_CACHE_ENABLED=FALSE

# Install dependencies for renv library
RUN apt-get -y update &&  \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install --no-install-recommends -y \
  pandoc \
  && rm -rf /var/lib/apt/lists/*

# Install renv to enable later package installation
RUN Rscript -e "install.packages('renv')"

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' && \
  rm -rf ~/.cache/R/renv && \
  rm -rf /tmp/downloaded_packages && \
  rm -rf /tmp/Rtmp*

# Install additional dependencies not captured in renv 
RUN R -e "install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')"

# Set CMD to bash to activate the environment for any commands
CMD ["/bin/bash"]
