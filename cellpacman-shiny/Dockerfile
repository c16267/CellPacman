# Serve the PACMAN Shiny application from a container.
#
# Build from the REPOSITORY ROOT (the directory that holds DESCRIPTION,
# cellpacman_<version>.tar.gz and cellpacman-shiny/):
#
#   R CMD build .                                   # only if the tarball is missing/stale
#   docker build -f cellpacman-shiny/Dockerfile -t cellpacman-shiny .
#   docker run --rm -p 3838:3838 -e CELLPACMAN_EXAMPLE_DATA=true cellpacman-shiny
FROM rocker/r-ver:4.4.3

RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
      libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
      libtiff5-dev libjpeg-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

# rocker images pin CRAN to a dated snapshot; the app needs current shiny/bslib.
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"))' \
      >> /usr/local/lib/R/etc/Rprofile.site

# 1. The analysis package (from the source tarball in the repository root).
COPY cellpacman_*.tar.gz /tmp/
RUN R -q -e "install.packages('remotes'); \
             tarball <- list.files('/tmp', pattern = '^cellpacman_.*[.]tar[.]gz$', full.names = TRUE); \
             remotes::install_local(tarball, dependencies = NA, upgrade = 'never')" \
    && rm -f /tmp/cellpacman_*.tar.gz

# 2. The application's own dependencies.
RUN R -q -e "install.packages(c('shiny', 'bslib', 'htmltools', 'reactable', 'shinycssloaders', \
                                'plotly', 'ggplot2', 'dplyr', 'tidyr', 'future', 'promises', \
                                'digest', 'ggpubr', 'ggrepel'))"

# 3. The application itself.
COPY cellpacman-shiny /srv/cellpacman-shiny
WORKDIR /srv/cellpacman-shiny

ENV CELLPACMAN_EXAMPLE_DATA=false \
    CELLPACMAN_WORKERS=1

EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/srv/cellpacman-shiny', host = '0.0.0.0', port = 3838)"]
