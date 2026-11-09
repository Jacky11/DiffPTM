# DiffPTM 2 — deployment image (ShinyProxy / standalone)
FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev \
    libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pandoc \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('remotes','BiocManager'), repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c('limma','sva','Biostrings','QFeatures','msqrob2','MSstatsPTM'), ask=FALSE, update=FALSE)" \
 && R -e "install.packages(c('shiny','data.table','ggplot2','cli','openxlsx','officer','rvg','DT','plotly','ggrepel','ggdendro','ggridges','imp4p','cp4p','AhoCorasickTrie','readxl','markdown'), repos='https://cloud.r-project.org')"

COPY . /opt/DiffPTM
RUN R -e "remotes::install_local('/opt/DiffPTM', dependencies = FALSE)"

EXPOSE 3838
CMD ["R", "-e", "options(shiny.port = 3838, shiny.host = '0.0.0.0'); DiffPTM::run_app()"]
