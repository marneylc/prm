FROM debian:12

LABEL maintaner="Luke Marney <marneyl@oregonstate.edu>"

# dependencies
RUN apt-get update && apt-get -y install \
    software-properties-common \
    htop \
	git \
	libnetcdf-dev \
	libxml2-dev \
	libssl-dev \
	libcurl4-openssl-dev \
	libfontconfig1-dev \
	libharfbuzz-dev \
	libfribidi-dev \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev

# timezone
ENV TZ=America/Vancouver
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install R and packages
RUN apt-get -y install r-base-dev

# Set a stable CRAN mirror
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> /usr/lib/R/etc/Rprofile.site

RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('ggplot2')"
RUN R -e "install.packages('reshape2')"
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('grid')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('lattice')"
RUN R -e "install.packages('tidyr')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('xcms')"

WORKDIR /code

# podman run -it -v ${PWD}:/code prm /usr/bin/R
