FROM debian:11

LABEL maintaner="Luke Marney <marneyl@oregonstate.edu>"

# dependencies
RUN apt-get update 
RUN apt-get -y install \
    software-properties-common \
    htop \
	git \
	libnetcdf-dev \
	libxml2-dev \
	libssl-dev

# timezone
ENV TZ=America/Vancouver
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# install R and packages
RUN apt-get -y install r-base-dev
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

# docker run -it -v ${PWD}:/code debianr_xcms /usr/bin/R
