FROM ubuntu:16.04
MAINTAINER Xiao Li

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-7-jdk \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
        pkg-config \
        libgd-dev \
        libperl-dev \
        libgsl0-dev \
        git \
        bwa \
        libgfortran5 \
	libopenblas-dev \
    && rm -rf /var/lib/apt/lists/*

#--------------------------
#cpan
#--------------------------
RUN  wget -O- http://cpanmin.us | perl - -l ~/perl5 App::cpanminus local::lib && \
     eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib` && \
     echo 'eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`' >> ~/.bashrc && \
     echo 'export MANPATH=$HOME/perl5/man:$MANPATH' >> ~/.bashrc && \
     cpanm Statistics::Descriptive && \
     cpanm GD::Graph::histogram

#-----------------------------
# Pipeline components
#-----------------------------

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 && \
    tar -xf htslib-1.8.tar.bz2 && rm htslib-1.8.tar.bz2 && cd htslib-1.8 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
    tar -xf samtools-1.8.tar.bz2 && rm samtools-1.8.tar.bz2 && cd samtools-1.8 && \
    ./configure --with-htslib=/opt/htslib-1.8 && make && make install && make clean

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# Picard tools
RUN mkdir /opt/picard-tools && \
    wget --no-check-certificate -P /opt/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.9.0/picard.jar
#vcftools
RUN apt-get install -y vcftools 
# Pilon
RUN mkdir /opt/pilon && \
    wget --no-check-certificate -P /opt/pilon/ https://github.com/broadinstitute/pilon/releases/download/1.23/pilon.jar

# fasttree
RUN apt-get install -y fasttree

# RAXML
RUN apt-get install -y raxml

# GATK
RUN cd /opt && \
    wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    tar -xf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && rm GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/ GATK-3.8/
    
# FastQC
RUN cd /opt && \
    wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && rm fastqc_v0.11.9.zip && cd FastQC && chmod 755 fastqc && \
    sudo ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc 

#bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar -xf bcftools-1.13.tar.bz2 && rm bcftools-1.13.tar.bz2 && cd bcftools-1.13 && \
    ./configure --enable-libgsl --enable-perl-filters && make && make install && make clean
ENV BCFTOOLS_PLUGINS /opt/bcftools-1.13/plugins:$BCFTOOLS_PLUGINS

#breakdancer
RUN cd /opt && \
    git clone --recursive https://github.com/genome/breakdancer.git && cd breakdancer && mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install && make clean
    
#snpeff
RUN mkdir /opt/snpEff && cd /opt && \
    wget --no-check-certificate https://sourceforge.net/projects/snpeff/files/snpEff_v4_5covid19_core.zip && \
    unzip snpEff_v4_5covid19_core.zip

#readseq
RUN mkdir /opt/readseq && \
    wget --no-check-certificate -P /opt/readseq/ https://sourceforge.net/projects/readseq/files/latest/download/readseq.jar

#gemma
RUN cd /opt && \
    wget --no-check-certificate https://github.com/genetics-statistics/GEMMA/archive/refs/tags/v0.98.5.tar.gz && \
    tar -xf v0.98.5.tar.gz && rm v0.98.5.tar.gz && cd GEMMA-0.98.5/ && make -j 4 && cd bin/ && cp gemma /usr/local/bin/ && \
    make clean

#phyml
RUN apt-get install -y phyml

#plink2
RUN cd /opt && \
    wget --no-check-certificate https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && unzip plink2_linux_x86_64.zip && \
    rm plink2_linux_x86_64.zip && cp plink2 /usr/local/bin/
    
# python modules
RUN pip3 install --upgrade pip
RUN pip3 install tables numpy pandas funpipe

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# scripts
COPY scripts scripts/
