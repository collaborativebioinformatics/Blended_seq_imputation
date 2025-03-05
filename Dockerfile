FROM ubuntu:24.04

LABEL org.opencontainers.image.created="2025-03-05"
LABEL org.opencontainers.image.url="https://github.com/odelaneau/GLIMPSE"
LABEL org.opencontainers.image.version="3.0.0"
LABEL org.opencontainers.image.licences="MIT"
LABEL org.opencontainers.image.title="glimpse2"
LABEL org.opencontainers.image.authors="simone.rubinacci@unil.ch"

WORKDIR /docker_build/

# Install required packages
RUN apt-get update -y && \
    apt-get install -y  build-essential curl wget make autoconf automake zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev libbz2-dev unzip git cmake libcurl4-openssl-dev openjdk-21-jdk parallel python3 python3-dev python3-pip libssl-dev libdeflate-dev automake perl zlib1g-dev libperl-dev libgsl0-dev sudo


ENV HTSLIB_CONFIGURE_OPTIONS="--enable-gcs --enable-libcurl"
ENV HTSLIB_LIBRARY_DIR=/usr/local/lib
ENV HTSLIB_INCLUDE_DIR=/usr/local/include
ENV HTSLIB_SOURCE=/htslib-1.21

# Download and build boost program_options and iostreams
RUN wget  https://archives.boost.io/release/1.78.0/source/boost_1_78_0.tar.gz && \
    tar -xf boost_1_78_0.tar.gz && \
    rm boost_1_78_0.tar.gz && \
    cd boost_1_78_0/ && \
    ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
    ./b2 install && \
    cd .. && \
    cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
    cp -r boost/include/boost/ /usr/include/ && \
    rm -r boost_1_78_0 boost

# Download and build htslib, samtools, bcftools, bwa and picard
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \ 
    tar -xf htslib-1.21.tar.bz2 && \
    rm htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    autoreconf -i  && \
    ./configure && \
    make && make install && \
    cd .. && \
    rm -r htslib-1.21 && \ 
    
    wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -xjf bcftools-1.21.tar.bz2 && \
    rm bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters && \
    make && make install && \
    cd .. && \
    rm -r bcftools-1.21 && \
    
    wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjf samtools-1.21.tar.bz2 && \
    rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21  && \
    autoreconf -i  && ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -r samtools-1.21 && \
    
    wget https://github.com/lh3/bwa/archive/refs/tags/v0.7.18.tar.gz && \
    tar -xzf v0.7.18.tar.gz && \
    rm v0.7.18.tar.gz && \
    cd bwa-0.7.18 && make CC='gcc -fcommon' && \
    cp bwa /usr/bin/ && \
    cd .. && \
    rm -r bwa-0.7.18 && \

    wget https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar && \
    mv picard.jar /usr/bin/

# Download and build GLIMPSE
RUN git clone https://github.com/odelaneau/GLIMPSE.git && \
    cd GLIMPSE && \
    make clean && \
    make COMPILATION_ENV=docker

RUN mv GLIMPSE/chunk/bin/GLIMPSE2_chunk GLIMPSE/split_reference/bin/GLIMPSE2_split_reference GLIMPSE/phase/bin/GLIMPSE2_phase GLIMPSE/ligate/bin/GLIMPSE2_ligate GLIMPSE/concordance/bin/GLIMPSE2_concordance /bin && \ 
chmod +x /bin/GLIMPSE2* && \
rm -r GLIMPSE

# nextflow
RUN curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow && \
    mkdir -p $HOME/.local/bin/ && \
    mv nextflow $HOME/.local/bin/

# user
RUN useradd -m -s /bin/bash bumblebee && \
    echo "bumblebee:password" | chpasswd && \
    usermod -aG sudo bumblebee
USER bumblebee

# workspace for blended genome analysis
RUN mkdir -p ~/analysis 
COPY . .

CMD ["/bin/bash"]
