ARG CUDA_VERSION=9.0

#######################
# BEGIN DEVELOP BUILD #
#######################
FROM nvidia/cuda:${CUDA_VERSION}-devel as build

# install libboost & HDF5
RUN /usr/bin/apt update && /usr/bin/apt install -y libboost-all-dev libhdf5-dev libdcmtk-dev libopencv-dev

# install build system + depends
RUN /usr/bin/apt install -y wget build-essential git && \
    /usr/bin/wget https://cmake.org/files/v3.14/cmake-3.14.0-Linux-x86_64.tar.gz && \
    tar -xvf cmake-3.14.0-Linux-x86_64.tar.gz && \
    rm cmake-3.14.0-Linux-x86_64.tar.gz && \
    ln -s /cmake-3.14.0-Linux-x86_64/bin/cmake /usr/local/bin/cmake

# compile dosecalc code
COPY . /src
RUN cd /src && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX="/dosecalc" && make -j8 && make install
COPY ./docker_run.sh /docker_run.sh
WORKDIR /dosecalc
ENV DOSECALC_DATA="/dosecalc/share"
CMD ["bash", "/docker_run.sh"]

#######################
# BEGIN RUNTIME BUILD #
#######################
# FROM nvidia/cuda:${CUDA_VERSION}-base
# COPy --from=build /dosecalc /dosecalc
# RUN /usr/bin/apt install -y libboost-all libhdf5 libdcmtk libopencv
# ENV PATH="/dosecalc/bin:${PATH}"

