FROM fedora:37

RUN dnf -y update \
    && dnf -y install \
        gcc-gfortran \
        gcc-c++ \
        make \
        metis-devel \
        lapack-devel \
        openblas-devel \
        cmake \
        valgrind \
    && dnf clean all

# Build the SuiteSparse libraries for sparse matrix support
# (-k included because of problem with SuiteSparse security certificate - 1 Aug 2021)
RUN curl -kLO http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.1.0.tar.gz \
    && tar -zxvf SuiteSparse-5.1.0.tar.gz \
    && export CXX=/usr/bin/cc \
    && cd SuiteSparse \
    && make install INSTALL=/usr/local BLAS="-L/lib64 -lopenblas"

# Install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/6.1.0.tar.gz \
    && tar -zxvf 6.1.0.tar.gz \
    && cd json-fortran-6.1.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && make install

# copy CVODE source
COPY cvode-3.4-alpha.tar.gz /cvode-3.4-alpha.tar.gz

# Install a modified version of CVODE
RUN tar -zxvf /cvode-3.4-alpha.tar.gz \
    && cd cvode-3.4-alpha \
    && mkdir build \
    && cd build \
    && cmake -D CMAKE_BUILD_TYPE=release \
             -D CMAKE_C_FLAGS_DEBUG="-g -pg" \
             -D CMAKE_EXE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_MODULE_LINKER_FLAGS_DEBUG="-pg" \
             -D CMAKE_SHARED_LINKER_FLAGS_DEBUG="-pg" \
             -D KLU_ENABLE:BOOL=TRUE \
             -D KLU_LIBRARY_DIR=/usr/local/lib \
             -D KLU_INCLUDE_DIR=/usr/local/include \
             .. \
    && make install

# copy CAMP source code and data
COPY . /camp

# Update environment variables
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/lib64:/usr/local/jsonfortran-gnu-6.1.0/lib"
ENV PATH="${PATH}:/usr/local/jsonfortran-gnu-6.1.0/lib"

# Build CAMP
 RUN mkdir build \
    && cd build \
    && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-6.1.0" \
    && cmake -D CMAKE_BUILD_TYPE=release \
             -D CMAKE_C_FLAGS_RELEASE="-g -O3" \
             -D CMAKE_Fortran_FLAGS_RELEASE="-g -O3 -fbacktrace" \
             -D CMAKE_MODULE_LINKER_FLAGS="-g" \
             -D CAMP_ENABLE_MEMCHECK:BOOL=TRUE \
             /camp \
    && make install

WORKDIR /build
