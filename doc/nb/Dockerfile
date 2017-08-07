FROM ubuntu:16.04

MAINTAINER Theodore Kisner <tskisner@lbl.gov>

# Use bash

SHELL ["/bin/bash", "-c"]

# Install system dependencies.

RUN apt-get update \
    && apt-get install -y curl procps build-essential gfortran git subversion \
    python libcairo2-dev libpixman-1-dev libgsl-dev flex pkg-config cmake \
    autoconf m4 libtool automake locales \
    && rm -fr /var/lib/apt/lists/*

# Set up locales, to workaround a pip bug

RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8 
ENV LANGUAGE en_US:en 
ENV LC_ALL en_US.UTF-8

# We install everything directly into /usr so that we do
# not need to modify the default library and executable
# search paths.  Shifter will manipulate LD_LIBRARY_PATH,
# so it is important not to use that variable.

# Create working directory for builds

RUN mkdir /usr/src/desi
WORKDIR /usr/src/desi

# Install conda root environment

ENV PYTHONPATH ""
ENV PYTHONSTARTUP ""
ENV PYTHONNOUSERSITE "1"
ENV PYTHONUSERBASE "/tmp"

RUN curl -SL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    -o miniconda.sh \
    && /bin/bash miniconda.sh -b -f -p /usr \
    && conda config --add channels intel \
    && conda config --remove channels intel \
    && conda install --copy --yes python=3.5 \
    && rm miniconda.sh \
    && rm -rf /usr/pkgs/*

ENV ACCEPT_INTEL_PYTHON_EULA yes

RUN conda config --add channels intel \
    && conda install --copy --yes $(echo intelpython3.5_core | sed -e "s#\.._#_#") \
    && rm -rf /usr/pkgs/*

# Install conda packages.

RUN conda install --copy --yes \
    nose \
    requests \
    future \
    cython \
    numpy \
    scipy \
    matplotlib \
    basemap \
    seaborn \
    pyyaml \
    astropy \
    hdf5 \
    h5py \
    psutil \
    ephem \
    psycopg2 \
    numba=0.28.1 \
    scikit-learn \
    scikit-image \
    && conda install --copy --yes -c defaults \
    ipython ipython-notebook \
    && rm -rf /usr/pkgs/*

# Install pip packages.

RUN pip install --no-binary :all: \
    speclite \
    hpsspy \
    photutils \
    https://github.com/esheldon/fitsio/archive/v0.9.12rc1.zip

# Install MKL headers

RUN curl -SL https://www.dropbox.com/s/pzps5m45mnmv6ap/mkl_2017.0.3_include.tar.gz?dl=1 \
    -o mkl_2017.0.3_include.tar.gz \
    && tar xzf mkl_2017.0.3_include.tar.gz -C /usr/include --strip=1 \
    && rm mkl_2017.0.3_include.tar.gz \
    && if [ "/usr" != "/usr" ]; then \
	for lb in $(ls /usr/lib/libmkl*); do \
	    lbname=$(basename ${lb}); \
	    linkpath="/usr/lib/${lbname}"; \
	    if [ ! -e ${linkpath} ]; then \
	        ln -s "${lb}" "${linkpath}"; \
	    fi; \
	done; \
    fi

# The conda TCL packages overwrite the system-installed regex.h.  So
# now we force reinstall of the package that provides that

RUN apt-get update \
    && apt-get install -y --reinstall libc6-dev \
    && rm -fr /var/lib/apt/lists/*

# Copy all patch files to current working directory

RUN mkdir ./rules
ADD rules/patch_* ./rules/

# Install MPICH 3.2 which is compatible with the external
# Cray MPICH which is prepended to LD_LIBRARY_PATH as part
# of shifter.

RUN curl -SL http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz \
    -o mpich-3.2.tar.gz \
    && tar -xzf mpich-3.2.tar.gz \
    && cd mpich-3.2 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" ./configure  --prefix=/usr \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf mpich-3.2*

# Install mpi4py.

RUN curl -SL https://pypi.python.org/packages/ee/b8/f443e1de0b6495479fc73c5863b7b5272a4ece5122e3589db6cd3bb57eeb/mpi4py-2.0.0.tar.gz#md5=4f7d8126d7367c239fd67615680990e3 \
    -o mpi4py-2.0.0.tar.gz \
    && tar xzf mpi4py-2.0.0.tar.gz \
    && cd mpi4py-2.0.0 \
    && python setup.py build --mpicc="mpicc" --mpicxx="mpicxx" \
    && python setup.py install --prefix=/usr \
    && cd .. \
    && rm -rf mpi4py*

# Install CFITSIO.

RUN curl -SL http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3410.tar.gz \
    -o cfitsio3410.tar.gz \
    && tar xzf cfitsio3410.tar.gz \
    && cd cfitsio \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" ./configure  --prefix=/usr \
    && make -j 4 && make shared \
    && make fpack && make funpack \
    && make install \
    && cd .. \
    && rm -rf cfitsio*

RUN curl -SL https://heasarc.gsfc.nasa.gov/docs/software/ftools/fitsverify/fitsverify-4.18.tar.gz \
    -o fitsverify-4.18.tar.gz \
    && tar xzf fitsverify-4.18.tar.gz \
    && cd fitsverify \
    && gcc -O3 -fPIC -pthread -I/usr/include -DSTANDALONE -o fitsverify ftverify.c \
    fvrf_data.c fvrf_file.c fvrf_head.c fvrf_key.c fvrf_misc.c \
    -L/usr/lib -lcfitsio -lm \
    && cp -a fitsverify "/usr/bin/" \
    && cd .. \
    && rm -rf fitsverify*

# Install FFTW.

RUN curl -SL http://www.fftw.org/fftw-3.3.5.tar.gz \
    -o fftw-3.3.5.tar.gz \
    && tar xzf fftw-3.3.5.tar.gz \
    && cd fftw-3.3.5 \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" ./configure --enable-threads  --prefix="/usr" \
    && make -j 4 && make install \
    && make clean \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" ./configure --enable-float --enable-threads  --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf fftw*

# Install BOOST.

RUN curl -SL https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.bz2/download \
    -o boost_1_63_0.tar.bz2 \
    && tar xjf boost_1_63_0.tar.bz2 \
    && cd boost_1_63_0 \
    && echo "" > tools/build/user-config.jam \
    && echo 'using gcc : : g++ : <cflags>"-O3 -fPIC -pthread" <cxxflags>"-O3 -fPIC -pthread" ;' >> tools/build/user-config.jam \
    && echo 'using mpi : mpicxx : <include>"/usr/include" <library-path>"/usr/lib" <find-shared-library>"mpichcxx" <find-shared-library>"mpich" ;' >> tools/build/user-config.jam \
    && BOOST_BUILD_USER_CONFIG=tools/build/user-config.jam ./bootstrap.sh \
    --with-toolset=gcc \
    --with-python=python3.5 \
    --prefix=/usr \
    && ./b2 --layout=tagged \
    $(python3.5-config --includes | sed -e 's/-I//g' -e 's/\([^[:space:]]\+\)/ include=\1/g') \
    variant=release threading=multi link=shared runtime-link=shared install \
    && cd .. \
    && rm -rf boost*

# Install healpix

RUN curl -SL https://github.com/tskisner/healpix-autotools/releases/download/v3.31.2/healpix-autotools-3.31.2.tar.gz \
    -o healpix-autotools-3.31.2.tar.gz \
    && tar -xzf healpix-autotools-3.31.2.tar.gz \
    && cd healpix-autotools-3.31.2 \
    && CC="gcc" CXX="g++" FC="gfortran" \
    CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" FCFLAGS="-O3 -fPIC -fexceptions -pthread" \
    ./configure  --with-cfitsio="/usr" --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf healpix*

# GalSim and dependencies

RUN curl -SL http://prdownloads.sourceforge.net/scons/scons-2.5.1.tar.gz \
    -o scons-2.5.1.tar.gz \
    && tar xzf scons-2.5.1.tar.gz \
    && cd scons-2.5.1 \
    && python2 setup.py install --prefix="/usr" \
    && cd .. \
    && rm -rf scons* \
    && perl -i -p -e 's/usr\/bin\/env python/usr\/bin\/env python2/g' /usr/bin/scons*

RUN curl -SL https://github.com/rmjarvis/tmv/archive/v0.74.tar.gz \
    -o tmv-0.74.tar.gz \
    && tar xzf tmv-0.74.tar.gz \
    && cd tmv-0.74 \
    && scons PREFIX="/usr" \
    CXX="g++" FLAGS="-O3 -fPIC -pthread" \
    EXTRA_INCLUDE_PATH="/usr/include" \
    EXTRA_LIB_PATH=$(echo -L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl | sed -e 's#-L\([^ ]*\).*#\1#') \
    LIBS="$(echo -L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl | sed -e 's#-L\([^ ]*\)##' | sed -e 's#-l\([^ ]*\)#\1#g' | sed -e 's#-[^ ]*##g')" \
    FORCE_FBLAS=true \
    && scons PREFIX="/usr" install \
    && cd .. \
    && rm -rf tmv*


RUN curl -SL https://github.com/GalSim-developers/GalSim/archive/v1.4.3.tar.gz \
    -o GalSim-1.4.3.tar.gz \
    && tar xzf GalSim-1.4.3.tar.gz \
    && cd GalSim-1.4.3 \
    && scons PREFIX="/usr" \
    CXX="g++" FLAGS="-O3 -fPIC -pthread -std=c++98" \
    TMV_DIR="/usr" \
    FFTW_DIR="/usr" \
    BOOST_DIR="/usr" \
    TMV_DIR="/usr" \
    PYTHON="/usr/bin/python" \
    EXTRA_INCLUDE_PATH="/usr/include" \
    EXTRA_LIB_PATH="/usr/lib:/usr/lib" \
    && scons PREFIX="/usr" install \
    && cd .. \
    && rm -rf GalSim*


# Image libraries needed by astrometry

RUN curl -SL http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz \
    -o bzip2-1.0.6.tar.gz \
    && tar xzf bzip2-1.0.6.tar.gz \
    && cd bzip2-1.0.6 \
    && make CC="gcc" CFLAGS="-O3 -fPIC -pthread" PREFIX="/usr" \
    && make CC="gcc" CFLAGS="-O3 -fPIC -pthread" PREFIX="/usr" install \
    && cd .. \
    && rm -rf bzip2*

RUN curl -SL http://zlib.net/zlib-1.2.11.tar.gz \
    -o zlib-1.2.11.tar.gz \
    && tar xzf zlib-1.2.11.tar.gz \
    && cd zlib-1.2.11 \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf zlib*

RUN curl -SL http://prdownloads.sourceforge.net/libpng/libpng-1.6.29.tar.gz?download \
    -o libpng-1.6.29.tar.gz \
    && tar xzf libpng-1.6.29.tar.gz \
    && cd libpng-1.6.29 \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" LDFLAGS="-L/usr/lib -lz" \
    ./configure  \
    --with-zlib-prefix="/usr" \
    --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf libpng*

RUN curl -SL http://sourceforge.net/projects/libjpeg/files/libjpeg/6b/jpegsrc.v6b.tar.gz/download \
    -o jpegsrc.v6b.tar.gz \
    && tar xzf jpegsrc.v6b.tar.gz \
    && cd jpeg-6b \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix="/usr" \
    && make -j 4 \
    && mkdir -p "/usr/man/man1" \
    && make install \
    && make install-lib \
    && make install-headers \
    && cd .. \
    && rm -rf jpeg*


RUN echo hello

RUN curl -SL http://download.osgeo.org/libtiff/tiff-4.0.7.tar.gz \
    -o tiff-4.0.7.tar.gz \
    && tar xzf tiff-4.0.7.tar.gz \
    && cd tiff-4.0.7 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf tiff*

# wcslib

#RUN curl  -SL ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib-5.16.tar.bz2 \
RUN curl  -SL https://launchpad.net/ubuntu/+archive/primary/+files/wcslib_5.16.orig.tar.bz2 \
    -o wcslib-5.16.tar.bz2 \
    && tar xjf wcslib-5.16.tar.bz2 \
    && cd wcslib-5.16 \
    && chmod -R u+w . \
    && patch -p1 < ../rules/patch_wcslib \
    && autoconf \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" \
    CPPFLAGS="-I/usr/include" \
    LDFLAGS="-L/usr/lib" \
    ./configure  \
    --disable-fortran \
    --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf wcslib*

# pcre

RUN curl -SL https://sourceforge.net/projects/pcre/files/pcre/8.40/pcre-8.40.tar.bz2/download \
    -o pcre-8.40.tar.bz2 \
    && tar xjf pcre-8.40.tar.bz2 \
    && cd pcre-8.40 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf pcre*

# swig

RUN curl -SL http://prdownloads.sourceforge.net/swig/swig-3.0.12.tar.gz \
    -o swig-3.0.12.tar.gz \
    && tar xzf swig-3.0.12.tar.gz \
    && cd swig-3.0.12 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix=/usr \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf swig*

# ceres

RUN curl -SL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz \
    -o metis-5.1.0.tar.gz \
    && tar xzf metis-5.1.0.tar.gz \
    && cd metis-5.1.0 && patch -p1 < ../rules/patch_metis \
    && cd .. \
    && mkdir metis_build && cd metis_build \
    && cmake \
    -D CMAKE_C_COMPILER="gcc" \
    -D CMAKE_C_FLAGS="-O3 -fPIC -pthread" \
    -D GKLIB_PATH="../metis-5.1.0/GKlib" \
    -D CMAKE_INSTALL_PREFIX="/usr" \
    ../metis-5.1.0 \
    && make \
    && make install \
    && cd .. \
    && rm -rf metis*

RUN curl -SL https://github.com/google/glog/archive/v0.3.4.tar.gz \
    -o glog-0.3.4.tar.gz \
    && tar xzf glog-0.3.4.tar.gz \
    && cd glog-0.3.4 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" ./configure  \
    --prefix="/usr" \
    && make -j 4 && make install \
    && cd .. \
    && rm -rf glog*

RUN curl -SL http://bitbucket.org/eigen/eigen/get/3.3.3.tar.bz2 \
    -o eigen-eigen-67e894c6cd8f.tar.bz2 \
    && tar xjf eigen-eigen-67e894c6cd8f.tar.bz2 \
    && mkdir eigen_build && cd eigen_build \
    && cmake \
    -D CMAKE_C_COMPILER="gcc" \
    -D CMAKE_CXX_COMPILER="g++" \
    -D CMAKE_C_FLAGS="-O3 -fPIC -pthread" \
    -D CMAKE_CXX_FLAGS="-O3 -fPIC -pthread" \
    -D CMAKE_INSTALL_PREFIX="/usr" \
    ../eigen-eigen-67e894c6cd8f \
    && make \
    && make install \
    && cd .. \
    && rm -rf eigen*

RUN curl -SL http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.5.4.tar.gz \
    -o SuiteSparse-4.5.4.tar.gz \
    && tar xzf SuiteSparse-4.5.4.tar.gz \
    && cd SuiteSparse \
    && patch -p1 < ../rules/patch_suitesparse \
    && make CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" AUTOCC=no \
    F77="gfortran" F77FLAGS="@FCFLAGS" \
    CFOPENMP="-fopenmp" LAPACK=" " BLAS="-L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl" \
    MY_METIS_INC="/usr/include" MY_METIS_LIB="-L/usr/lib -lmetis" \
    && make install CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" AUTOCC=no \
    F77="gfortran" F77FLAGS="@FCFLAGS" \
    CFOPENMP="-fopenmp" LAPACK=" " BLAS="-L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl" \
    MY_METIS_INC="/usr/include" MY_METIS_LIB="-L/usr/lib -lmetis" \
    INSTALL="/usr" \
    && cd .. \
    && rm -rf SuiteSparse*

RUN curl -SL http://ceres-solver.org/ceres-solver-1.12.0.tar.gz \
    -o ceres-solver-1.12.0.tar.gz \
    && tar xzf ceres-solver-1.12.0.tar.gz \
    && mkdir ceres_build && cd ceres_build \
    && cmake \
    -D CMAKE_C_COMPILER="gcc" \
    -D CMAKE_CXX_COMPILER="g++" \
    -D CMAKE_C_FLAGS="-O3 -fPIC -pthread" \
    -D CMAKE_CXX_FLAGS="-O3 -fPIC -pthread" \
    -D EIGEN_INCLUDE_DIR="/usr/include/eigen3" \
    -D SUITESPARSE_INCLUDE_DIR_HINTS="/usr/include" \
    -D SUITESPARSE_LIBRARY_DIR_HINTS="/usr/lib" \
    -D GLOG_INCLUDE_DIR="/usr/include" \
    -D GLOG_LIBRARY="/usr/lib/libglog.so" \
    -D BUILD_SHARED_LIBS=ON \
    -D BUILD_TESTING=FALSE \
    -D BUILD_EXAMPLES=FALSE \
    -D CMAKE_INSTALL_PREFIX="/usr" \
    ../ceres-solver-1.12.0 \
    && make \
    && make install \
    && cd .. \
    && rm -rf ceres*

# astrometry.net
# This package and tractor install python tools to a "lib/python"
# directory rather than site-packages.  Make a symlink here.

RUN ln -s "/usr/lib/python3.5/site-packages" "/usr/lib/python"

RUN curl -SL https://github.com/dstndstn/astrometry.net/releases/download/0.72/astrometry.net-0.72.tar.gz \
    -o astrometry.net-0.72.tar.gz \
    && tar xzf astrometry.net-0.72.tar.gz \
    && cd astrometry.net-0.72 \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" \
    LDFLAGS="-L/usr/lib -lz" make \
    WCSLIB_INC="-I/usr/include/wcslib" WCSLIB_LIB="-L/usr/lib -lwcs" \
    JPEG_INC="-I/usr/include" JPEG_LIB="-L/usr/lib -ljpeg" \
    CFITS_INC="-I/usr/include" CFITS_LIB="-L/usr/lib -lcfitsio -lm" \
    && CC="gcc" CXX="g++" CFLAGS="-O3 -fPIC -pthread" CXXFLAGS="-O3 -fPIC -pthread" \
    LDFLAGS="-L/usr/lib -lz" make \
    WCSLIB_INC="-I/usr/include/wcslib" WCSLIB_LIB="-L/usr/lib -lwcs" \
    JPEG_INC="-I/usr/include" JPEG_LIB="-L/usr/lib -ljpeg" \
    CFITS_INC="-I/usr/include" CFITS_LIB="-L/usr/lib -lcfitsio -lm" \
    extra \
    && make install INSTALL_DIR="/usr" \
    && cd .. \
    && rm -rf astrometry*

# tractor

RUN curl -SL https://github.com/dstndstn/tractor/archive/dr5.2.tar.gz \
    -o tractor-dr5.2.tar.gz \
    && tar xzf tractor-dr5.2.tar.gz \
    && cd tractor-dr5.2 \
    && CERES_LIB="-L/usr/lib -lceres" \
    python setup.py install --with-ceres --prefix="/usr" \
    && cd .. \
    && rm -rf tractor*

# Astromatic toolkit pieces

RUN curl -SL http://www.astromatic.net/download/sextractor/sextractor-2.19.5.tar.gz \
    -o sextractor-2.19.5.tar.gz \
    && tar xzf sextractor-2.19.5.tar.gz \
    && cd sextractor-2.19.5 \
    && patch -p1 < ../rules/patch_sextractor \
    && chmod +x autogen.sh \
    && ./autogen.sh \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" CPPFLAGS="-I/usr/include" ./configure  \
    --with-lapacke=" -L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl" --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf sextractor*

RUN curl -SL http://www.astromatic.net/download/psfex/psfex-3.17.1.tar.gz \
    -o psfex-3.17.1.tar.gz \
    && tar xzf psfex-3.17.1.tar.gz \
    && cd psfex-3.17.1 \
    && patch -p1 < ../rules/patch_psfex \
    && chmod +x autogen.sh \
    && ./autogen.sh \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" CPPFLAGS="-I/usr/include" ./configure  \
    --with-lapacke=" -L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl" --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf psfex*

RUN curl -SL http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient-3.84.tar.gz \
    -o cdsclient-3.84.tar.gz \
    && tar xzf cdsclient-3.84.tar.gz \
    && cd cdsclient-3.84 \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" \
    ./configure  \
    --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf cdsclient*

RUN curl -SL http://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz \
    -o scamp-2.0.4.tar.gz \
    && tar xzf scamp-2.0.4.tar.gz \
    && cd scamp-2.0.4 \
    && patch -p1 < ../rules/patch_scamp \
    && chmod +x autogen.sh \
    && ./autogen.sh \
    && CC="gcc" CFLAGS="-O3 -fPIC -pthread" CPPFLAGS="-I/usr/include" ./configure  \
    --with-lapacke=" -L/usr/lib -lmkl_rt -fopenmp -lpthread -lm -ldl" --prefix="/usr" \
    && make && make install \
    && cd .. \
    && rm -rf scamp*

# Remove the temporary patches

RUN rm -rf ./rules

# Precompile all python modules.  Ignore errors.

RUN python -m compileall -f "/usr/lib/python3.5/site-packages"; exit 0

# Create a fake home directory so that packages can create
# config files if needed

RUN mkdir /home/desi
RUN mkdir /home/desi/.astropy

WORKDIR /home/desi
ENV HOME /home/desi

RUN python -c "import astropy"
RUN python -c "import matplotlib.font_manager as fm; f = fm.FontManager"

# Set the entrypoint and default command

ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

