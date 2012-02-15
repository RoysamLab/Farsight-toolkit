#!/bin/sh
#Repeat scripts 1 & 2 twice
LD_LIBRARY_PATH='/data/research/gcc/lib:/data/research/gcc/libexec/gcc/x86_64-unknown-linux-gnu/4.6.2'
export LD_LIBRARY_PATH
LD_RUN_PATH='/data/research/gcc/lib:/data/research/gcc/libexec/gcc/x86_64-unknown-linux-gnu/4.6.2'
export LD_RUN_PATH
PATH=/data/research/gcc/bin:$PATH
export PATH
#gmp
cd /data
chmod research 755
cd research
wget ftp://ftp.gmplib.org/pub/gmp-5.0.4/gmp-5.0.4.tar.bz2
tar -xvjf gmp-5.0.4.tar.bz2
cd gmp-5.0.4
./configure --prefix=/data/research/gcc --enable-cxx > output.log 2>&1 && make -j69 install >> output.log 2>&1
cd .. && rm gmp-5.0.4.tar.bz2 && rm -fr gmp-5.0.4

#mpfr
wget http://www.mpfr.org/mpfr-current/mpfr-3.1.0.tar.bz2
tar -xvjf mpfr-3.1.0.tar.bz2
cd mpfr-3.1.0
./configure --prefix=/data/research/gcc --with-gmp=/data/research/gcc > output.log 2>&1 && make -j69 install >> output.log 2>&1
cd .. && rm mpfr-3.1.0.tar.bz2 && rm -fr mpfr-3.1.0

#mpc
wget http://www.multiprecision.org/mpc/download/mpc-0.9.tar.gz
tar -xvzf mpc-0.9.tar.gz
cd mpc-0.9
./configure --prefix=/data/research/gcc --with-gmp=/data/research/gcc --with-mpfr=/data/research/gcc > output.log 2>&1 && make -j69 install >> output.log 2>&1
cd .. && rm mpc-0.9.tar.gz && rm -fr mpc-0.9

#ppl
wget ftp://ftp.cs.unipr.it/pub/ppl/releases/0.11.2/ppl-0.11.2.tar.bz2
tar -xvjf ppl-0.11.2.tar.bz2
cd ppl-0.11.2
## in file config on line 11156 insert: MAX_PAIN='-I/data/research/gcc/include'
## change two subsequent variables
## ac_compile='$CXX -c $CXXFLAGS $CPPFLAGS $MAX_PAIN conftest.$ac_ext >&5'
## ac_link='$CXX -o conftest$ac_exeext $CXXFLAGS $CPPFLAGS $MAX_PAIN $LDFLAGS conftest.$ac_ext $LIBS >&5'
#cd /data/research/ppl-0.11.2
#./configure --prefix=/data/research/gcc --with-gmp-prefix=/data/research/gcc > output.log 2>&1 && make -j69 install >> output.log 2>&1
#cd .. && rm ppl-0.11.2.tar.bz2 && rm -fr ppl-0.11.2

##cloog
#wget http://www.bastoul.net/cloog/pages/download/count.php3?url=./cloog-0.17.0.tar.gz
#tar -xvzf cloog-0.17.0.tar.gz
#cd cloog-0.17.0
#./configure --prefix=/data/research/gcc --with-gmp=system --with-gmp-prefix=/data/research/gcc > output.log 2>&1 && make -j69 install >> output.log 2>&1
#cd .. && rm cloog-0.17.0.tar.gz && rm -fr cloog-0.17.0

##gcc
#wget http://www.netgull.com/gcc/releases/gcc-4.6.2/gcc-4.6.2.tar.bz2
#tar -xvjf gcc-4.6.2.tar.bz2
#mkdir gcc-4.6.2-bin
#cd gcc-4.6.2-bin
#/data/research/gcc-4.6.2/configure --enable-threads=posix --prefix=/data/research/gcc --enable-cloog-backend=isl --with-gmp=/data/research/gcc --with-mpfr=/data/research/gcc --with-mpc=/data/research/gcc --with-ppl=/data/research/gcc --with-cloog=/data/research/gcc > output.log 2>&1 && make -j75 >> output.log 2>&1 && make install >> output.log 2>&1
#cd .. && rm gcc-4.6.2.tar.bz2 && rm -fr gcc-4.6.2 gcc-4.6.2-bin
#libtool --finish /data/research/gcc/libexec/gcc/x86_64-unknown-linux-gnu/4.6.2
