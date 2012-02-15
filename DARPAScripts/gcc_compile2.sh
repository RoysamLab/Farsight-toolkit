#!/bin/sh
## in file config on line 11156 insert: MAX_PAIN='-I/data/research/gcc/include'
## change two subsequent variables
## ac_compile='$CXX -c $CXXFLAGS $CPPFLAGS $MAX_PAIN conftest.$ac_ext >&5'
## ac_link='$CXX -o conftest$ac_exeext $CXXFLAGS $CPPFLAGS $MAX_PAIN $LDFLAGS conftest.$ac_ext $LIBS >&5'
cd /data/research/ppl-0.11.2
./configure --prefix=/data/research/gcc --with-gmp-prefix=/data/research/gcc > output.log 2>&1 && make -j69 install >> output.log 2>&1
cd .. && rm ppl-0.11.2.tar.bz2 && rm -fr ppl-0.11.2

##cloog
wget ftp://gcc.gnu.org/pub/gcc/infrastructure/cloog-0.16.1.tar.gz
tar -xvzf cloog-0.16.1.tar.gz
cd cloog-0.16.1
./configure --prefix=/data/research/gcc --with-gmp=system --with-gmp-prefix=/data/research/gcc > output.log 2>&1 && make -j75 install >> output.log 2>&1
cd .. && rm cloog-0.16.1.tar.gz && rm -fr cloog-0.16.1

##gcc
wget http://www.netgull.com/gcc/releases/gcc-4.6.2/gcc-4.6.2.tar.bz2
tar -xvjf gcc-4.6.2.tar.bz2
mkdir gcc-4.6.2-bin
cd gcc-4.6.2-bin
/data/research/gcc-4.6.2/configure --enable-threads=posix --disable-multilib --disable-shared --prefix=/data/research/gcc --enable-cloog-backend=isl --with-gmp=/data/research/gcc --with-mpfr=/data/research/gcc --with-mpc=/data/research/gcc --with-ppl=/data/research/gcc --with-cloog=/data/research/gcc > output.log 2>&1 && make -j75 >> output.log 2>&1 && make install >> output.log 2>&1
cd .. && rm gcc-4.6.2.tar.bz2 && rm -fr gcc-4.6.2 gcc-4.6.2-bin
libtool --finish /data/research/gcc/libexec/gcc/x86_64-unknown-linux-gnu/4.6.2
