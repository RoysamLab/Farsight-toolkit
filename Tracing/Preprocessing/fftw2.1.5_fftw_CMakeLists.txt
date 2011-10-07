
FILE(GLOB SourceFiles *.c)
FILE(GLOB HeaderFiles *.h)


ADD_LIBRARY(libfftw215 ${SourceFiles} ${HeaderFiles})
INSTALL(TARGETS libfftw215 ARCHIVE DESTINATION lib)
INSTALL(FILES fftw.h DESTINATION include)
