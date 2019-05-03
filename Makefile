CC = g++
DOUBLE_PRECISION = defined
ifdef DOUBLE_PRECISION
PREC = -DUSE_DOUBLE_PRECISION
FFTW = -lfftw3
else
PREC = 
FFTW = -lfftw3f
endif
CFLAGS = -O3 -std=c++0x $(PREC) -DHAVE_TIFF -DHAVE_HDF5 -I/usr/include/hdf5/serial -I/home/optics/Software/fftw/include
#CFLAGS = -O3 -std=c++0x $(PREC) -DHAVE_TIFF -DHAVE_HDF5 -I/usr/include/hdf5/serial
LDFLAGS =  -lm -ltiff -lpng -lhdf5 $(FFTW)
src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
all: $(obj)
	$(CC) -o main  $(obj) $(LDFLAGS)

%.o: %.cpp 
	$(CC) $(CFLAGS) -c -o $@ $< 

.PHONY: clean

clean:
	rm $(obj)
