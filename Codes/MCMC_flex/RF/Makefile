INST_DIR = $(HOME)/bin
BIN = 

FC = gfortran
CC = g++
#FC = g77
#CC = gcc

FFLAGS = -O2 -Wall -ffixed-line-length-none
CFLAGS = -O2 -Wall

LDLIBS = -L/uufs/chpc.utah.edu/sys/installdir/fftw/3.3.6/lib -lfftw3 #emberg
#LDLIBS =

FOBJS = four1.o qlayer.o

theo.o: theo.f $(FOBJS)
	$(FC) $(FFLAGS) -c -o theo.o theo.f $(FOBJS) 

clean : 
	rm -f $(BIN) core $(FOBJS)

