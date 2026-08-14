INST_DIR = $(HOME)/bin
BIN = do_MC_Para
BIN3 = do_MC_Para_Post_process_v3
BIN1 = $(BIN).o
#FC = /usr/bin/gfortran
#CC = /usr/bin/g++
FFLAGS = -O2 -Wall -ffixed-line-length-none
#CFLAGS = -O2 -Wall -fopenmp -I ~/CODE/boost_1_55_0
#LDLIBS = -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4
FC = gfortran
CC = g++
CFLAGS = -O2 -Wall -fopenmp -std=c++0x
# Point the linker/runtime loader at the active conda env's libgfortran
# (needed when g++ is the system compiler but gfortran comes from conda).
# CONDA_PREFIX is set automatically when the MCMC env is activated; if make
# is run without the env active it is unset/empty, so fall back to the known
# MCMC env path (otherwise -L$(CONDA_PREFIX)/lib collapses to -L/lib and the
# linker fails with "cannot find -lgfortran").
CONDA_PREFIX ?= $(HOME)/anaconda3/envs/MCMC
ifeq ($(strip $(CONDA_PREFIX)),)
CONDA_PREFIX := $(HOME)/anaconda3/envs/MCMC
endif
LDLIBS = -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
INC = -I../boost-1.55.0-gcc630/include



FOBJS = do_MC_Para.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

FOBJS1 = do_MC_Para_Post_process.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

FOBJS1 = do_MC_Para_Post_process_v1.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

FOBJS2 = do_MC_Para_Post_process_v1_M.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

FOBJS3 = do_MC_Para_Post_process_v2.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

FOBJS4 = do_MC_Para_Post_process_v3.o \
	./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o


do_MC_Para_Post_process : do_MC_Para_Post_process.o
	$(CC) $(CFLAGS) $(FOBJS1) -o do_MC_Para_Post_process ${LDLIBS} -fopenmp -lgfortran

do_MC_Para_Post_process_v1 : do_MC_Para_Post_process_v1.o
	$(CC) $(CFLAGS) $(FOBJS1) -o do_MC_Para_Post_process_v1 ${LDLIBS} -fopenmp -lgfortran

do_MC_Para_Post_process_v1_M : do_MC_Para_Post_process_v1_M.o
	$(CC) $(CFLAGS) $(FOBJS2) -o do_MC_Para_Post_process_v1_M ${LDLIBS} -fopenmp -lgfortran

do_MC_Para_Post_process_v2 : do_MC_Para_Post_process_v2.o
	$(CC) $(CFLAGS) $(FOBJS3) -o do_MC_Para_Post_process_v2 ${LDLIBS} -fopenmp -lgfortran

do_MC_Para_Post_process_v3 : do_MC_Para_Post_process_v3.o
	$(CC) $(CFLAGS) $(FOBJS4) -o do_MC_Para_Post_process_v3 ${LDLIBS} -fopenmp -lgfortran

$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN) ${LDLIBS} -fopenmp -lgfortran

do_MC_Para.o : do_MC_Para.C
	$(CC) $(CFLAGS) -c $^ $(INC)

do_MC_Para_Post_process.o : do_MC_Para_Post_process.C
	$(CC) $(CFLAGS) -c $^ $(INC)

do_MC_Para_Post_process_v1.o : do_MC_Para_Post_process_v1.C
	$(CC) $(CFLAGS) -c $^ $(INC)

do_MC_Para_Post_process_v1_M.o : do_MC_Para_Post_process_v1_M.C
	$(CC) $(CFLAGS) -c $^ $(INC)

do_MC_Para_Post_process_v2.o : do_MC_Para_Post_process_v2.C
	$(CC) $(CFLAGS) -c $^ $(INC)

do_MC_Para_Post_process_v3.o : do_MC_Para_Post_process_v3.C
	$(CC) $(CFLAGS) -c $^ $(INC)

clean ::
	rm -f $(BIN) do_MC_Para.o
	rm -f $(BIN3) do_MC_Para_Post_process_v3.o
		
