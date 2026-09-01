INST_DIR = $(HOME)/bin
BIN = do_MC_Para
BIN3 = do_MC_Para_Post_process_v3
BIN1 = $(BIN).o
#FC = /usr/bin/gfortran
#CC = /usr/bin/g++
FFLAGS = -O2 -Wall -ffixed-line-length-none
#CFLAGS = -O2 -Wall -fopenmp -I ~/CODE/boost_1_55_0
#LDLIBS = -L/usr/lib/gcc/x86_64-redhat-linux/4.4.4

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
# ===================== macOS (Apple Silicon / Intel) =====================
# Apple clang ships no OpenMP and no libgfortran, so build with the Homebrew
# GNU toolchain:   brew install gcc
# Override the versions on the command line if yours differ, e.g.
#     make CC=g++-14 FC=gfortran-14
FC = gfortran-15
CC = g++-15
CFLAGS = -O2 -Wall -fopenmp -std=c++0x
FFLAGS = -O2 -ffixed-line-length-none -std=legacy
# libgfortran / libgomp / libstdc++ live in the Homebrew gcc lib dir; the rpath
# lets the binaries run without DYLD_LIBRARY_PATH and without a conda env.
BREW_GCC_LIB := $(shell brew --prefix gcc 2>/dev/null)/lib/gcc/current
LDLIBS = -L$(BREW_GCC_LIB) -Wl,-rpath,$(BREW_GCC_LIB)
else
# ============================== Linux ====================================
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
endif

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

# ---------------------------------------------------------------------------
# Fortran forward-model objects (DISP2/ surface waves, RF/ receiver functions).
# These are shipped in the repo as prebuilt LINUX ELF .o files and the original
# Makefile had no rule to rebuild them.  This rule recompiles them from source,
# which is mandatory on macOS (Linux ELF objects cannot be linked there).
# ---------------------------------------------------------------------------
%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

FOBJS_F = ./DISP2/surfa.o ./DISP2/flat1.o ./DISP2/init.o ./DISP2/calcul.o ./DISP2/fast_surf.o \
	./RF/theo.o ./RF/qlayer.o ./RF/four1.o

ifeq ($(UNAME_S),Darwin)
# ---------------------------------------------------------------------------
# macOS only: normalise inconsistent COMMON-block dimensions AT BUILD TIME.
#
# DISP2/fast_surf.f declares  nper=200  and DISP2/calcul.f declares  ndep=20,
# while init.f / surfa.f declare nper=2000 / ndep=100 for the *same* COMMON
# blocks /o/ and /derivd/.  GNU ld sizes a common block to the LARGEST
# definition, so the inconsistency is harmless on Linux.  The macOS Mach-O
# linker instead takes the SMALLEST one: /o/ gets 32,800 bytes instead of
# 328,000, so init.f's write to t(1) (at offset 160,000) lands outside the
# block and the run dies with "Bus error: 10".
#
# The repo sources are deliberately left UNTOUCHED: a normalised copy is
# generated under macbuild/ and compiled instead, so the Linux box keeps
# building from exactly the same files.
# ---------------------------------------------------------------------------
MACSRC = macbuild

$(MACSRC)/fast_surf.f : DISP2/fast_surf.f
	@mkdir -p $(MACSRC)
	sed 's/nper=200,/nper=2000,/' $< > $@

$(MACSRC)/calcul.f : DISP2/calcul.f
	@mkdir -p $(MACSRC)
	sed 's/ndep=20,/ndep=100,/' $< > $@

./DISP2/fast_surf.o : $(MACSRC)/fast_surf.f
	$(FC) $(FFLAGS) -c $< -o $@

./DISP2/calcul.o : $(MACSRC)/calcul.f
	$(FC) $(FFLAGS) -c $< -o $@
endif

clean ::
	rm -f $(BIN) do_MC_Para.o
	rm -f $(BIN3) do_MC_Para_Post_process_v3.o

# Remove the compiled Fortran objects.  Run this once when moving a checkout
# between Linux and macOS, otherwise the stale objects of the other platform
# are reused and the link fails with "unknown file type".
clean-fortran ::
	rm -f $(FOBJS_F)
	rm -rf macbuild
		
