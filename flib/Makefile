# Makefile for flib

include ../make.sys

OBJS = \
avrec.o \
atomic_number.o \
bachel.o \
capital.o \
dost.o \
erf.o \
functionals.o \
lsda_functionals.o \
more_functionals.o \
iglocal.o \
inpfile.o \
int_to_char.o \
invmat.o \
invmat_complex.o \
latgen.o \
linpack.o \
matches.o \
recips.o \
remove_tot_torque.o \
simpsn.o \
sort.o \
sph_bes.o \
sph_dbes.o \
transto.o \
date_and_tim.o \
sort_gvec.o \
volume.o \
dylmr2.o \
ylmr2.o \
cryst_to_car.o

POBJS = \
flush_unit.o \
gridsetup.o \
localdim.o \
localindex.o

all: flib.a ptools.a

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.a *.o *.mod *.i *.F90 core* *.L

include make.depend
