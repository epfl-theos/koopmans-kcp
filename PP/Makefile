# Makefile for Koopmans utils

include	../make.sys

#MODULES = $(PWOBJS) $(QEMODS)

all: merge_evc.x

merge_evc.x: merge_evc.o
	$(LD) $(LDFLAGS) -o $@ merge_evc.o
	- ( cd ../bin ; ln -fs ../PP/$@ . )

#epsilon.x: epsilon.o $(MODULES) $(LIBOBJS)
#	$(LD) $(LDFLAGS) -o $@ epsilon.o $(MODULES) \
#	$(LIBOBJS) $(QELIBS)
#	- ( cd ../bin ; ln -fs ../src/$@ . )

clean:
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L *genmod*

include make.depend
