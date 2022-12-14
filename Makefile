PREFIX ?= /usr/local

.PHONY: install

default :
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure'
	@echo '  make all'
	@echo '  make install'

all : kcp pp

kcp : bindir mods libs libiotk afclib
	if test -d CPV ; then \
	( cd CPV ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= kcp ; \
	else $(MAKE) $(MFLAGS) TLDEPS= kcp ; fi ) ; fi

pp : bindir mods libs libiotk
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

libiotk :
	if test -d iotk ; then \
	( cd iotk ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= lib+util ; \
	else $(MAKE) $(MFLAGS) TLDEPS= lib+util ; fi ) ; fi

mods : libiotk
	( cd Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )

afclib:
	( if test -d AFC90/src; then cd AFC90/src ; \
        if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= libafc90.a ;\
        else $(MAKE) $(MFLAGS) TLDEPS= libafc90.a ; fi; fi )

libs : mods
	( cd clib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
	( cd flib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )

bindir :
	test -d bin || mkdir bin

# remove object files and executables
clean :
	touch make.sys 
	for dir in \
		CPV PP Modules clib flib iotk AFC90 \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done
	- /bin/rm -rf bin/*.x tmp
	#- cd tests; /bin/rm -rf CRASH *.out *.out2 

# remove configuration files too
distclean veryclean : clean
	- /bin/rm -rf make.sys \
		      config.log configure.msg config.status autom4te.cache \
		      espresso.tar.gz Modules/version.h ChangeLog* \
		      intel.pcl */intel.pcl
	- if test -d GUI ; then \
	( cd GUI ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= veryclean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= veryclean ; fi ) \
	  fi

tar :
	@if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /CVS/ -e /results/ -e'/\.' -e'\.o$$' \
             -e'\.mod$$' -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'\.F90$$' -e'\.x$$' \
	     -e'~$$' -e'\./GUI' | xargs tar rvf espresso.tar
	gzip espresso.tar

links : bindir
	( cd bin/ ; \
	for exe in \
	    ../CPV/kcp.x \
	    ../CPV/cppp.x \
	; do \
	      if test -f $$exe ; then ln -fs $$exe . ; fi \
	done \
	)

depend:
	@echo 'Checking dependencies...'
	- ( if test -x ./makedeps.sh ; then ./makedeps.sh ; fi)

install:
	mkdir -p $(PREFIX)/bin ; \
	for x in `find . -name *.x -type f` ; do \
	cp -v $$x $(PREFIX)/bin/ ; done
	@echo -e '\nkoopmans-kcp binaries are installed in $(PREFIX)/bin\n'

check:
	( cd tests/ ; make run-tests )

# DO NOT DELETE
