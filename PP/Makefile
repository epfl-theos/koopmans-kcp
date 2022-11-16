# Makefile for post-processing utilities
include	../make.sys

MODULES = \
../Modules/atom.o \
../Modules/autopilot.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/dspev_drv.o \
../Modules/electrons_base.o \
../Modules/energies.o \
../Modules/error_handler.o \
../Modules/fft_base.o \
../Modules/fft_parallel.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/griddim.o \
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/ions_nose.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/metagga.o \
../Modules/parallel_types.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/radial_grids.o \
../Modules/random_numbers.o \
../Modules/read_upf_v1.o \
../Modules/read_upf_v2.o \
../Modules/read_cards.o \
../Modules/read_namelists.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/sic.o \
../Modules/smallbox.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/twin_types.o \
../Modules/upf_to_internal.o \
../Modules/upf.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/wavefunctions.o \
../Modules/wave_base.o \
../Modules/write_upf_v2.o \
../Modules/xml_input.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o \
../Modules/wannier_new.o \
../Modules/wrappers.o \
../Modules/compute_dipole.o

all: merge_evc.x bin2xml.x

bin2xml.x: bin2xml.o
	$(LD) $(LDFLAGS) -o $@ bin2xml.o $(MODULES) $(LIBOBJS) $(LIBS)
		- ( cd ../bin ; ln -fs ../PP/$@ . )

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
