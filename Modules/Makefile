# Makefile for Modules

include ../make.sys

MODULES = \
atom.o \
autopilot.o \
berry_phase.o \
bfgs_module.o \
cell_base.o  \
check_stop.o  \
clocks.o \
constants.o \
constraints_module.o \
control_flags.o \
descriptors.o \
dspev_drv.o \
electrons_base.o \
energies.o \
error_handler.o \
exc_t.o \
fft_base.o \
fft_parallel.o \
fft_scalar.o \
fft_types.o \
functionals.o \
griddim.o \
input_parameters.o \
io_files.o \
io_global.o  \
ions_base.o \
ions_nose.o \
kind.o \
metagga.o \
mm_dispersion.o \
mp.o \
mp_base.o \
mp_global.o \
mp_wave.o \
path_formats.o \
path_io_routines.o  \
path_opt_routines.o \
path_reparametrisation.o \
path_variables.o \
parallel_include.o \
parallel_types.o \
parameters.o \
parser.o \
paw_variables.o \
pseudo_types.o \
printout_base.o \
ptoolkit.o \
radial_grids.o \
random_numbers.o \
read_namelists.o \
read_ncpp.o \
read_upf_v1.o \
read_upf_v2.o \
read_uspp.o \
recvec.o \
shmem_include.o \
sic.o \
smallbox.o \
splinelib.o \
stick_base.o \
task_groups.o \
timestep.o \
twin_types.o \
version.o \
upf.o \
upf_to_internal.o \
uspp.o \
vxc_t.o \
vxcgc.o \
wave_base.o \
wavefunctions.o \
wannier.o \
write_upf_v2.o \
xml_input.o \
xml_io_base.o \
zhpev_drv.o \
wannier_new.o \
wrappers.o\
read_cards.o \
read_oddalpha_file.o \
compute_dipole.o #added by giovanni, from trunk version July 2012 

all : $(MODULES)

clean :
	- /bin/rm -f *.o *.d *.i *~ *.F90 *.mod *.L

version.o : version.h

version.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> version.h

include make.depend
