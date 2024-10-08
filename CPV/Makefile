# Makefile for CP/FPMD

include	../make.sys

FOBJS = \
atoms_type.o \
bessel.o \
berryion.o \
bforceion.o \
centers_and_spreads.o \
cg.o \
cg_sub.o \
cglib.o \
chargedensity.o \
chargemix.o \
chi2.o \
compute_fes_grads.o \
compute_scf.o \
cp_autopilot.o \
cp_emass.o \
cp_fpmd.o \
cp_interfaces.o \
cp_restart.o \
cplib.o \
cpr_mod.o \
cpr.o \
cprsub.o \
dealloc.o \
dforceb.o \
eelib.o \
efermi.o \
efield.o \
eigs0.o \
electrons.o \
emptystates.o \
ensemble_dft.o \
environment.o \
exch_corr.o \
exx_divergence.o \
fft.o \
forces.o \
fromscra.o \
gradrho.o \
gtable.o \
init.o \
init_run.o \
inner_loop.o \
inner_loop_cold.o \
inner_loop_smear.o \
input.o \
io_pot_sic_xml.o \
ions_positions.o \
ksstates.o \
main.o \
mainvar.o \
main_loops.o \
cplib_meta.o \
metadyn_base.o \
metadyn_io.o \
metadyn_vars.o \
metaxc.o \
modules.o \
move_electrons.o \
nksiclib.o \
nksic_corrections.o \
hflib.o \
nl_base.o \
nlcc.o \
ortho_base.o \
ortho.o \
ortho_check.o \
para.o \
path_base.o \
path_routines.o \
phasefactor.o \
polarization.o \
potentials.o \
pres_ai_mod.o \
print_out.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
pseudopot_sub.o \
qmatrixd.o \
qqberry.o \
read_pseudo.o \
restart.o \
restart_sub.o \
runcp.o \
spharmonic.o \
spline.o \
stress.o \
stop_run.o \
symm_wannier.o \
turbo.o \
util.o \
vanderwaals.o \
vol_clu.o \
wannier_base.o \
wannier.o \
wave.o \
wave_types.o \
wf.o \
write_hamiltonian.o \
writetofile.o \
gram_swap.o \
odd_alpha.o \
cg_empty_sub.o \
wave_init_wannier.o \
inner_loop_generalize.o \
makov_payne.o \
pc3nc_fixed.o \
empty_koopmans_pp.o \
perturbing_pot.o \
print_real_space_orbital_density.o

# dgradcorr.o 
# setup_gga.o

LOBJS = \
adjef.o \
entropy.o \
forceconv.o \
geninv.o \
indices.o

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

AFCLIB=../AFC90/src/libafc90.a

TLDEPS= bindir mods libs libiotk afclib

all : kcp
kcp : tldeps libcp.a kcp.x cppp.x

kcp.x : cprstart.o libcp.a $(LIBOBJS) $(AFCLIB)
	$(LD) $(LDFLAGS) -o kcp.x \
		cprstart.o $(MODULES) libcp.a \
                $(LIBOBJS) $(AFCLIB) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/kcp.x . )

cp_test : tldeps libcp.a cp_test.x cppp.x

cp_test.x : cprstart.o libcp.a $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cp_test.x \
		cprstart.o $(MODULES) libcp.a $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp_test.x . )

cp_bad_sd : tldeps libcp.a cp_bad_sd.x cppp.x

cp_bad_sd.x : cprstart.o libcp.a $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cp_bad_sd.x \
		cprstart.o $(MODULES) libcp.a $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp_bad_sd.x . )

libcp.a : $(FOBJS) $(WAN90) $(LOBJS) 
	 $(AR) $(ARFLAGS) $@ $?
	 $(RANLIB) $@

PPOBJS = \
../Modules/xml_io_base.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/mp_global.o \
../Modules/io_global.o \
../Modules/io_files.o \
../Modules/mp_wave.o \
../Modules/parser.o \
../Modules/kind.o \
../Modules/control_flags.o \
../Modules/parameters.o \
../Modules/parallel_include.o \
../Modules/error_handler.o \
../Modules/constants.o \
../Modules/wrappers.o

cppp.x : fpmdpp.o $(PPOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cppp.x fpmdpp.o $(PPOBJS) $(LIBOBJS) $(LIBS)
	- (cd ../bin ; ln -fs ../CPV/cppp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f cppp.x *.o *.mod cpver.h *.i core* *.F90 fort* \
		*.cpp *.d *.L *.a *.s kcp.x

include make.depend
