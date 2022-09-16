set(src_flib
  avrec.f90
  atomic_number.f90
  bachel.f90
  capital.f90
  dost.f90
  erf.f90
  functionals.f90
  lsda_functionals.f90
  more_functionals.f90
  iglocal.f90
  inpfile.f90
  int_to_char.f90
  invmat.f90
  invmat_complex.f90
  latgen.f90
  linpack.f90
  matches.f90
  recips.f90
  remove_tot_torque.f90
  simpsn.f90
  sort.f90
  sph_bes.f90
  sph_dbes.f90
  transto.f90
  date_and_tim.f90
  sort_gvec.f90
  volume.f90
  dylmr2.f90
  ylmr2.f90
  cryst_to_car.f90)

qe_add_library(qe_flib ${src_flib})

target_link_libraries(qe_flib
  PRIVATE
    qe_modules)

qe_install_targets(qe_flib)
