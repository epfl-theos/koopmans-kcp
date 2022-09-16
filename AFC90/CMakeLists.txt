set(src_afc90
  src/afcc.f90
  src/besseli.f90
  src/besselk.f90
  src/cylharm.f90
  src/cylharm0.f90
  src/cylharmasympt.f90
  src/cylharmlogslope.f90
  src/cylharmrk.f90
  src/cylharmseries.f90
  src/cylharmslope.f90
  src/ei.f90
  src/eimlog.f90
  src/fft1d.f90
  src/fft2d.f90
  src/fft3d.f90
  src/gaussian.f90
  src/gaussiank.f90
  src/gaussianl.f90
  src/interp.f90
  src/main.f90
  src/nfft.f90
  src/phi0d.f90
  src/phi1d.f90
  src/phi2d.f90
  src/phi3d.f90
  src/reciprocal.f90
  src/rtoaxis.f90
  src/rungekutta.f90
  src/steprk.f90
  src/vectorproduct.f90
  src/volume.f90
  src/ydot.f90)

qe_add_library(qe_afc90 ${src_afc90})

qe_install_targets(qe_afc90)
