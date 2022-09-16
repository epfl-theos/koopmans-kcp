set(src_clib
  stack.c
  c_mkdir.c
  cptimer.c
  eval_infix.c
  fft_stick.c
  indici.c
  memstat.c
  qsort.c)

qe_add_library(qe_clib ${src_clib})

qe_install_targets(qe_clib)
