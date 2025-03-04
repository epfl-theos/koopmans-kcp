
#!/bin/sh
# check whether ECHO has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

. ../ENVIRONMENT 

$ECHO "  checking that needed directories and files exist...\c"

# check for directories
for DIR in "$BIN_DIR" "$PSEUDO_DIR" ; do
    if test ! -d $DIR ; then
        $ECHO
        $ECHO "ERROR: $DIR not existent or not a directory"
	$ECHO "       check the ENVIRONMENT file"
        $ECHO "       Aborting"
        exit 1
    fi
done

BIN_LIST="kcp.x"

# check for executables
for FILE in $BIN_LIST ; do
    if test ! -x $BIN_DIR/$FILE ; then
        $ECHO
        $ECHO "ERROR: $BIN_DIR/$FILE not existent or not executable"
	$ECHO "       check the ENVIRONMENT file"
        $ECHO "       Aborting"
        exit 1
    fi
done
$ECHO " done"


echo ""
echo "Running kcp.x as" $PARA_PREFIX $BIN_DIR/kcp.x
echo ""
cat > dft.cpi <<EOF
&CONTROL
   calculation      = "cp"
   verbosity        = "low"
   restart_mode     = "from_scratch"
   iprint           = 1
   outdir           = "./TMP-CP/"
   prefix           = "kc"
   disk_io          = "high"
   pseudo_dir       = "$PSEUDO_DIR"
   ndr              = 98
   ndw              = 98
   write_hr         = .false.
   print_real_space_density = .false.
/
&SYSTEM
   ibrav            = 0
   nbnd             = 4
   tot_charge       = 0
   ecutwfc          = 35.0
   nspin            = 1
   do_orbdep        = .false.
   fixed_state      = .false.
   do_ee            = .true.
   nelec            = 8
   do_wf_cmplx      = .false.
   ntyp             = 1
   nat              = 1
/
&ELECTRONS
   conv_thr         = 1.8000000000000002e-08
   ortho_para       = 1
   maxiter          = 300
   electron_dynamics = "cg"
   passop           = 2.0
   do_outerloop     = .true.
   do_outerloop_empty = .true.
/
&IONS
   ion_dynamics     = "none"
   ion_damping      = 0.5 
/
&CELL
/
&EE
   which_compensation = "tcc"
/
&NKSIC
   do_innerloop     = .false.
   innerloop_cg_nreset = 20
   innerloop_cg_nsd = 2
   innerloop_init_n = 3
   hartree_only_sic = .false.
   esic_conv_thr    = 1.8000000000000002e-08
   do_innerloop_cg  = .true.
   innerloop_nmax   = 100
   do_innerloop_empty = .false.
/

ATOMIC_SPECIES
Ar 1.000 Ar_ONCV_PBE-1.2.upf

CELL_PARAMETERS angstrom
8.0000000000000 0.00000000000000 0.00000000000000
0.00000000000000 8.00000000000000 0.00000000000000
0.00000000000000 0.00000000000000 8.00000000000000

ATOMIC_POSITIONS angstrom
Ar 4.000 4.000 4.000
EOF

$PARA_PREFIX $BIN_DIR/kcp.x -in dft.cpi > dft.cpo 


cat > kipz.cpi <<EOF
&CONTROL
   calculation      = "cp"
   verbosity        = "low"
   restart_mode     = "restart"
   iprint           = 1
   outdir           = "./TMP-CP/"
   prefix           = "kc"
   disk_io          = "high"
   pseudo_dir       = "$PSEUDO_DIR"
   ndr              = 98
   ndw              = 99
   write_hr         = .false.
   print_real_space_density = .false.
   nstep = 100
   dt    = 50 
/
&SYSTEM
   ibrav            = 0
   nbnd             = 4
   tot_charge       = 0
   ecutwfc          = 35.0
   nspin            = 1
   do_orbdep        = .true.
   fixed_state      = .false.
   do_ee            = .true.
   nelec            = 8
   do_wf_cmplx      = .false.
   ntyp             = 1
   nat              = 1
/
&ELECTRONS
   conv_thr         = 1.8000000000000002e-08
   ortho_para       = 1
   maxiter          = 10
   electron_dynamics = "cg"
   passop           = 2.0
   do_outerloop     = .true.
   do_outerloop_empty = .true.
/
&IONS
   ion_dynamics     = "none"
   ion_damping = 0.5
/
&CELL
/
&EE
   which_compensation = "tcc"
/
&NKSIC
   which_orbdep     = "nkipz"
   nkscalfact       = 1.0 
   do_innerloop     = .false.
   innerloop_cg_nreset = 20
   innerloop_cg_nsd = 2
   innerloop_init_n = 3
   hartree_only_sic = .false.
   esic_conv_thr    = 1.8000000000000002e-08
   do_innerloop_cg  = .true.
   innerloop_nmax   = 100
   do_innerloop_empty = .false.
   l_group_minimization = .true.
   group_dimensions(1) = 1
   group_dimensions(2) = 3
/

ATOMIC_SPECIES
Ar 1.000 Ar_ONCV_PBE-1.2.upf

CELL_PARAMETERS angstrom
8.0000000000000 0.00000000000000 0.00000000000000
0.00000000000000 8.00000000000000 0.00000000000000
0.00000000000000 0.00000000000000 8.00000000000000

ATOMIC_POSITIONS angstrom
Ar 4.000 4.000 4.000
EOF

$PARA_PREFIX $BIN_DIR/kcp.x -in kipz.cpi > kipz.cpo 

