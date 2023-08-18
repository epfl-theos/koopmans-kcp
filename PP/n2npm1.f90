! DONE:
! 1) this is the desired comand with arguments from command line
!      <binary name> --fill  <spin_channel> <index> tmpdir_in tmpdir_out
!      <binary name> --empty <spin_channel ><index> tmpdir_in tmpdir_out
!    where <index> is the index of the orbital of the N-electron calculation 
!    we want to add/remove (NB: need a convention here on how to count empty states
!    and spin. I think it would be convenient to give the spin channel from input as
!    we have different evc files for different spin ) 
! 2) Add reading routine for empty state (only if --fill) 
! TODO: 
! 3) Check Makefile and the need of mpif90 wrapper (it loooks like we need it 
!    as iotk as some dependendence on mpi library. Also we might need a different 
!    rule to compile this program (no need for any precompilation flags like e.g. -D__MPI)
! RELEVANT observations:
! 1) for the N-1 it seems to me that no matter what the strcuture of the evc files for N-1
!    is identical to that of N. This is because the number of bands is determined by the 
!    majority spin channel (that does not change in structure going from N to N-1). 
!    It should be possible to restart a N-1 calculation from the N wfcs. Still it would 
!    be better to correctly initialize the last band in the minority wfc to ZERO (I guess 
!    this is what the code will do -it will read only N-1 out of N and set to zero the last
!    one- but better to be safe)
!
!--------------------------------------------------------------------------
PROGRAM n2npm1 
  !------------------------------------------------------------------------
  !
  IMPLICIT NONE
  CHARACTER(LEN=256)   :: task, dir_in, filename_in, dir_out, filename_out
  INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
  COMPLEX(DP), ALLOCATABLE :: evc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_empty(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_out(:,:,:)
  INTEGER :: ngw, nbnd(2), ispin, nspin, igwx, ibnd, ibnd_, spin_channel, ig
  INTEGER :: nbnd_emp(2)
  INTEGER :: nbnd_out(2)
  REAL(DP) :: nel(2)
  LOGICAL :: is_cmplx, l_fill, gamma_only, exst
  INTEGER :: index
  REAL(DP) :: scalef
  INTEGER :: iuni
  !
  !Defauls
  task="fill"             ! FIXME this will be dependent on the input  
  index=1                 ! FIXME this will be dependent on the input  
  dir_in='./'             ! FIXME this will be dependent on the input path
  dir_out='./'            ! FIXME this will be dependent on the input path
  spin_channel = 1        ! FIXME: this will track on which spin channel we need to add/remove the electron
  !                       ! Requires a way determine which spin channel from input (either explicite input or
  !                       ! inferred from index which should go from 1 to nbnd(1)+nbnd(2)
  !
  CALL read_command_line (task, index, spin_channel, dir_in, dir_out)
  !
  ! l_fill determine what to do: remove from an occupied one (l_fill=F)
  ! or add to empty (l_fill=T)
  l_fill=.false.
  IF (task=="fill") l_fill=.true.
  !
  IF (l_fill) THEN 
    WRITE(*, '(/, 3X, 2(A, I5))') "TASK: N --> N+1; Adding &
            one electron taken from empty orbital= ", index, " spin= ", spin_channel
  ELSE
    WRITE(*, '(/, 3X, 2(A, I5))') "TASK: N --> N-1; Removing &
            one electron from occupied orbital= ", index, " spin= ", spin_channel
  ENDIF
    WRITE(*,'(/, 5X, "SUMMARY: task     =", A)')  TRIM(task)
    WRITE(*,'(   5X, "         index    =", I5)') index
    WRITE(*,'(   5X, "         spin     =", I5)') spin_channel
    WRITE(*,'(   5X, "         dir_in   =", A)')  TRIM(dir_in)
    WRITE(*,'(   5X, "         dir_out  =", A)') TRIM(dir_out)
  !
  IF (dir_in == dir_out) THEN 
     WRITE(*,'(3X, 3A)') "ERROR: dir_out = dir_in"
     STOP
  ENDIF
  WRITE(*,'(/ 5X, A)') "Reading OCC manifold ..."
  ngw = 0; nbnd=0; nspin=0; nel=0.D0
  ! READ the shape of the wfc object (# of PW, # spin, complx/real, # bands) 
  filename_in=TRIM(dir_in)//'evc01.dat'
  CALL read_sizes (filename_in, ngw, nbnd(1), nspin, is_cmplx, igwx, scalef, gamma_only)
  !
  IF (nspin==2) THEN
     ! If needed READ the shape of the wfc object for spin down
     filename_in=TRIM(dir_in)//'evc02.dat'
     CALL read_sizes (filename_in, ngw, nbnd(2), nspin, is_cmplx, igwx, scalef, gamma_only)
  ENDIF
  !
  WRITE(*,'(/, 7X, "SUMMARY: igwx     =", I8)') igwx
  WRITE(*,'(   7X, "         ngw      =", I8)') ngw
  WRITE(*,'(   7X, "         nspin    =", I8)') nspin
  WRITE(*,'(   7X, "         nbnd(1)  =", I8)') nbnd(1)
  WRITE(*,'(   7X, "         nbnd(2)  =", I8)') nbnd(2)
  WRITE(*,'(   7X, "         is_cmplx =", L8)') is_cmplx
  WRITE(*,'(   7X, "         G_trick  =", L8, /)') gamma_only
  !
  ALLOCATE (evc(igwx, max(nbnd(1), nbnd(2)), nspin))
  !
  ! reset filename 
  filename_in=TRIM(dir_in)//'evc01.dat'
  DO ispin = 1, nspin
    IF (ispin==2 ) filename_in=TRIM(dir_in)//'evc02.dat'
    CALL read_wf_occ (filename_in, evc(:,:,ispin), ngw, nbnd(ispin)) 
    ! TODO: if KS state are used at initialization, also empty states are written to evc
    !       then this countig does not make sense
    CALL check_nele (ngw, nbnd(ispin), evc(:,:,ispin), nel(ispin))
    WRITE(*,'(7X, "CHECK: ispin =", I5, " nel =", F18.12, " INT(nel)=", I5)') &
            ispin, nel(ispin), NINT(nel(ispin))
    !WRITE(*,) ispin, nel(ispin), NINT(nel(ispin))
  ENDDO
  !
  ! At this stage I read the empty manifold if needed 
  IF (l_fill) THEN
     WRITE(*,'(/ 5X, A)') "Reading EMP manifold ..."
     ! Here we need to read also empty states if needed.
     ! First the shape of the arrey ngw, nbnd_emp for spin up
     iuni=789
     filename_in=TRIM(dir_in)//'evc0_empty1.dat'
     INQUIRE (FILE=TRIM(filename_in), EXIST=EXST)
     IF (.NOT. exst) THEN 
        WRITE(*,'(3A)') "File ", TRIM(filename_in), " do not exist"
        STOP
     ENDIF
     OPEN (UNIT=iuni, FILE=TRIM(filename_in), status='unknown', FORM='UNFORMATTED')
     READ (iuni) ngw, nbnd_emp(1)
     ! then for spin down channel
     IF (nspin == 2) THEN 
       iuni=iuni+1
       filename_in=TRIM(dir_in)//'evc0_empty2.dat'
       INQUIRE (FILE=TRIM(filename_in), EXIST=EXST)
       IF (.NOT. exst) THEN
          WRITE(*,'(3A)') "File ", TRIM(filename_in), " do not exist"
          STOP
       ENDIF
       OPEN (UNIT=iuni, FILE=TRIM(filename_in), status='unknown', FORM='UNFORMATTED')
       READ (iuni) ngw, nbnd_emp(2)
     ENDIF
     WRITE(*,'(/, 7X, "SUMMARY: ngw      =", I8)') ngw
     WRITE(*,'(   7X, "         nbnd(1)  =", I8)') nbnd_emp(1)
     WRITE(*,'(   7X, "         nbnd(2)  =", I8)') nbnd_emp(2)
     ! allocation
     ALLOCATE ( evc_empty (ngw, MAX(nbnd_emp(1), nbnd_emp(2)),nspin) )
     !
     ! read the empty-state wfcs
     iuni=789
     filename_in=TRIM(dir_in)//'evc0_empty1.dat'
     DO ispin = 1, nspin
       IF (ispin == 2) filename_in=TRIM(dir_in)//'evc0_empty2.dat'
       IF (ispin == 2) iuni=iuni+1
       DO ibnd = 1, nbnd_emp(ispin)
          READ (iuni) (evc_empty(ig,ibnd,ispin), ig=1, ngw)
       ENDDO
     ENDDO
  ENDIF
  !
  ! check if the index of the orbital is within a valid range
  IF (l_fill) THEN
    !
    IF (spin_channel == 1 .AND. index .gt. nbnd_emp(1)) THEN
       WRITE(*,*) "index out of range ", index, nbnd_emp(1);  STOP
    ENDIF
    IF (spin_channel == 2 .AND. index .gt. nbnd_emp(2)) THEN
       WRITE(*,*) "index out of range ", index, nbnd_emp(2);  STOP
    ENDIF
    !
  ELSE
    !
    IF (spin_channel == 1 .AND. index .gt. nbnd(1)) THEN
       WRITE(*,*) "index out of range ", index, nbnd(1);  STOP
    ENDIF
    IF (spin_channel == 2 .AND. index .gt. nbnd(2)) THEN
       WRITE(*,*) "index out of range ", index, nbnd(2);  STOP
    ENDIF
    !
  ENDIF
  !
  ! From now on we reshape the wfc object
  WRITE(*,'(/ 5X, A)') "Re-shaping OCC manifold ..."
  nbnd_out = nbnd
  IF (l_fill) THEN 
     nbnd_out(spin_channel) = nbnd(spin_channel) + 1
     nel(spin_channel) = nel(spin_channel)+1
     IF (nel(1) .lt. nel(2)) THEN
       WRITE(*,*) "nel up must be greater than or equal to nel dw ", NINT(nel(:));  STOP
     ENDIF
     ALLOCATE (evc_out(igwx, max(nbnd_out(1), nbnd_out(2)), nspin))
     evc_out = CMPLX(0.D0, 0.D0, kind =DP)
     DO ispin = 1, nspin
       DO ibnd = 1,nbnd(ispin)
         evc_out(:,ibnd,ispin) = evc(:,ibnd,ispin)
       ENDDO
     ENDDO
     evc_out(:,nbnd_out(spin_channel),spin_channel) = evc_empty(:,index,spin_channel) ! 
     !
  ELSE
     nbnd_out(spin_channel) = nbnd(spin_channel) - 1
     nel(spin_channel) = nel(spin_channel)-1
     IF (nel(1) .lt. nel(2)) THEN
       WRITE(*,*) "nel up must be greater than or equal to nel dw ", NINT(nel(:));  STOP
    ENDIF
     ALLOCATE (evc_out(igwx, max(nbnd_out(1), nbnd_out(2)), nspin))
     evc_out = CMPLX(0.D0, 0.D0, kind =DP)
     DO ispin = 1, nspin
       ibnd_=0
       DO ibnd = 1,nbnd(ispin)
         IF (ispin == spin_channel .AND. ibnd == index) CYCLE
         ibnd_=ibnd_+1
         evc_out(:,ibnd_,ispin) = evc(:,ibnd,ispin)
       ENDDO
       WRITE(*,*) ispin,  evc_out(1,1,ispin), evc_out(igwx,MAX(nbnd_out(1),nbnd_out(2)),ispin)! check debug
     ENDDO
  ENDIF
  !
  WRITE(*,'(/, 7X, "SUMMARY: igwx     =", I8)') igwx
  WRITE(*,'(   7X, "         ngw      =", I8)') ngw
  WRITE(*,'(   7X, "         nspin    =", I8)') nspin
  WRITE(*,'(   7X, "         nbnd(1)  =", I8)') MAX(nbnd_out(1),nbnd_out(2))
  WRITE(*,'(   7X, "         nbnd(2)  =", I8)') MAX(nbnd_out(1),nbnd_out(2))
  WRITE(*,'(   7X, "         is_cmplx =", L8)') is_cmplx
  WRITE(*,'(   7X, "         G_trick  =", L8, /)') gamma_only
  !
  WRITE(*,'(/ 5X, A)') "Writing new WFCs to disk ..."
  !
  filename_out=TRIM(dir_out)//'evc01.dat'
  DO ispin = 1, nspin
    IF (ispin==2 ) filename_out=TRIM(dir_out)//'evc02.dat'
    CALL write_wf (filename_out, evc_out(:,:,ispin), nspin, ispin, ngw, igwx, &
            MAX(nbnd_out(1), nbnd_out(2)), is_cmplx, scalef, gamma_only)
  ENDDO
  !
  WRITE(*,'(/ 3X, A, /)') "GAME OVER" 
  CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE read_command_line (task, index, spin_channel, dir_in, dir_out) 
     !------------------------------------------------------------------
     !     arg#      1      2      3       4        5
     ! n2npm1.x --<task> <index> <spin> <dir_in> <dir_out>
     !
     IMPLICIT NONE 
     INTEGER :: index, spin_channel
     CHARACTER(LEN=256) :: task, dir_in, dir_out
     INTEGER nargs, iarg, counter
     CHARACTER(LEN=256) :: arg
     !
     nargs = command_argument_count()
     !WRITE(*,*) nargs
     !   
     IF ( nargs /=  5 ) THEN
        WRITE(*,'(/, 5X, A, I5)') "Wrong number of arguments: STOP", nargs
        WRITE(*,'(   5X, A, /)') "Usage: n2npm1.x --<task> <index> <spin_channel> <dir_in> <dir_out>"
        STOP
     ENDIF
     !
     iarg = 1
     counter = 0
     !
     CALL get_command_argument( iarg, arg )
     iarg = iarg + 1
     !
     SELECT CASE ( trim(arg) )
     CASE ( '--fill' )
        task="fill"
     CASE ( '--empty')
        task="empty"
     CASE DEFAULT
       WRITE(*,'(A)') 'unrecognised argument option'
       STOP
     END SELECT
     !
     CALL get_command_argument( iarg, arg )
     READ(arg,*) index
     iarg = iarg + 1
     CALL get_command_argument( iarg, arg )
     READ(arg,*) spin_channel
     iarg = iarg + 1
     CALL get_command_argument( iarg, dir_in )
     iarg = iarg + 1
     CALL get_command_argument( iarg, dir_out )
     !
     RETURN
     !
  END SUBROUTINE
  !
  !----------------------------------------------------------------------
  SUBROUTINE write_wf(filename, wf, nspin, ispin, ngw, igwx, nbnd, is_cmplx, scalef, gamma_only)
    !--------------------------------------------------------------------
    !
    USE iotk_module
    !
    IMPLICIT NONE
    INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
    CHARACTER(LEN=256), INTENT(IN)   :: filename
    COMPLEX(DP) :: wf(ngw, nbnd) 
    INTEGER, INTENT(IN) :: ngw, nbnd, nspin, igwx, ispin
    LOGICAL, INTENT(IN) :: is_cmplx, gamma_only
    REAL(DP), INTENT(IN) :: scalef
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: j
    COMPLEX(DP), ALLOCATABLE :: wtmp(:)
    INTEGER                  :: iuni
    !
    iuni=789
    !
    CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
    !
    CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
    CALL iotk_write_attr( attr, "igwx",         igwx )
    CALL iotk_write_attr( attr, "do_wf_cmplx",   is_cmplx ) !added:giovanni
    CALL iotk_write_attr( attr, "gamma_only",   gamma_only.and..not.is_cmplx )
    CALL iotk_write_attr( attr, "nbnd",         nbnd )
    CALL iotk_write_attr( attr, "ik",           ispin )
    CALL iotk_write_attr( attr, "nk",           nspin )
    CALL iotk_write_attr( attr, "ispin",        ispin )
    CALL iotk_write_attr( attr, "nspin",        nspin )
    CALL iotk_write_attr( attr, "scale_factor", scalef )
    !
    CALL iotk_write_empty( iuni, "INFO", attr )
    !
    ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
    wtmp = 0.0_DP
    DO j = 1, nbnd
       wtmp(:)=wf(:,j)
      CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
    ENDDO
    !
    CALL iotk_close_write( iuni )
    !
  END SUBROUTINE
  !
  !----------------------------------------------------------------------
  SUBROUTINE check_nele (ngw, nbnd, evc, nele)
    !--------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
    INTEGER, INTENT(IN) :: ngw, nbnd
    COMPLEX(DP), INTENT(IN) :: evc(ngw,nbnd)
    REAL(DP), INTENT(OUT) :: nele
    INTEGER :: ibnd
    !
    nele=0.D0
    DO ibnd = 1, nbnd
      nele = nele + SUM(evc(:,ibnd)*CONJG(evc(:,ibnd)))
      !WRITE(*,'(I5,3x , 2F20.15, 3x, F10.6)') ibnd, SUM(evc(:,ibnd)*CONJG(evc(:,ibnd))), nele
    ENDDO
    RETURN
    !
  END SUBROUTINE 
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_sizes(filename, ngw, nbnd, nspin, is_cmplx, igwx, scalef, gamma_only)
    !------------------------------------------------------------------------
    !
    USE iotk_module
    !
    IMPLICIT NONE
    CHARACTER(LEN=256), INTENT(IN)   :: filename
    INTEGER, INTENT(OUT) :: ngw, nbnd,  nspin, igwx
    LOGICAL, INTENT(OUT) :: is_cmplx, gamma_only
  
    INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: j, ispin
    INTEGER                  :: ierr
    INTEGER                  :: ik, nk
    INTEGER                  :: iuni
    REAL(DP)                 :: scalef 
    !
    ierr = 0
    iuni=789
    !
    CALL iotk_open_read( iuni, FILE = filename, BINARY = .TRUE., IERR = ierr )
    !
    IF (ierr /= 0) THEN 
       WRITE(*,*) "Something wrong while reading file = ", filename
       STOP
    ENDIF
    !
    CALL iotk_scan_empty( iuni, "INFO", attr )
    !
    !WRITE(*,*) "NICOLA", ngw, nbnd, ik, nk, nspin, igwx, is_cmplx
    CALL iotk_scan_attr( attr, "ngw",          ngw )
    CALL iotk_scan_attr( attr, "nbnd",         nbnd )
    CALL iotk_scan_attr( attr, "ik",           ik )
    CALL iotk_scan_attr( attr, "nk",           nk )
    CALL iotk_scan_attr( attr, "nspin",        nspin )
    CALL iotk_scan_attr( attr, "ispin",        ispin )
    CALL iotk_scan_attr( attr, "igwx",         igwx )
    CALL iotk_scan_attr( attr, "do_wf_cmplx",  is_cmplx )
    CALL iotk_scan_attr( attr, "gamma_only",  gamma_only )
    CALL iotk_scan_attr( attr, "scale_factor", scalef )
    !
    !WRITE(*,*) "NICOLA", ispin, ngw, nbnd, ik, nk, nspin, igwx, is_cmplx, scalef
    CALL iotk_close_read( iuni )
    RETURN
    !
  END subroutine
  !
  !------------------------------------------------------------------------
  SUBROUTINE read_wf_occ(filename, wf, ngw, nbnd)
    !------------------------------------------------------------------------
    !
    USE iotk_module
    !
    IMPLICIT NONE
    INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
    CHARACTER(LEN=256), INTENT(IN)   :: filename
    COMPLEX(DP) :: wf(ngw, nbnd) 
    INTEGER, INTENT(IN) :: ngw, nbnd
    CHARACTER(iotk_attlenx)  :: attr
    INTEGER                  :: j
    COMPLEX(DP), ALLOCATABLE :: wtmp(:)
    INTEGER                  :: ierr
    INTEGER                  :: igwx, ngw_, nbnd_
    INTEGER                  :: iuni
    LOGICAL                  :: is_cmplx
    !
    ierr = 0
    iuni=789
    !
    CALL iotk_open_read( iuni, FILE = filename, &
                            BINARY = .TRUE., IERR = ierr )
    !
    IF (ierr /= 0) THEN 
       WRITE(*,*) "Something wrong while reading file = ", filename
       STOP
    ENDIF
    !
    CALL iotk_scan_empty( iuni, "INFO", attr )
    !
    ngw_=0; nbnd_=0; 
    !WRITE(*,*) "NICOLA", ngw_, nbnd_
    CALL iotk_scan_attr( attr, "ngw",          ngw_ )
    CALL iotk_scan_attr( attr, "nbnd",         nbnd_ )
    CALL iotk_scan_attr( attr, "igwx",         igwx )
    CALL iotk_scan_attr( attr, "ispin",        ispin )
    CALL iotk_scan_attr( attr, "do_wf_cmplx",  is_cmplx )
    !
    !WRITE(*,*) "NICOLA", TRIM(filename), ngw_, nbnd_, ispin, igwx, is_cmplx
    !
    IF (nbnd_ /= nbnd .OR. ngw_ /= ngw) THEN 
       WRITE(*,*) "Size Mismatch. STOP"
       STOP
    ENDIF
    
    ALLOCATE( wtmp( igwx ) )
    !
    DO j = 1, nbnd
       CALL iotk_scan_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
       wf(:,j) = wtmp (:) 
    ENDDO
    !
    DEALLOCATE (wtmp)
    !
    CALL iotk_close_read( iuni )
    RETURN
    !
  END subroutine
  !
END PROGRAM
