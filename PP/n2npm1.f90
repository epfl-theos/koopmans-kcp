! TO DO:
! 1) this is the desired comand with arguments from command line
!      <binary name> --fill  <spin_channel> <index> tmpdir_in tmpdir_out
!      <binary name> --empty <spin_channel ><index> tmpdir_in tmpdir_out
!    where <index> is the index of the orbital of the N-electron calculation 
!    we want to add/remove (NB: need a convention here on how to count empty states
!    and spin. I think it would be convenient to give the spin channel from input as
!    we have different evc files for different spin ) 
! 2) Add reading routine for empty state (only if --fill) 
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
  INTEGER :: ngw, nbnd(2), ispin, nspin, igwx, ibnd, ibnd_, spin_channel
  INTEGER :: nbnd_out(2)
  REAL(DP) :: nel(2)
  LOGICAL :: is_cmplx, l_fill, gamma_only
  INTEGER :: index
  REAL(DP) :: scalef
  !
  task="empty"            ! FIXME this will be dependent on the input  
  index=1                 ! FIXME this will be dependent on the input  
  dir_in='./'             ! FIXME this will be dependent on the input path
  dir_out='./'            ! FIXME this will be dependent on the input path
  spin_channel = 1        ! FIXME: this will track on which spin channel we need to add/remove the electron
  !                       ! Requires a way determine which spin channel from input (either explicite input or
  !                       ! inferred from index which should go from 1 to nbnd(1)+nbnd(2)
  !
  !CALL get_command_argument( 1, arg )
  !CALL parse_args( io_files, nrtot, output_exst )
  !
  ! l_fill determine what to do: remove from an occupied one (l_fill=F)
  ! or add to empty (l_fill=T)
  l_fill=.false.
  IF (task=="fill") l_fill=.true.
  !
  IF (l_fill) THEN 
    WRITE(*, '(/, 3X, 2(A, I5))') "TASK: Going to add &
            one electron taken from empty orbital= ", index, " spin= ", spin_channel
  ELSE
    WRITE(*, '(/, 3X, 2(A, I5))') "TASK: Going to remove &
            one electron from occupied orbital= ", index, " spin= ", spin_channel
  ENDIF
  !
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
  WRITE(*,*) igwx, ngw, nbnd(:), nspin
  !
  ALLOCATE (evc(igwx, max(nbnd(1), nbnd(2)), nspin))
  !
  ! reset filename 
  filename_in=TRIM(dir_in)//'evc01.dat'
  DO ispin = 1, nspin
    IF (ispin==2 ) filename_in=TRIM(dir_in)//'evc02.dat'
    CALL read_wf_occ (filename_in, evc(:,:,ispin), ngw, nbnd(ispin)) 
    CALL check_nele (ngw, nbnd(ispin), evc(:,:,ispin), nel(ispin))
    WRITE(*,*) ispin, nel(ispin), NINT(nel(ispin)), evc(1,1,ispin), evc(igwx,nbnd(ispin),ispin)! check debug
  ENDDO
  !
  ! At this stage I read the empty manifold if needed 
  IF (l_fill) THEN
     ! Here we need to read also empty states if 
     WRITE(*,*) "Case l_fill = T NOT IMPLEMENTED YET: STOP"
     STOP
  ENDIF
  !
  ! check if the index of the orbital is within a valid range
  IF (l_fill) THEN
    !
  ELSE
    IF (index .gt. nbnd(1)+nbnd(2)) THEN
       WRITE(*,*) "index out of range ", index, nbnd(1)+nbnd(2)
       STOP
    ENDIF
  ENDIF
  !
  ! From now on we reshape the wfc object
  nbnd_out = nbnd
  IF (l_fill) THEN 
     nbnd_out(spin_channel) = nbnd(spin_channel) + 1
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
  filename_out=TRIM(dir_out)//'evc01_.dat'
  DO ispin = 1, nspin
    IF (ispin==2 ) filename_out=TRIM(dir_out)//'evc02_.dat'
    CALL write_wf (filename_out, evc_out(:,:,ispin), nspin, ispin, ngw, igwx, &
            MAX(nbnd(1), nbnd(2)), is_cmplx, scalef, gamma_only)
  ENDDO
  CONTAINS
  !
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
    WRITE(*,*) "NICOLA", ispin, ngw, nbnd, ik, nk, nspin, igwx, is_cmplx, scalef
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
    WRITE(*,*) "NICOLA", ngw_, nbnd_
    CALL iotk_scan_attr( attr, "ngw",          ngw_ )
    CALL iotk_scan_attr( attr, "nbnd",         nbnd_ )
    CALL iotk_scan_attr( attr, "igwx",         igwx )
    CALL iotk_scan_attr( attr, "ispin",        ispin )
    CALL iotk_scan_attr( attr, "do_wf_cmplx",  is_cmplx )
    !
    WRITE(*,*) "NICOLA", TRIM(filename), ngw_, nbnd_, ispin, igwx, is_cmplx
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
