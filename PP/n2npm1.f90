! TO DO:
! 1) this is the desired comand with arguments from command line
!    <binary name> --fill  <index> tmpdir_in tmpdir_out
!    <binary name> --empty <index> tmpdir_in tmpdir_out
! 2) Add reading routine for empty state (only if --empty) 
! 3) Check Makefile and the need of mpif90 wrapper (it loooks like we need it 
!    as iotk as some dependendence on mpi library. Also we might need a different 
!    rule to compile this program (no need for any precompilation flags like e.g. -D__MPI)
PROGRAM n2npm1 
  !------------------------------------------------------------------------
  !
  IMPLICIT NONE
  CHARACTER(LEN=256)   :: filename
  INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
  COMPLEX(DP), ALLOCATABLE :: evc(:,:,:)
  INTEGER :: ngw, nbnd(2), ispin, nspin
  LOGICAL :: is_cmplx
  !
  ngw = 0; nbnd=0; nspin=0
  filename='evc01.dat' ! FIXME this will be dependent on the input path
  !
  CALL read_sizes (filename, ngw, nbnd(1), nspin, is_cmplx)
  !
  IF (nspin==2) THEN
     filename='evc02.dat'
     CALL read_sizes (filename, ngw, nbnd(2), nspin, is_cmplx)
  ENDIF
  !
  WRITE(*,*) ngw, nbnd(:), nspin
  ALLOCATE (evc(ngw, max(nbnd(1), nbnd(2)), nspin))
  !
  ! reset filename
  filename='evc01.dat'
  DO ispin = 1, nspin
    IF (ispin==2 ) filename="evc02.dat"
    CALL read_wf_occ (filename, evc(:,:,ispin), ngw, nbnd(ispin)) 
    WRITE(*,*) ispin, evc(1:3,1,ispin) ! check debug
  ENDDO
  !
  CONTAINS
  !
  !
  SUBROUTINE read_sizes(filename, ngw, nbnd, nspin, is_cmplx)
  !------------------------------------------------------------------------
  !
  USE iotk_module
  !
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN)   :: filename
  INTEGER, INTENT(OUT) :: ngw, nbnd,  nspin
  LOGICAL, INTENT(OUT) :: is_cmplx

  INTEGER, PARAMETER   :: DP = selected_real_kind(14,200)
  CHARACTER(iotk_attlenx)  :: attr
  INTEGER                  :: j
  INTEGER                  :: ierr
  INTEGER                  :: igwx, ik, nk
  INTEGER                  :: iuni
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
  WRITE(*,*) "NICOLA", ngw, nbnd, ik, nk, nspin, igwx, is_cmplx
  CALL iotk_scan_attr( attr, "ngw",          ngw )
  CALL iotk_scan_attr( attr, "nbnd",         nbnd )
  CALL iotk_scan_attr( attr, "ik",           ik )
  CALL iotk_scan_attr( attr, "nk",           nk )
  CALL iotk_scan_attr( attr, "nspin",        nspin )
  CALL iotk_scan_attr( attr, "igwx",         igwx )
  CALL iotk_scan_attr( attr, "do_wf_cmplx",  is_cmplx )
  !
  WRITE(*,*) "NICOLA", ngw, nbnd, ik, nk, nspin, igwx, is_cmplx
  CALL iotk_close_read( iuni )
  RETURN
  !
  END subroutine
  !
  !
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
  WRITE(*,*) "NICOLA", ngw_, nbnd_, ispin, igwx, is_cmplx
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
