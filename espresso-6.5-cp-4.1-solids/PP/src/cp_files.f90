!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!---------------------------------------------------------------------
MODULE cp_files
  !-------------------------------------------------------------------
  !
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: write_wannier_cp
  !
  CONTAINS
  !
  !-------------------------------------------------------------------
  SUBROUTINE write_wannier_cp( evc, typ )
    !-----------------------------------------------------------------
    !
    ! ...  This routine takes the Wannier functions in input and 
    ! ...  writes them into a file readable by the CP-Koopmans code
    ! 
    ! ...  typ = 'occupied' ---> output file: 'evc_occupied.dat'
    ! ...  typ = 'empty'    ---> output file: 'evc0_empty.dat'
    !
    USE kinds,               ONLY : DP
    USE io_global,           ONLY : ionode, ionode_id
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE mp_wave,             ONLY : mergewf
    USE mp_world,            ONLY : mpime, nproc
    USE mp,                  ONLY : mp_sum
    USE noncollin_module,    ONLY : npol
    USE fft_supercell,       ONLY : npwxcp, ig_l2g_cp
    USE read_wannier,        ONLY : num_bands, num_kpts 
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: evc(:,:,:)        ! input wfc
    CHARACTER(LEN=*), INTENT(IN) :: typ         ! state type
    !
    INTEGER :: io_level = 1
    INTEGER :: cp_unit = 124
    INTEGER :: npw_g           ! global number of PWs
    INTEGER :: ir, ibnd, ibnd_, ipw
    COMPLEX(DP), ALLOCATABLE :: evc_g(:)
    CHARACTER(LEN=33) :: filename
    !
    !
    ! ... defining output file name for occ/emp states
    !
    IF ( trim(typ) .eq. 'occupied' ) THEN
      filename = 'evc_occupied.dat'
    ELSE IF ( trim(typ) .eq. 'empty' ) THEN
      filename = 'evc0_empty.dat'
    ELSE
      CALL errore( 'write_wannier_cp', 'Invalid value for input variable typ', 1 )
    ENDIF
    !
    !
    npw_g = npwxcp
    CALL mp_sum( npw_g, intra_bgrp_comm )
    ALLOCATE( evc_g(npw_g) )
    !
    IF ( ionode ) THEN
      OPEN( UNIT=cp_unit, FILE=trim(filename), STATUS='unknown', FORM='unformatted' )
      WRITE( cp_unit ) npw_g, num_bands*2
    ENDIF
    !
    ! ... here we gather the wfc from all the processes
    ! ... and we write it to file
    !
    DO ir = 1, num_kpts
      DO ibnd = 1, num_bands*2
        !
        ! ... force the spin symmetry (CP wfc will be nspin=2 !!!)
        !
        IF ( ibnd .gt. num_bands ) THEN
          ibnd_ = ibnd - num_bands
        ELSE
          ibnd_ = ibnd
        ENDIF
        !
        evc_g(:) = ( 0.D0, 0.D0 )
        CALL mergewf( evc(:,ibnd_,ir), evc_g, npwxcp, ig_l2g_cp, &
                      mpime, nproc, ionode_id, intra_bgrp_comm )
        !
        IF ( ionode ) THEN
          !
          WRITE( cp_unit ) ( evc_g(ipw), ipw=1,npw_g )
          !
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    IF ( ionode ) CLOSE( cp_unit )
    !
    !
  END SUBROUTINE write_wannier_cp
  !
  !
END MODULE cp_files
