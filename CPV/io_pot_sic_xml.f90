!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE io_pot_sic_xml
   !----------------------------------------------------------------------------
   !
   USE kinds, ONLY: DP
   USE xml_io_base, ONLY: create_directory, write_pot_xml, read_pot_xml, &
                          restart_dir
   !
   PRIVATE
   !
   PUBLIC :: write_pot_sic, read_pot_sic

   INTERFACE write_pot_sic
      MODULE PROCEDURE write_pot_sic_realspace
      MODULE PROCEDURE write_pot_sic_reciprocal
   END INTERFACE

   INTERFACE read_pot_sic
      MODULE PROCEDURE read_pot_sic_realspace
      MODULE PROCEDURE read_pot_sic_reciprocal
   END INTERFACE

CONTAINS

   !------------------------------------------------------------------------
   SUBROUTINE write_pot_sic_realspace(pot, extension, field_specifier)
      !------------------------------------------------------------------------
      !
      ! ... this routine writes the charge-density in xml format into the
      ! ... '.save' directory
      ! ... the '.save' directory is created if not already present
      !
      USE io_files, ONLY: outdir, prefix
      USE gvecp, ONLY: ngm
      USE fft_base, ONLY: dfftp
      USE io_global, ONLY: ionode
      USE mp_global, ONLY: intra_pool_comm, inter_pool_comm
      USE control_flags, ONLY: ndw
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)           :: pot(dfftp%nnrx)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: field_specifier
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      REAL(DP), ALLOCATABLE :: potaux(:)
      !
      !
      ext = ' '
      !
      !dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      dirname = restart_dir(outdir, ndw)
      !
      CALL create_directory(dirname)
      !
      IF (PRESENT(extension)) ext = '.'//TRIM(extension)
      !
      IF (PRESENT(extension)) THEN
         file_base = TRIM(dirname)//'/'//TRIM(field_specifier)//TRIM(ext)
      ELSE
         file_base = TRIM(dirname)//'/sic_potential'//TRIM(ext)
      END IF
      !
      CALL write_pot_xml(file_base, pot(:), dfftp%nr1, dfftp%nr2, &
                         dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                         ionode, intra_pool_comm, inter_pool_comm)
      RETURN
      !
   END SUBROUTINE write_pot_sic_realspace
   !
   SUBROUTINE write_pot_sic(pot, extension, field_specifier)

      ! This routine writes a potential (stored in reciprocal space) to file (in real space)

      USE fft_base, ONLY: dfftp
      use gvecp, only: ngm
      use cp_interfaces, only: invfft

      implicit none

      ! Arguments
      COMPLEX(DP), INTENT(IN)           :: pot(ngm)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: field_specifier
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension

      ! Local variables
      COMPLEX(DP)           :: psi(dfftp%nnrx)
      REAL(DP)              :: realpot(dfftp%nnrx)

      ! Transform potential to real-space
      call rho2psi('Dense', psi, dfftp%nnr, pot, ngm)
      call invfft('Dense', psi, dfftp)
      realpot = dble(psi)

      ! Write to file
      call write_pot_sic_realspace(realpot, extension, field_specifier)

   END SUBROUTINE write_pot_sic

   !------------------------------------------------------------------------
   SUBROUTINE read_pot_sic_realspace(pot, extension)
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the effective potential in xml format from the
      ! ... files saved into the '.save' directory
      !
      USE io_files, ONLY: tmp_dir, prefix
      USE fft_base, ONLY: dfftp
      USe gvecp, ONLY: ngm
      USE io_global, ONLY: ionode
      USE mp_global, ONLY: intra_pool_comm, inter_pool_comm
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(OUT)          :: pot(dfftp%nnrx)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension
      !
      CHARACTER(LEN=256)    :: dirname, file_base
      CHARACTER(LEN=256)    :: ext
      !
      ext = ' '
      !
      dirname = TRIM(tmp_dir)//TRIM(prefix)//'.save'
      !
      IF (PRESENT(extension)) ext = '.'//TRIM(extension)
      !
      file_base = TRIM(dirname)//'/sic-potential'//TRIM(ext)
      !
      CALL read_pot_xml(file_base, pot(:), dfftp%nr1, dfftp%nr2, &
                        dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%ipp, dfftp%npp, &
                        ionode, intra_pool_comm, inter_pool_comm)
      !
      RETURN
      !
   END SUBROUTINE read_pot_sic_realspace
   !
   SUBROUTINE read_pot_sic(pot, extension)

      ! This routine reads a potential from file (in real space) and stores it (in reciprocal space)

      USE fft_base, ONLY: dfftp
      use gvecp, only: ngm
      use cp_interfaces, only: fwfft

      implicit none

      ! Arguments
      COMPLEX(DP), INTENT(OUT)          :: pot(ngm)
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extension

      ! Local variables
      COMPLEX(DP)           :: psi(dfftp%nnrx)
      REAL(DP)              :: realpot(dfftp%nnrx)

      ! Read from file
      call read_pot_sic_realspace(realpot, extension)

      ! Transform to reciprocal space
      psi = realpot(:)
      call fwfft('Dense', psi, dfftp)
      call psi2rho('Dense', psi, dfftp%nnr, pot, ngm)

   END SUBROUTINE read_pot_sic

END MODULE io_pot_sic_xml
