!
! Copyright (C) 2007-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Non-Koopmans method
! Developed and implemented by I. Dabo
!      (Universite Paris-Est, Ecole des Ponts, ParisTech)
! Further developed and optimized by Andrea Ferretti
!      (MIT, University of Oxford)
!
#include "f_defs.h"

!-----------------------------------------------------------------------
   subroutine nksic_potential(nbsp, nx, c, f_diag, bec, becsum, &
                              deeq_sic, ispin, iupdwn, nupdwn, &
                              rhor, rhoc, wtot_realspace, wtot_reciprocal, vsic_realspace, vsic_reciprocal, do_wxd_, pink, nudx, &
                              wfc_centers, wfc_spreads, &
                              icompute_spread, is_empty)
!-----------------------------------------------------------------------
!
! ....calculate orbital dependent potentials,
!     following the Non-Koopmans' (NK) scheme,
!     but also Perdew-Zunger (PZ),
!     Non-Koopmans' integral definition (NKI),
!     Non-Joopmans on Perdew Zunger (PZNK)
!
      use kinds, only: dp
      use gvecp, only: ngm
      use gvecw, only: ngw
      use grid_dimensions, only: nnrx
      USE electrons_base, ONLY: nspin
      use funct, only: dft_is_gradient
      use nksic, only: orb_rhor, wxdsic_realspace, wxdsic_reciprocal, &
                       wrefsic, rhoref, rhobar, &
                       do_nk, do_nki, do_pz, do_nkpz, &
                       do_nkipz, do_pz_renorm, &
                       grhobar, fion_sic, pzalpha => odd_alpha, &
                       kfact, upsilonkin, upsilonw, edens, &
                       taukin, tauw, valpsi, odd_alpha, nkscalfact
      use nksic, only: epsi2 => epsi2_cutoff_renorm
      use nksic_corrections, only: nksic_correction_nk, nksic_correction_nki, &
                       nksic_correction_nkpz, nksic_correction_nkipz, &
                       nksic_correction_pz
      use ions_base, only: nat
      use control_flags, only: gamma_only, do_wf_cmplx
      use uspp_param, only: nhm
      use cp_interfaces, only: nksic_get_orbitalrho
      use input_parameters, only: draw_pot, pot_number, odd_nkscalfact  !added:linh draw vsic_realspace potentials
      use io_pot_sic_xml, only: write_pot_sic  !added:linh draw vsic_realspace potentials
      USE io_global, ONLY: stdout
      use core, ONLY: nlcc_any
      use twin_types
      !
      implicit none
      !
      integer, intent(in)  :: nbsp, nx, nudx
      complex(dp), intent(in)  :: c(ngw, nx)
      type(twin_matrix), intent(in)  :: bec!(nkb,nbsp) !modified:giovanni
      real(dp), intent(in)  :: becsum(nhm*(nhm + 1)/2, nat, nspin)
      integer, intent(in)  :: ispin(nx)
      integer, intent(in)  :: iupdwn(nspin), nupdwn(nspin)
      real(dp), intent(in)  :: f_diag(nx)
      real(dp)                 :: rhor(nnrx, nspin)
      real(dp), intent(in)  :: rhoc(nnrx)
      real(dp), intent(out) :: vsic_realspace(nnrx, nx), wtot_realspace(nnrx, 2)
      complex(dp), intent(out) :: vsic_reciprocal(ngm, nx), wtot_reciprocal(ngm, 2)
      real(dp), intent(out) :: deeq_sic(nhm, nhm, nat, nx)
      logical, intent(in)  :: do_wxd_
      real(dp), intent(out) :: pink(nx)
      logical  :: icompute_spread
      real(DP) :: wfc_centers(4, nudx, nspin)
      real(DP) :: wfc_spreads(nudx, nspin, 2)
      logical  :: is_empty
      !
      ! local variables
      !
      integer  :: i, j, jj, ibnd, ir
      real(dp) :: focc, pinkpz, shart
      real(dp), allocatable :: vsicpz_realspace(:), rhor_nocc(:, :)
      complex(dp), allocatable :: vsicpz_reciprocal(:)
      complex(dp), allocatable :: rhobarg(:, :)
      logical  :: lgam, is_empty_
      !
      ! main body
      !
      CALL start_clock('nksic_drv')
      lgam = gamma_only .and. .not. do_wf_cmplx
      !
      is_empty_ = is_empty
      !
      ! compute potentials
      !
      if (dft_is_gradient()) then
         allocate (rhobarg(ngm, 2))
         !write(6,*) "allocated rhobarg"
      else
         allocate (rhobarg(1, 1))
      end if
      !
      if (nlcc_any) then
         !
         allocate (rhor_nocc(nnrx, nspin))
         rhor_nocc(:, :) = rhor(:, :)
         !
         ! add core charge
         !
         call add_cc_rspace(rhoc, rhor)
         !
      end if
      !
      if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
         wtot_realspace = 0.0_dp
         wtot_reciprocal = 0.0_dp
      end if
      !
      if (do_nkpz .or. do_nkipz) then
         allocate (vsicpz_realspace(nnrx))
         allocate (vsicpz_reciprocal(ngm))
         vsicpz_realspace = 0.0_dp
         vsicpz_reciprocal = 0.0_dp
      end if
      !
      pink = 0.0_dp
      vsic_realspace = 0.0_dp
      vsic_reciprocal = 0.0_dp
      !
      !
      ! if using pz_renorm factors, compute here tauw and upsilonw
      !
      if (do_pz_renorm) THEN
         !
         edens = 0.d0
         taukin = 0.d0
         tauw = 0.d0
         !
      END IF
      !
      ! loop over bands (2 ffts at the same time)
      !
      !
      do j = 1, nbsp, 2
         !
         ! compute orbital densities
         ! n odd => c(:,n+1) is already set to zero
         !
         call nksic_get_orbitalrho(ngw, nnrx, bec, ispin, nbsp, &
                                   c(:, j), c(:, j + 1), orb_rhor, j, j + 1, lgam) !warning:giovanni need modification
         !
         ! compute centers and spreads of nksic or pz
         ! minimizing orbitals
         !
         if (icompute_spread) then
            !
            call compute_nksic_centers(nnrx, nx, nudx, nbsp, nspin, iupdwn, &
                                       nupdwn, ispin, orb_rhor, wfc_centers, wfc_spreads, j, j + 1)
            !
         end if
         !
         shart = 0.d0
         !
         ! compute orbital potentials
         !
         inner_loop: do jj = 1, 2
            !
            i = j + jj - 1
            !
            ! this condition is important when n is odd
            !
            if (i > nbsp) exit inner_loop
            !
            ibnd = i
            !
            if (nspin == 2) then
               !
               if (i >= iupdwn(2)) ibnd = i - iupdwn(2) + 1
               !
            end if
            !
            ! note: iupdwn(2) is set to zero if nspin = 1
            !
            focc = f_diag(i)*DBLE(nspin)/2.0d0
            !
            ! compute parameters needed for PZ-renormalization
            !
            IF (do_pz_renorm) THEN
               !
               !call nksic_get_taukin_pz( focc, nspin, ispin(i), orb_rhor(:,jj), &
               !taukin, ibnd, 1)
               !
               IF (ibnd == 1) THEN
                  !
                  IF (nspin == 1) THEN
                     !
                     !call nksic_get_taukin_pz( 0.5d0, nspin, ispin(i), &
                     !                   rhor(:,1), tauw, ibnd, nupdwn(ispin(i)))
                     !
                  ELSE IF (nspin == 2) THEN
                     !
                     !call nksic_get_taukin_pz( 1.d0, nspin, ispin(i), &
                     !                   rhor(:,ispin(i)), tauw, ibnd, nupdwn(ispin(i)))
                     !
                  END IF
                  !
               END IF
               !
            END IF
            !
            ! define rhoref and rhobar
            !
            call nksic_get_rhoref(i, nnrx, ispin(i), nspin, &
                                  focc, rhor, orb_rhor(:, jj), &
                                  rhoref, rhobar, rhobarg, grhobar)

            !
            ! compute nk pieces to build the potentials and the energy
            !
            if (do_nk .or. do_nkpz) then
               !
               call nksic_correction_nk(focc, ispin(i), orb_rhor(:, jj), &
                                        rhor, rhoref, rhobar, rhobarg, grhobar, &
                                        vsic_realspace(:, i), wxdsic_realspace, wrefsic, do_wxd_, &
                                        pink(i), ibnd, shart)
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
               ! here information is accumulated over states
               ! (wtot_realspace is added in the next loop)
               !
               wtot_realspace(1:nnrx, 1:2) = wtot_realspace(1:nnrx, 1:2) + wxdsic_realspace(1:nnrx, 1:2)
               !
               ! ths sic potential is partly updated here to save some memory
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + wrefsic(1:nnrx) &
                                 - wxdsic_realspace(1:nnrx, ispin(i))
               !
            end if
            !
            ! compute nkpz pieces to build the potential and the energy
            !
            if (do_nkpz) then
               !
               call nksic_correction_nkpz(focc, orb_rhor(:, jj), vsicpz_realspace, &
                                          wrefsic, pinkpz, ibnd, ispin(i))
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + vsicpz_realspace(1:nnrx) &
                                 + wrefsic(1:nnrx)
               !
               pink(i) = pink(i) + pinkpz
               !
            end if
            !
            ! compute pz potentials and energy
            !
            if (do_pz) then
               !
               call nksic_correction_pz(focc, ispin(i), orb_rhor(:, jj), &
                                        vsic_realspace(:, i), vsic_reciprocal(:, i), pink(i), pzalpha(i), ibnd, shart)
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
               if (do_pz_renorm) then
                  !
                  do ir = 1, nnrx
                     !
                     edens(ir, ispin(i)) = edens(ir, ispin(i)) + pink(i)*(orb_rhor(ir, jj) + epsi2)**(kfact + 1.)
                     !
                  end do
                  !
               end if
               !
            end if
            !
            ! compute nki pieces to build the potentials and the energy
            !
            if (do_nki .or. do_nkipz) then
               !
               call nksic_correction_nki(focc, ispin(i), orb_rhor(:, jj), &
                                         rhor, rhoref, rhobar, rhobarg, grhobar, &
                             vsic_realspace(:, i), vsic_reciprocal(:, i), wxdsic_realspace, wxdsic_reciprocal, do_wxd_, pink(i), ibnd, shart, is_empty_)
               !
               ! here information is accumulated over states
               ! (wtot_realspace is added in the next loop)
               !
               wtot_realspace(1:nnrx, 1:2) = wtot_realspace(1:nnrx, 1:2) + wxdsic_realspace(1:nnrx, 1:2)
               wtot_reciprocal(1:ngm, 1:2) = wtot_reciprocal(1:ngm, 1:2) + wxdsic_reciprocal(1:ngm, 1:2)
               !
               ! ths sic potential is partly updated here to save some memory
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) - wxdsic_realspace(1:nnrx, ispin(i))
               vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i) - wxdsic_reciprocal(1:ngm, ispin(i))
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
            end if

            if (do_nkipz) then
               !
               call nksic_correction_nkipz(focc, ispin(i), orb_rhor(:, jj), vsicpz_realspace, vsicpz_reciprocal, &
                                           pinkpz, ibnd, shart, is_empty_)
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + vsicpz_realspace(1:nnrx)
               vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i) + vsicpz_reciprocal(1:ngm)
               !
               pink(i) = pink(i) + pinkpz
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
            end if
            !
            ! take care of spin symmetry
            !
            if (.not. do_pz_renorm) then
               !
               if (.not. is_empty_) then
                  !
                  pink(i) = 2.d0*pink(i)/nspin
                  !
               else
                  !
                  pink(i) = 2.d0*pink(i)/nspin
                  !
               end if
               !
            end if
            !
            if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
               !
               if (nspin == 1) then
                  !
                  wtot_realspace(1:nnrx, 1) = wtot_realspace(1:nnrx, 1) + wxdsic_realspace(1:nnrx, 2)
                  wtot_reciprocal(1:ngm, 1) = wtot_reciprocal(1:ngm, 1) + wxdsic_reciprocal(1:ngm, 2)
                  !
                  wtot_realspace(1:nnrx, 2) = wtot_realspace(1:nnrx, 2) + wxdsic_realspace(1:nnrx, 1)
                  wtot_reciprocal(1:ngm, 2) = wtot_reciprocal(1:ngm, 2) + wxdsic_reciprocal(1:ngm, 1)
                  !
               end if
               !
            end if
            !
         end do inner_loop
         !
      end do
      !
      ! Switch off the icompute_spread flag if present
      !
      IF (icompute_spread) THEN
         !
         icompute_spread = .false.
         !
      END IF
      !
      ! now wtot_realspace is completely built and can be added to vsic
      !
      if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
         !
         do i = 1, nbsp
            !
            vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + wtot_realspace(1:nnrx, ispin(i))
            vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i) + wtot_reciprocal(1:ngm, ispin(i))
            !
         end do
         !
      end if
      !
      ! computing orbital dependent alpha
      !
      if (odd_nkscalfact) then
         !
         do j = 1, nbsp, 2
            !
            inner_loop_odd_alpha: do jj = 1, 2
               !
               i = j + jj - 1
               !
               if (i > nbsp) exit inner_loop_odd_alpha
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i)*odd_alpha(i)/nkscalfact
               vsic_reciprocal(:, i) = vsic_reciprocal(:, i)*odd_alpha(i)/nkscalfact
               !
               valpsi(i, :) = valpsi(i, :)*pink(i)/nkscalfact
               !
               pink(i) = pink(i)*odd_alpha(i)/nkscalfact
               !
            end do inner_loop_odd_alpha
            !
         end do
         !
      end if
      !
      if (do_pz_renorm) then
         !
         do j = 1, nbsp, 2
            !
            call nksic_get_orbitalrho(ngw, nnrx, bec, ispin, nbsp, &
                                      c(:, j), c(:, j + 1), orb_rhor, j, j + 1, lgam)
            !
            inner_loop_renorm: do jj = 1, 2
               !
               i = j + jj - 1
               !
               if (i > nbsp) exit inner_loop_renorm
               !
               ibnd = i
               !
               focc = f_diag(i)*DBLE(nspin)/2.0d0
               !
               if (nspin == 2) then
                  !
                  if (i >= iupdwn(2)) ibnd = i - iupdwn(2) + 1
                  !
               end if
               !
               call nksic_get_pz_factor(nspin, ispin(i), orb_rhor(:, jj), rhor, &
                                        taukin, tauw, pzalpha(i), ibnd, kfact)
               !
               !
               ! update vsic_realspace with factor here: it works for pz, will it work for
               ! nk-type functionals?
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i)*pzalpha(i)
               vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i)*pzalpha(i)
               !
               call nksic_get_pzfactor_potential(focc, nspin, ispin(i), rhor, orb_rhor(:, jj), &
                                            pink(i), taukin, tauw, edens, upsilonkin, upsilonw, vsic_realspace(:, i), &
                                            vsic_reciprocal(:, i), pzalpha(i), ibnd, kfact)
               !
               pink(i) = pink(i)*pzalpha(i)
               !
               if (.not. is_empty_) then
                  !
                  pink(i) = f_diag(i)*pink(i)
                  !
               else
                  !
                  pink(i) = 2.d0*pink(i)/nspin
                  !
               end if
               !
            end do inner_loop_renorm
            !
         end do
         !
      end if
      !
      if (draw_pot) then !added:linh draw vsic_realspace potentials
         !
         write (stdout, *) "I am writing out vsic", nbsp
         !
         do i = 1, nbsp
            !
            if (i == pot_number) call write_pot_sic(vsic_realspace(:, i))
            !
         end do
         !
      end if !added:linh draw vsic_realspace potentials
      !
      if (allocated(vsicpz_realspace)) deallocate (vsicpz_realspace)
      if (allocated(vsicpz_reciprocal)) deallocate (vsicpz_reciprocal)
      !
      ! USPP:
      ! compute corrections to the D coefficients of the pseudopots
      ! due to vsic_realspace(r, i) in the case of orbital dependent functionals.
      ! The corresponding contributions to the forces are computed.
      !
      ! IMPORTANT: the following call makes use of newd.
      !            It must be done before we call newd for the
      !            total potentials, because deeq is overwritten at every call
      !
      fion_sic(:, :) = 0.0d0
      !
      IF (nhm > 0) then
         !
         deeq_sic(:, :, :, :) = 0.0d0
         !
         DO i = 1, nbsp
            !
            CALL nksic_newd(i, nnrx, ispin(i), nspin, vsic_realspace(:, i), vsic_reciprocal(:, i), nat, nhm, &
                            becsum, fion_sic, deeq_sic(:, :, :, i)) !this is for ultrasoft! watch out! warning:giovanni this has to be modified in order to run ultrasoft
            !
         END DO
         !
      END IF
      !
      deallocate (rhobarg)
      !
      if (nlcc_any) then
         !
         rhor(:, :) = rhor_nocc(:, :)
         deallocate (rhor_nocc)
         !
      end if
      call ebl_check(vsic_realspace(:, 1), vsic_reciprocal(:, 1))
      !
      CALL stop_clock('nksic_drv')
      return
      !
!-----------------------------------------------------------------------
   end subroutine nksic_potential
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_get_orbitalrho_real(ngw, nnrx, bec, ispin, nbsp, &
                                        c1, c2, orb_rhor, i1, i2)
!-----------------------------------------------------------------------
!
! Computes orbital densities on the real (not smooth) grid
!
      use kinds, only: dp
      use constants, only: ci
      use cp_interfaces, only: fwfft, invfft, calrhovan
      use fft_base, only: dffts, dfftp
      use cell_base, only: omega
      use gvecp, only: ngm
      use gvecs, only: ngs, nps, nms
      use recvecs_indexes, only: np, nm
      use smooth_grid_dimensions, only: nnrsx
      use cp_main_variables, only: irb, eigrb
      use uspp_param, only: nhm
      use electrons_base, only: nspin
      use ions_base, only: nat
      use uspp, only: okvan, nkb
      !
      implicit none

      !
      ! input/output vars
      !
      integer, intent(in) :: ngw, nnrx, i1, i2
      integer, intent(in) :: nbsp, ispin(nbsp)
      real(dp), intent(in) :: bec(nkb, nbsp)
      complex(dp), intent(in) :: c1(ngw), c2(ngw)
      real(dp), intent(out) :: orb_rhor(nnrx, 2)

      !
      ! local vars
      !
      character(20) :: subname = 'nksic_get_orbitalrho'
      integer       :: ir, ig, ierr
      real(dp)      :: sa1
      complex(dp)   :: fm, fp
      complex(dp), allocatable :: psis(:), psi(:)
      complex(dp), allocatable :: orb_rhog(:, :)
      real(dp), allocatable :: orb_rhos(:)
      real(dp), allocatable :: rhovan(:, :, :)
      real(dp), allocatable :: rhovanaux(:, :, :)
      !
      !====================
      ! main body
      !====================
      !
      call start_clock('nk_orbrho')

      !
      if (okvan) then
         !
         allocate (rhovan(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovan', abs(ierr))
         allocate (rhovanaux(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovanaux', abs(ierr))
         !
      end if
      !
      allocate (psi(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi', abs(ierr))
      !
      allocate (orb_rhog(ngm, 2), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating orb_rhog', abs(ierr))

      sa1 = 1.0d0/omega

      !
      ! check whether it is necessary to
      ! deal with the smooth and dense grids separately
      !
      if (nnrsx == nnrx) then
         !
         ! This case should be the case when using NCPP
         !
         CALL c2psi(psi, nnrx, c1, c2, ngw, 2)
         !
         !CALL invfft('Wave', psi, dffts )
         CALL invfft('Dense', psi, dfftp)
         !
         ! computing the orbital charge in real space on the full grid
         !
         do ir = 1, nnrx
            !
            orb_rhor(ir, 1) = sa1*(DBLE(psi(ir)))**2
            orb_rhor(ir, 2) = sa1*(AIMAG(psi(ir)))**2
            !
         end do
         !
      else
         !
         ! this is the general case,
         ! normally used with USPP
         !
         allocate (psis(nnrsx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psis', abs(ierr))
         allocate (orb_rhos(2), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating orb_rhos', abs(ierr))
         !
         CALL c2psi(psis, nnrsx, c1, c2, ngw, 2)
         !
         CALL invfft('Wave', psis, dffts)
         !
         ! computing the orbital charge
         ! in real space on the smooth grid
         !
         do ir = 1, nnrsx
            !
            orb_rhos(1) = sa1*(DBLE(psis(ir)))**2
            orb_rhos(2) = sa1*(AIMAG(psis(ir)))**2
            !
            psis(ir) = CMPLX(orb_rhos(1), orb_rhos(2))
         end do
         !
         ! orbital charges are taken to the G space
         !
         CALL fwfft('Smooth', psis, dffts)
         !
         do ig = 1, ngs
            !
            fp = psis(nps(ig)) + psis(nms(ig))
            fm = psis(nps(ig)) - psis(nms(ig))
            orb_rhog(ig, 1) = 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            orb_rhog(ig, 2) = 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
            !
         end do
         !
         psi(:) = (0.d0, 0.d0)
         do ig = 1, ngs
            !
            psi(nm(ig)) = CONJG(orb_rhog(ig, 1)) + ci*CONJG(orb_rhog(ig, 2))
            psi(np(ig)) = orb_rhog(ig, 1) + ci*orb_rhog(ig, 2)
            !
         end do
         !
         call invfft('Dense', psi, dfftp)
         !
         do ir = 1, nnrx
            !
            orb_rhor(ir, 1) = DBLE(psi(ir))
            orb_rhor(ir, 2) = AIMAG(psi(ir))
         end do

         deallocate (psis)
         deallocate (orb_rhos)

      end if

      !
      ! add Vanderbilt contribution to orbital density
      !
      if (okvan) then
         !
         rhovan(:, :, :) = 0.0d0
         !
         if (nspin == 2) then
            !
            if (i1 <= nbsp) then
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1))
            end if
            !
            if (i2 <= nbsp) then
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 2) = rhovanaux(:, :, ispin(i2))
            end if
            !
            call rhov(irb, eigrb, rhovan, orb_rhog, orb_rhor, .true.)
         else
            !
            if (i1 <= nbsp) then
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1))*0.5d0 ! 1/2 factor since rhovanaux is counted twice in the case nspin=2
               !
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 1), orb_rhor(:, 1), .true.)
               !
            end if
            !
            if (i2 <= nbsp) then
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i2))*0.5d0 ! 1/2 factor since rhovanaux is counted twice in the case nspin=2
               !
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 2), orb_rhor(:, 2), .true.)
               !
            end if
            !
         end if
         !
      end if
      !
      deallocate (psi)
      deallocate (orb_rhog)
      !
      if (okvan) then
         deallocate (rhovan)
         deallocate (rhovanaux)
      end if
      !
      call stop_clock('nk_orbrho')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_get_orbitalrho_real
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_get_orbitalrho_twin_non_ortho(ngw, nnrx, bec, becdual, ispin, nbsp, &
                                                  c1, c2, c1dual, c2dual, orb_rhor, i1, i2, lgam)
!-----------------------------------------------------------------------
!
! Computes orbital densities on the real (not smooth) grid
!
      use kinds, only: dp
      use constants, only: ci
      use cp_interfaces, only: fwfft, invfft, calrhovan
      use fft_base, only: dffts, dfftp
      use cell_base, only: omega
      use gvecp, only: ngm
      use gvecs, only: ngs, nps, nms
      use recvecs_indexes, only: np, nm
      use smooth_grid_dimensions, only: nnrsx
      use cp_main_variables, only: irb, eigrb
      use uspp_param, only: nhm
      use electrons_base, only: nspin
      use ions_base, only: nat
      use uspp, only: okvan
      use twin_types
      !
      implicit none

      !
      ! input/output vars
      !
      integer, intent(in) :: ngw, nnrx, i1, i2
      integer, intent(in) :: nbsp, ispin(nbsp)
      type(twin_matrix) :: bec, becdual !(nkb, nbsp)
      complex(dp), intent(in) :: c1(ngw), c2(ngw), c1dual(ngw), c2dual(ngw)
      real(dp), intent(out) :: orb_rhor(nnrx, 2)
      logical :: lgam
      !
      ! local vars
      !
      character(20) :: subname = 'nksic_get_orbitalrho'
      integer       :: ir, ig, ierr
      real(dp)      :: sa1
      complex(dp)   :: fm, fp
      complex(dp), allocatable :: psis1(:), psis2(:), psi1(:), psi2(:), &
                                  psi1d(:), psi2d(:)
      complex(dp), allocatable :: orb_rhog(:, :)
      real(dp), allocatable :: orb_rhos(:)
      real(dp), allocatable :: rhovan(:, :, :)
      real(dp), allocatable :: rhovanaux(:, :, :)
      !
      !====================
      ! main body
      !====================
      !
      call start_clock('nksic_orbrho')

      !
      if (okvan) then
         !
         allocate (rhovan(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovan', abs(ierr))
         allocate (rhovanaux(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovanaux', abs(ierr))
         !
      end if
      !
      allocate (psi1(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi1', abs(ierr))

      allocate (psi1d(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi1d', abs(ierr))

      !
      if (.not. lgam) then
         allocate (psi2(nnrx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psi2', abs(ierr))
         allocate (psi2d(nnrx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psi2d', abs(ierr))
      end if
      !
      allocate (orb_rhog(ngm, 2), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating orb_rhog', abs(ierr))
      sa1 = 1.0d0/omega
      !
      ! check whether it is necessary to
      ! deal with the smooth and dense grids separately
      !
      if (nnrsx == nnrx) then
         !
         ! This case should be the one when using NCPP
         !
         if (lgam) then
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 2)
            CALL c2psi(psi1d, nnrx, c1dual, c2dual, ngw, 2)
         else
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 0)
            CALL c2psi(psi2, nnrx, c2, c1, ngw, 0)
            CALL c2psi(psi1d, nnrx, c1dual, c2dual, ngw, 0)
            CALL c2psi(psi2d, nnrx, c2dual, c1dual, ngw, 0)
         end if
         !
         CALL invfft('Dense', psi1, dfftp)
         CALL invfft('Dense', psi1d, dfftp)
         !
         !
         if (.not. lgam) then
            CALL invfft('Dense', psi2, dfftp)
            CALL invfft('Dense', psi2d, dfftp)
         end if
         !
         ! computing the orbital charge in real space on the full grid
         !
         if (lgam) then
            do ir = 1, nnrx
               !
               orb_rhor(ir, 1) = sa1*DBLE(psi1(ir))*DBLE(psi1d(ir))
               orb_rhor(ir, 2) = sa1*AIMAG(psi1(ir))*AIMAG(psi1d(ir))
               !
            end do
         else
            do ir = 1, nnrx
               !
               orb_rhor(ir, 1) = sa1*DBLE(CONJG(psi1d(ir))*psi1(ir))
               orb_rhor(ir, 2) = sa1*DBLE(CONJG(psi2d(ir))*psi2(ir))
               !
            end do
         end if
         !
      else
         !
         ! this is the general case,
         ! normally used with USPP
         !

         allocate (psis1(nnrsx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psis1', abs(ierr))
         if (.not. lgam) then
            allocate (psis2(nnrsx), stat=ierr)
            if (ierr /= 0) call errore(subname, 'allocating psis2', abs(ierr))
         end if

         allocate (orb_rhos(2), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating orb_rhos', abs(ierr))
         !
         if (lgam) then
            CALL c2psi(psis1, nnrsx, c1, c2, ngw, 2)
         else
            CALL c2psi(psis1, nnrsx, c1, c2, ngw, 0)
            CALL c2psi(psis2, nnrsx, c2, c1, ngw, 0)
         end if
         !

         CALL invfft('Wave', psis1, dffts)
         !
         if (.not. lgam) then
            CALL invfft('Wave', psis2, dffts)
         end if
         !
         ! computing the orbital charge
         ! in real space on the smooth grid
         !
         if (lgam) then
            do ir = 1, nnrsx
               !
               orb_rhos(1) = sa1*((DBLE(psis1(ir)))**2)
               orb_rhos(2) = sa1*((AIMAG(psis1(ir)))**2)
               !
               psis1(ir) = CMPLX(orb_rhos(1), orb_rhos(2))
            end do
         else
            do ir = 1, nnrsx
               !
               orb_rhos(1) = sa1*((DBLE(psis1(ir)))**2 + (AIMAG(psis1(ir)))**2)
               orb_rhos(2) = sa1*((DBLE(psis2(ir)))**2 + (AIMAG(psis2(ir)))**2)
               !
               psis1(ir) = CMPLX(orb_rhos(1), orb_rhos(2)) !!!### comment for k points
!                 psis1( ir )  = cmplx( orb_rhos(1), 0.d0) !!!### uncomment for k points
!                 psis2( ir )  = cmplx( orb_rhos(2), 0.d0) !!!### uncomment for k points
            end do
         end if
!           write(6,*) "psis", psis1 !added:giovanni:debug
         !
         ! orbital charges are taken to the G space
         !

         CALL fwfft('Smooth', psis1, dffts)
!           IF(.not.lgam) THEN !  !!!### uncomment for k points
!               CALL fwfft('Smooth',psis2, dffts ) !!!### uncomment for k points
!           ENDIF !!!### uncomment for k points
         !
!         IF(lgam) then !!!### uncomment for k points
         do ig = 1, ngs
            !
            fp = psis1(nps(ig)) + psis1(nms(ig))
            fm = psis1(nps(ig)) - psis1(nms(ig))
            orb_rhog(ig, 1) = 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            orb_rhog(ig, 2) = 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
            !
         end do
!           else !!!### uncomment for k points
!             do ig = 1, ngs !!!### uncomment for k points
         !
!                 fp=psis1(nps(ig)) !!!### uncomment for k points
!                 fm=psis2(nps(ig)) !!!### uncomment for k points
!                 orb_rhog(ig,1)=fp !!!### uncomment for k points
!                 orb_rhog(ig,2)=fm !!!### uncomment for k points
         !
!             enddo !!!### uncomment for k points
!           endif !!!### uncomment for k points
         !
         psi1(:) = CMPLX(0.d0, 0.d0)
!           if(lgam) then !!!### uncomment for k points
         do ig = 1, ngs
            !
            psi1(nm(ig)) = CONJG(orb_rhog(ig, 1)) &
                           + ci*CONJG(orb_rhog(ig, 2))
            psi1(np(ig)) = orb_rhog(ig, 1) + ci*orb_rhog(ig, 2)
            !
         end do
!           else !!!### uncomment for k points
!             do ig=1,ngs !!!### uncomment for k points
         !
!                 psi1(nm(ig)) = conjg( orb_rhog(ig,1) ) &
!                               +ci*conjg( orb_rhog(ig,2) )
!                 psi1(np(ig)) = orb_rhog(ig,1) +ci*orb_rhog(ig,2)  !!!### uncomment for k points
         !
!             enddo !!!### uncomment for k points
!           endif !!!### uncomment for k points
         !
         call invfft('Dense', psi1, dfftp)
         !
         do ir = 1, nnrx
            !
            orb_rhor(ir, 1) = DBLE(psi1(ir))
            orb_rhor(ir, 2) = AIMAG(psi1(ir))
         end do

         deallocate (psis1)
         if (.not. lgam) then
            deallocate (psis2)
         end if

         deallocate (orb_rhos)

      end if
!       write(6,*) "orb_rhog", orb_rhog !added:giovanni:debug
      !
      ! add Vanderbilt contribution to orbital density
      !
      if (okvan) then
         !
         rhovan(:, :, :) = 0.0d0
         !
         if (nspin == 2) then
            !
            if (i1 <= nbsp) then
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1))
            end if
            !
            if (i2 <= nbsp) then
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 2) = rhovanaux(:, :, ispin(i2))
            end if
            !
            call rhov(irb, eigrb, rhovan, orb_rhog, orb_rhor, lgam)
         else
            !
            if (i1 <= nbsp) then
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1))
               !
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 1), orb_rhor(:, 1), lgam)
            end if
            !
            if (i2 <= nbsp) then
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i2))
               !
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 2), orb_rhor(:, 2), lgam)
            end if
            !
         end if
         !
      end if
      !
!       write(6,*) "rhovan", rhovan(:,:,1) !added:giovanni:debug
!       stop
      deallocate (psi1)
      if (allocated(psi2)) then
         deallocate (psi2)
      end if
      if (allocated(psi1d)) then
         deallocate (psi1d)
      end if
      if (allocated(psi2d)) then
         deallocate (psi2d)
      end if

      deallocate (orb_rhog)
      !
      if (okvan) then
         deallocate (rhovan)
         deallocate (rhovanaux)
      end if
      !
      do ir = 1, nnrx
         if (orb_rhor(ir, 1) .lt. -1.d-3 .or. orb_rhor(ir, 2) .lt. -1.d-3) then
            write (6, *) "warning, negative density", orb_rhor(ir, 1), orb_rhor(ir, 2)
         end if
      end do
      !
      call stop_clock('nksic_orbrho')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_get_orbitalrho_twin_non_ortho
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_get_orbitalrho_twin(ngw, nnrx, bec, ispin, nbsp, &
                                        c1, c2, orb_rhor, i1, i2, lgam)
!-----------------------------------------------------------------------
!
! Computes orbital densities on the real (not smooth) grid
!
      use kinds, only: dp
      use constants, only: ci
      use cp_interfaces, only: fwfft, invfft, calrhovan
      use fft_base, only: dffts, dfftp
      use cell_base, only: omega
      use gvecp, only: ngm
      use gvecs, only: ngs, nps, nms
      use recvecs_indexes, only: np, nm
      use smooth_grid_dimensions, only: nnrsx
      use cp_main_variables, only: irb, eigrb
      use uspp_param, only: nhm
      use electrons_base, only: nspin
      use ions_base, only: nat
      use uspp, only: okvan
      use twin_types
      !
      implicit none
      !
      ! input/output vars
      !
      integer, intent(in) :: ngw, nnrx, i1, i2
      integer, intent(in) :: nbsp, ispin(nbsp)
      type(twin_matrix) :: bec !(nkb, nbsp)
      complex(dp), intent(in) :: c1(ngw), c2(ngw)
      real(dp), intent(out) :: orb_rhor(nnrx, 2)
      logical :: lgam
      !
      ! local vars
      !
      character(20) :: subname = 'nksic_get_orbitalrho'
      integer       :: ir, ig, ierr
      real(dp)      :: sa1
      complex(dp)   :: fm, fp
      complex(dp), allocatable :: psis1(:), psis2(:), psi1(:), psi2(:)
      complex(dp), allocatable :: orb_rhog(:, :)
      real(dp), allocatable :: orb_rhos(:)
      real(dp), allocatable :: rhovan(:, :, :)
      real(dp), allocatable :: rhovanaux(:, :, :)
      !
      !====================
      ! main body
      !====================
      !
      call start_clock('nksic_orbrho')

      !
      if (okvan) then
         !
         allocate (rhovan(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovan', abs(ierr))
         allocate (rhovanaux(nhm*(nhm + 1)/2, nat, nspin), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating rhovanaux', abs(ierr))
         !
      end if
      !
      allocate (psi1(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi1', abs(ierr))
      !
      if (.not. lgam) then
         allocate (psi2(nnrx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psi2', abs(ierr))
      end if
      !
      allocate (orb_rhog(ngm, 2), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating orb_rhog', abs(ierr))

      sa1 = 1.0d0/omega
      !
      ! check whether it is necessary to
      ! deal with the smooth and dense grids separately
      !
      if (nnrsx == nnrx) then
         !
         ! This case should be the one when using NCPP
         !
         if (lgam) then
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 2)
         else
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 0)
            CALL c2psi(psi2, nnrx, c2, c1, ngw, 0)
         end if
         !
         CALL invfft('Dense', psi1, dfftp)
         !
         if (.not. lgam) then
            CALL invfft('Dense', psi2, dfftp)
         end if
         !
         ! computing the orbital charge in real space on the full grid
         !
         if (lgam) then
            !
            do ir = 1, nnrx
               !
               orb_rhor(ir, 1) = sa1*((DBLE(psi1(ir)))**2)
               orb_rhor(ir, 2) = sa1*((AIMAG(psi1(ir)))**2)
               !
            end do
            !
         else
            !
            do ir = 1, nnrx
               !
               orb_rhor(ir, 1) = sa1*((abs(psi1(ir))))**2
               orb_rhor(ir, 2) = sa1*((abs(psi2(ir))))**2
               !
            end do
            !
         end if
         !
      else
         !
         ! this is the general case,
         ! normally used with USPP
         !

         allocate (psis1(nnrsx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psis1', abs(ierr))
         if (.not. lgam) then
            !
            allocate (psis2(nnrsx), stat=ierr)
            if (ierr /= 0) call errore(subname, 'allocating psis2', abs(ierr))
            !
         end if
         !
         allocate (orb_rhos(2), stat=ierr)
         !
         if (ierr /= 0) call errore(subname, 'allocating orb_rhos', abs(ierr))
         !
         if (lgam) then
            !
            CALL c2psi(psis1, nnrsx, c1, c2, ngw, 2)
            !
         else
            !
            CALL c2psi(psis1, nnrsx, c1, c2, ngw, 0)
            CALL c2psi(psis2, nnrsx, c2, c1, ngw, 0)
            !
         end if
         !
         CALL invfft('Wave', psis1, dffts)
         !
         if (.not. lgam) then
            !
            CALL invfft('Wave', psis2, dffts)
            !
         end if
         !
         ! computing the orbital charge
         ! in real space on the smooth grid
         !
         if (lgam) then
            !
            do ir = 1, nnrsx
               !
               orb_rhos(1) = sa1*((DBLE(psis1(ir)))**2)
               orb_rhos(2) = sa1*((AIMAG(psis1(ir)))**2)
               !
               psis1(ir) = CMPLX(orb_rhos(1), orb_rhos(2))
               !
            end do
            !
         else
            !
            do ir = 1, nnrsx
               !
               orb_rhos(1) = sa1*(abs(psis1(ir)))**2
               orb_rhos(2) = sa1*(abs(psis2(ir)))**2
               !
               psis1(ir) = CMPLX(orb_rhos(1), orb_rhos(2)) !!!### comment for k points
               !psis1( ir )  = cmplx( orb_rhos(1), 0.d0) !!!### uncomment for k points
               !psis2( ir )  = cmplx( orb_rhos(2), 0.d0) !!!### uncomment for k points
            end do
            !
         end if
!           write(6,*) "psis", psis1 !added:giovanni:debug
         !
         ! orbital charges are taken to the G space
         !

         CALL fwfft('Smooth', psis1, dffts)
!           IF(.not.lgam) THEN !  !!!### uncomment for k points
!               CALL fwfft('Smooth',psis2, dffts ) !!!### uncomment for k points
!           ENDIF !!!### uncomment for k points
         !
!           IF(lgam) then !!!### uncomment for k points
         do ig = 1, ngs
            !
            fp = psis1(nps(ig)) + psis1(nms(ig))
            fm = psis1(nps(ig)) - psis1(nms(ig))
            orb_rhog(ig, 1) = 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            orb_rhog(ig, 2) = 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
            !
         end do
!           else !!!### uncomment for k points
!               do ig = 1, ngs !!!### uncomment for k points
         !
!                   fp=psis1(nps(ig)) !!!### uncomment for k points
!                   fm=psis2(nps(ig)) !!!### uncomment for k points
!                   orb_rhog(ig,1)=fp !!!### uncomment for k points
!                   orb_rhog(ig,2)=fm !!!### uncomment for k points
         !
!               enddo !!!### uncomment for k points
!           endif !!!### uncomment for k points
         !
         psi1 = CMPLX(0.d0, 0.d0)
!           if(lgam) then !!!### uncomment for k points
         do ig = 1, ngs
            !
            psi1(nm(ig)) = CONJG(orb_rhog(ig, 1)) &
                           + ci*CONJG(orb_rhog(ig, 2))
            psi1(np(ig)) = orb_rhog(ig, 1) + ci*orb_rhog(ig, 2)
            !
         end do
!           else !!!### uncomment for k points
!               do ig=1,ngs !!!### uncomment for k points
         !
!                   psi1(nm(ig)) = conjg( orb_rhog(ig,1) ) &
!                                 +ci*conjg( orb_rhog(ig,2) )
!                   psi1(np(ig)) = orb_rhog(ig,1) +ci*orb_rhog(ig,2)  !!!### uncomment for k points
         !
!               enddo !!!### uncomment for k points
!           endif !!!### uncomment for k points
         !
         call invfft('Dense', psi1, dfftp)
         !
         do ir = 1, nnrx
            !
            orb_rhor(ir, 1) = DBLE(psi1(ir))
            orb_rhor(ir, 2) = AIMAG(psi1(ir))
            !
         end do

         deallocate (psis1)

         if (.not. lgam) then
            !
            deallocate (psis2)
            !
         end if

         deallocate (orb_rhos)

      end if
!       write(6,*) "orb_rhog", orb_rhog !added:giovanni:debug
      !
      ! add Vanderbilt contribution to orbital density
      !
      if (okvan) then
         !
         rhovan(:, :, :) = 0.0d0
         !
         if (nspin == 2) then
            !
            if (i1 <= nbsp) then
               !
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1))
               !
            end if
            !
            if (i2 <= nbsp) then
               !
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 2) = rhovanaux(:, :, ispin(i2))
               !
            end if
            !
            call rhov(irb, eigrb, rhovan, orb_rhog, orb_rhor, lgam)
            !
         else
            !
            if (i1 <= nbsp) then
               !
               call calrhovan(rhovanaux, bec, i1)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i1)) ! 0.5 to divide the factor f=2 which accounts for spin multiplicity inside calrhovan
               !
               write (6, *) "calling rhov inside nksic_get_orbitalrho"
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 1), orb_rhor(:, 1), lgam) !JUST-FOR-NOW ... do we need a factor of 0.5?
               !
            end if
            !
            if (i2 <= nbsp) then
               !
               call calrhovan(rhovanaux, bec, i2)
               rhovan(:, :, 1) = rhovanaux(:, :, ispin(i2)) ! 0.5 to divide the factor f=2 which accounts for spin multiplicity inside calrhovan
               !
               call rhov(irb, eigrb, rhovan, orb_rhog(:, 2), orb_rhor(:, 2), lgam) !JUST-FOR-NOW ... do we need a factor of 0.5?
               !
            end if
            !
         end if
         !
      end if
      !
!       if(okvan) write(131,*) "rhovan-calrhovan", rhovan(:,:,1) !added:giovanni:debug
!       stop
      deallocate (psi1)
      !
      if (allocated(psi2)) then
         !
         deallocate (psi2)
         !
      end if

      deallocate (orb_rhog)
      !
      if (okvan) then
         !
         deallocate (rhovan)
         deallocate (rhovanaux)
         !
      end if
      !
      call stop_clock('nksic_orbrho')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_get_orbitalrho_twin
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_get_rhoref(i, nnrx, ispin, nspin, f, &
                               rhor, orb_rhor, &
                               rhoref_, rhobar_, rhobarg, grhobar_)
!-----------------------------------------------------------------------
!
! Computes rhoref and rhobar
!
      use kinds, only: dp
      use gvecp, only: ngm
      use funct, only: dft_is_gradient
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use recvecs_indexes, only: np, nm
      use nksic, only: fref, rhobarfact
      use control_flags, only: gamma_only, do_wf_cmplx
      !
      implicit none

      !
      ! input/output vars
      !
      integer, intent(in)  :: i, nnrx
      integer, intent(in)  :: ispin, nspin
      real(dp), intent(in)  :: f
      real(dp), intent(in)  :: rhor(nnrx, nspin)
      real(dp), intent(in)  :: orb_rhor(nnrx)
      real(dp), intent(out) :: rhoref_(nnrx, 2)
      real(dp), intent(out) :: rhobar_(nnrx, 2)
      complex(dp)                :: rhobarg(ngm, 2)
      real(dp), intent(out) :: grhobar_(nnrx, 3, 2)
      !
      integer      :: ig
      complex(dp)  :: fp, fm
      complex(dp), allocatable :: psi(:)
      logical :: lgam

      !
      ! main body
      !
      call start_clock('nksic_get_rhoref')

      lgam = gamma_only .and. .not. do_wf_cmplx
      !write(6,*) ubound(rhobarg)
      !write(6,*) ubound(grhobar_)

      !
      ! define rhobar_i = rho - f_i * rho_i
      !
      if (nspin == 1) then
         rhobar_(:, 1) = rhor(:, 1)*0.5_dp
         rhobar_(:, 2) = rhor(:, 1)*0.5_dp
      else
         rhobar_(:, 1:2) = rhor(:, 1:2)
      end if
      !
      rhobar_(:, ispin) = rhobar_(:, ispin) - f*orb_rhor(:)
      !
      ! probably obsolete
      if (rhobarfact < 1.0d0) then
         rhobar_ = rhobar_*rhobarfact
      end if

      !
      ! define rhoref = rho + (f_ref -f_i) rho_i = rhobar_i + f_ref * rho_i
      ! build rhoref from scratch
      !
      rhoref_(:, 1:2) = rhobar_(:, 1:2)
      rhoref_(:, ispin) = rhoref_(:, ispin) + fref*orb_rhor(:)
      !

      !
      ! compute the gradient of rhobar if needed
      !
      if (dft_is_gradient()) then
         !
         ! allocate( rhobarg(ngm,2) ) modified:giovanni rhobarg became an argument of the subroutine
         allocate (psi(nnrx))
         !
         psi(:) = CMPLX(rhobar_(:, 1), rhobar_(:, 2))
         !
         call fwfft('Dense', psi, dfftp)
         !
         do ig = 1, ngm
            fp = psi(np(ig)) + psi(nm(ig))
            fm = psi(np(ig)) - psi(nm(ig))
            !
            rhobarg(ig, 1) = 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
            rhobarg(ig, 2) = 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
         end do
         !
         call fillgrad(2, rhobarg, grhobar_, lgam)
         !
         deallocate (psi)
         !
      end if
      !
      call stop_clock('nksic_get_rhoref')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_get_rhoref
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_newd(i, nnrx, ispin, nspin, vsic_realspace, vsic_reciprocal, nat, nhm, &
                         becsum, fion, deeq_sic)
!-----------------------------------------------------------------------
!
! computes the deeq coefficients (contributions to the D coeff of USPP)
! for the given orbital i. Coefficients are sotred in deeq_sic
!
      use kinds, only: dp
      use uspp, only: okvan, deeq
      use gvecp, only: ngm
      use cp_main_variables, only: irb, eigrb
      !
      implicit none

      !
      ! input/output vars
      !
      integer, intent(in)    :: i, nnrx, nat, nhm
      integer, intent(in)    :: ispin, nspin
      real(dp), intent(in)    :: vsic_realspace(nnrx)
      complex(dp), intent(in) :: vsic_reciprocal(ngm)
      real(dp), intent(in)    :: becsum(nhm*(nhm + 1)/2, nat, nspin)
      real(dp), intent(inout) :: fion(3, nat)
      real(dp), intent(out)   :: deeq_sic(nhm, nhm, nat)
      !
      ! local vars
      !
      real(dp), allocatable   :: vsic_aux(:, :)

      !
      ! main body
      !
      if (.not. okvan) then
         deeq_sic(:, :, :) = 0.0d0
         return
      end if
      !
      call start_clock('nk_newd')
      !
      allocate (vsic_aux(nnrx, nspin))

      !
      ! fion are updated
      ! deeq coefficients are overwritten
      !
      vsic_aux = 0.0d0
      vsic_aux(:, ispin) = vsic_realspace(:)
      !
      call newd(vsic_aux, irb, eigrb, becsum, fion)
      !
      deeq_sic(:, :, :) = deeq(:, :, :, ispin)

      deallocate (vsic_aux)
      !
      call stop_clock('nk_newd')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_newd
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_get_pz_factor(nspin, ispin, orb_rhor, rhor, &
                                  taukin, tauw, alpha, ibnd, kfact)
!---------------------------------------------------------------
!
! ... sum up the kinetic energy-density taukin ... this works both for summing
!     the orbital-resolved kinetic energy densities, and for the Weizsacker kinetic
!     energy density (involving the total density).
!
      use kinds, only: dp
      use cell_base, only: omega
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use funct, only: dft_is_gradient
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use control_flags, only: gamma_only, do_wf_cmplx
      use nksic, only: epsi2 => epsi2_cutoff_renorm
      !
      implicit none
      !
      integer, intent(in)  :: ispin, ibnd, nspin
      real(dp), intent(in)  :: orb_rhor(nnrx), taukin(nnrx, nspin), tauw(nnrx, nspin), rhor(nnrx, nspin)
      real(dp), intent(out)    :: alpha
      real(dp), intent(in)    :: kfact
      !
      INTEGER :: ir
      LOGICAL :: lgam
      real(dp) :: fact, temp, aidfract, norm, aidspin
!       real(dp), parameter :: epsi=1.d-3
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      temp = 0.d0
      norm = 0.d0
      !
      IF (nspin == 1) THEN
         aidspin = 0.5d0
      ELSE
         aidspin = 1.d0
      END IF
      !
      do ir = 1, nnrx
         !
! ! !          IF((tauw(ir,ispin)**2+taukin(ir,ispin)**2.gt.epsi2**4)) THEN !
         IF (aidspin*rhor(ir, ispin) .gt. epsi2) THEN !
            !
! ! !             aidfract=((tauw(ir,ispin)+epsi2)/(taukin(ir,ispin)+epsi2))**kfact
            aidfract = ((orb_rhor(ir) + epsi2)/(aidspin*rhor(ir, ispin) + epsi2))**kfact
            !
            IF (1.d0 - abs(aidfract) .lt. epsi2) THEN
               !
               aidfract = 1.d0
               !
            END IF
            !
            temp = temp + orb_rhor(ir)*aidfract
            !
         ELSE
            !
            temp = temp + orb_rhor(ir)
            !
         END IF
         !
!          norm=norm+orb_rhor(ir)
         !
      end do
      !
      call mp_sum(temp, intra_image_comm)
!       call mp_sum(norm,intra_image_comm)
      !
      temp = temp*fact
!       norm=norm*fact
!       write(6,*) "checknorm", norm
      !
      alpha = temp
      !
   end subroutine nksic_get_pz_factor

!---------------------------------------------------------------
   subroutine nksic_get_pzfactor_potential(f, nspin, ispin, rhor, orb_rhor, &
                                         pink, taukin, tauw, edens, upsilonkin, upsilonw, vsic_realspace, vsic_reciprocal, alpha, ibnd, kfact)
!---------------------------------------------------------------
!
! ... sum up the kinetic energy-density taukin ... this works both for summing
!     the orbital-resolved kinetic energy densities, and for the Weizsacker kinetic
!     energy density (involving the total density).
!
      use kinds, only: dp
      use cell_base, only: omega
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use gvecp, only: ngm
      use funct, only: dft_is_gradient
      use fft_base, only: dfftp
      use mp, only: mp_sum
      use control_flags, only: gamma_only, do_wf_cmplx
      use nksic, only: epsi2 => epsi2_cutoff_renorm

      !
      implicit none
      !
      integer, intent(in)  :: ispin, ibnd, nspin
      real(dp), intent(in)  :: kfact, f, orb_rhor(nnrx), taukin(nnrx, nspin)
      real(dp), intent(in)  :: tauw(nnrx, nspin), edens(nnrx, nspin), rhor(nnrx, nspin)
      real(dp), intent(inout) :: upsilonkin(nnrx, 3, nspin), upsilonw(nnrx, 3, nspin)
      real(dp), intent(in)    :: alpha
      real(dp), intent(inout) :: pink
      real(dp), intent(out)   :: vsic_realspace(nnrx)
      complex(dp), intent(out):: vsic_reciprocal(ngm)
      !
      INTEGER :: ir, j
      LOGICAL :: lgam
      real(dp) :: fact, temp, tempw, aidtau, aidfrac, aidspin
      complex(dp), allocatable :: rhog_dummy(:, :)
      real(dp), allocatable :: upsilonh(:, :, :), vsicaux(:, :)
      complex(dp), allocatable :: vtmp(:)
!       real(dp), parameter :: epsi=1.d-3
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      allocate (upsilonh(nnrx, 3, nspin))
      allocate (vsicaux(nnrx, nspin))
      allocate (rhog_dummy(1, 1))
      !
      upsilonh = 0.d0
      !
      vsicaux = 0.d0
      !
      IF (nspin == 1) THEN
         aidspin = 0.5d0
      ELSE
         aidspin = 1.d0
      END IF
!       write(6,*) "checkall", ibnd, ispin
      !
! ! !       call nksic_get_upsilon_pz( f, nspin, ispin, orb_rhor, &
! ! !                                       upsilonkin, ibnd)
      if (ibnd == 1) THEN !compute also upsilonw
         !
! !          call nksic_get_upsilon_pz( 1.d0, nspin, ispin, rhor(:,ispin), &
!                                       upsilonw, ibnd)
         !
      END IF
      !
      !
      upsilonh = 0.d0
      !
!       vsicaux(:,ispin)=vsicaux(:,ispin)
      !
      do ir = 1, nnrx
         !
         temp = 0.d0
         tempw = 0.d0
         !
         do j = 1, 3
            !
            temp = temp + upsilonkin(ir, j, ispin)**2.
            tempw = tempw + upsilonw(ir, j, ispin)**2.
            !
         end do
         !
         !
!          temp=sqrt(abs(temp))
!          tempw=sqrt(abs(tempw))
!                              write(6,*) "checktau", taukin(ir,ispin),tauw(ir,ispin),ispin
         !
! ! !          IF((tauw(ir,ispin)**2+taukin(ir,ispin)**2.gt.epsi2**4)) THEN
         IF ((aidspin*rhor(ir, ispin) .gt. epsi2)) THEN ! ! THEN
            !
! ! !             aidtau=0.5d0*(temp/(taukin(ir,ispin)+epsi2)-tempw/(tauw(ir,ispin)+epsi2))
            aidtau = -edens(ir, ispin)/(aidspin*rhor(ir, ispin) + epsi2)**(kfact + 1.d0)
            !
! ! !             aidfrac=((tauw(ir,ispin)+epsi2)/(taukin(ir,ispin)+epsi2))**kfact
            aidfrac = ((orb_rhor(ir) + epsi2)/(aidspin*rhor(ir, ispin) + epsi2))**kfact

            !
            IF (1.d0 - abs(aidfrac) .lt. epsi2) THEN
!                !
               aidfrac = 1.d0
               aidtau = 0.d0
               !
            END IF
! !             !
            IF (abs(aidtau) .lt. epsi2) THEN
               !
               aidtau = 0.d0
               !
            END IF
            !
            vsicaux(ir, ispin) = vsicaux(ir, ispin) &
                                 + pink/f*(-alpha + aidfrac)
            !
! ! !             vsicaux(ir,ispin) = vsicaux(ir,ispin)+kfact*edens(ir,ispin)*aidfrac*aidtau
            vsicaux(ir, ispin) = vsicaux(ir, ispin) + kfact*aidfrac*pink/f + aidtau*kfact
            !
            do j = 1, 3
               !
               !aidtau=0.5d0*(upsilonkin(ir,j,ispin)/(taukin(ir,ispin)+epsi2)-upsilonw(ir,j,ispin)/(tauw(ir,ispin)+epsi2))
               !
               IF (abs(aidfrac - 1.d0) .lt. epsi2) THEN !abs(aidtau).lt.epsi2**2
                  !
                  !aidtau=0.d0
                  !
               END IF
               !
               !upsilonh(ir,j,ispin) = upsilonh(ir,j,ispin) - kfact*edens(ir,ispin)*aidfrac*aidtau
               !
            end do
            !
         ELSE IF (abs(1.d0 - alpha) .gt. epsi2**2) THEN
            !
            vsicaux(ir, ispin) = vsicaux(ir, ispin) &
                                 + pink/f*(-alpha + 1.d0)
            !
         END IF
         !
      end do
      !
      ! Now we need to use fft's to add the gradient part to vsic... check the sign of this expression
      !
!       call gradh( 1, upsilonh(:,:,ispin:ispin), rhog_dummy, vsicaux(:,ispin:ispin), dexc_dummy, lgam )
      !
      do ir = 1, nnrx
         !
         vsic_realspace(ir) = vsic_realspace(ir) + vsicaux(ir, ispin)
         !
      end do
      !
      ! Store the vsic_realspace potential in reciprocal space
      allocate (vtmp(ngm))
      vtmp = vsic_realspace(:)
      call fwfft('Dense', vtmp, dfftp)
      call psi2rho('Dense', vtmp, dfftp%nnr, vsic_reciprocal, ngm)
      deallocate (vtmp)
      !
      deallocate (upsilonh, vsicaux, rhog_dummy)
      !
   end subroutine nksic_get_pzfactor_potential

!---------------------------------------------------------------
   subroutine add_up_taukin(nnrx, taukin, grhoraux, orb_rhor, f)
!---------------------------------------------------------------
      !
      USE kinds, only: DP
      use nksic, only: epsi => epsi_cutoff_renorm

      !
      INTEGER, INTENT(IN) :: nnrx
      REAL(DP) :: taukin(nnrx), orb_rhor(nnrx), f, grhoraux(nnrx, 3)
      !
      REAL(DP) :: temp_gradient, temp_rho
!         REAL(DP), PARAMETER :: epsi2=1.e-11
      INTEGER :: ir
      !

      do ir = 1, nnrx
         !
         temp_gradient = grhoraux(ir, 1)**2 + grhoraux(ir, 2)**2 + grhoraux(ir, 3)**2
         temp_rho = orb_rhor(ir)

         IF ((temp_gradient .lt. epsi**2)) THEN!(temp_rho.lt.epsi.or.temp_gradient.lt.epsi**2) THEN
            temp_gradient = 0.d0
            temp_rho = 1.d0
         ELSE
            taukin(ir) = taukin(ir) + f/(2.)*temp_gradient
         END IF
         !
      end do

   end subroutine add_up_taukin

!---------------------------------------------------------------
   subroutine nksic_get_taukin_pz(f, nspin, ispin, orb_rhor, &
                                  taukin, ibnd, mult)
!---------------------------------------------------------------
!
! ... sum up the kinetic energy-density taukin ... this works both for summing
!     the orbital-resolved kinetic energy densities, and for the Weizsacker kinetic
!     energy density (involving the total density).
!
      use kinds, only: dp
!       use nksic,                only : add_up_taukin
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use recvecs_indexes, only: np
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dft_is_gradient
      use control_flags, only: gamma_only, do_wf_cmplx
      use nksic, only: epsi => epsi_cutoff_renorm

      !
      implicit none
      !
      integer, intent(in)  :: ispin, ibnd, nspin, mult
      real(dp), intent(in)  :: f, orb_rhor(nnrx)
      real(dp), intent(inout) :: taukin(nnrx, nspin)
      !
      INTEGER :: ig, ir
      complex(dp), allocatable :: rhogaux(:, :)
      real(dp), allocatable :: grhoraux(:, :, :)
      complex(dp), allocatable :: vhaux(:)
      LOGICAL :: lgam
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      !
      allocate (rhogaux(ngm, 2))
      allocate (vhaux(nnrx))
      !
      IF (ibnd == 1) THEN !first band: initialize taukin for this spin_loop
         !
         taukin(1:nnrx, ispin) = 0.d0
         !
      END IF
      !
      rhogaux = 0.0_dp
      !
      do ir = 1, nnrx
         !
         vhaux(ir) = sqrt(abs(orb_rhor(ir) + mult*epsi))
         !
      end do
      !
      call fwfft('Dense', vhaux, dfftp)
      !
      do ig = 1, ngm
         rhogaux(ig, ispin) = vhaux(np(ig))
      end do
      !
!       call enkin_dens( rhogaux(:,ispin), ngm, f)
      !
      allocate (grhoraux(nnrx, 3, 2))

      grhoraux = 0.0_dp

      call fillgrad(1, rhogaux(:, ispin:ispin), grhoraux(:, :, ispin:ispin), lgam)

      call add_up_taukin(nnrx, taukin(:, ispin), grhoraux(:, :, ispin), orb_rhor(:), f)
!       vhaux=0.d0
!       do ig=1,ngm
!           vhaux( np(ig) )= rhogaux(ig,ispin)
!           vhaux( nm(ig) )= CONJG(rhogaux(ig,ispin))
!       enddo
!       !
!       call invfft('Dense',vhaux,dfftp )
!       !
!       do ir=1,nnrx
!          !
!           taukin(ir,ispin) = DBLE(vhaux(ir))
!          !
!       enddo
      !
      deallocate (vhaux, rhogaux, grhoraux)
      !
   end subroutine nksic_get_taukin_pz

!---------------------------------------------------------------
   subroutine nksic_get_upsilon_pz(f, nspin, ispin, orb_rhor, &
                                   upsilon, ibnd)
!---------------------------------------------------------------
!
! ... sum up the kinetic energy-density taukin ... this works both for summing
!     the orbital-resolved kinetic energy densities, and for the Weizsacker kinetic
!     energy density (involving the total density).
!
      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use recvecs_indexes, only: np
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dft_is_gradient
      use control_flags, only: gamma_only, do_wf_cmplx
      use nksic, only: epsi => epsi_cutoff_renorm
      !
      implicit none
      !
      integer, intent(in)  :: ispin, ibnd, nspin
      real(dp), intent(in)  :: f, orb_rhor(nnrx)
      real(dp), intent(out) :: upsilon(nnrx, 3, nspin)
      !
      INTEGER :: ig, ir, j
      complex(dp), allocatable :: rhogaux(:, :)
      real(dp), allocatable :: grhoraux(:, :, :)
      complex(dp), allocatable :: vhaux(:)
      real(dp) :: temp(3), tempnorm
!       real(dp), parameter :: epsi=1.d-3
      LOGICAL :: lgam
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      !
      allocate (rhogaux(ngm, 2))
      allocate (vhaux(nnrx))
      !
      rhogaux = 0.0_dp
      !
      do ir = 1, nnrx
         !
         vhaux(ir) = log(abs(orb_rhor(ir)))
         !
      end do
      !
      call fwfft('Dense', vhaux, dfftp)
      !
      do ig = 1, ngm
         rhogaux(ig, ispin) = vhaux(np(ig))
      end do
      !
      allocate (grhoraux(nnrx, 3, 2))
      !
      grhoraux = 0.0_dp
      !
      call fillgrad(1, rhogaux(:, ispin:ispin), grhoraux(:, :, ispin:ispin), lgam)
      !
      upsilon(1:nnrx, 1:3, ispin) = 0.d0
      !
      do ir = 1, nnrx
         !
         IF (.true.) THEN
            !
            tempnorm = 0.d0
            !
            do j = 1, 3
               !
               temp(j) = grhoraux(ir, j, ispin)!/(2.*(orb_rhor(ir)+epsi))
               tempnorm = tempnorm + temp(j)**2
               !
            end do
            !
            IF (tempnorm .gt. epsi) THEN
               !
               upsilon(ir, :, ispin) = temp(:)
               !
            END IF
            !
         END IF
         !
      end do
      !
      deallocate (vhaux, rhogaux, grhoraux)
      !
   end subroutine nksic_get_upsilon_pz

!-----------------------------------------------------------------------
   subroutine nksic_eforce(i, nbsp, nx, vsic_realspace, vsic_reciprocal, deeq_sic, bec, ngw, c1, c2, vsicpsi, lgam)
!-----------------------------------------------------------------------
!
! Compute vsic_realspace potential for orbitals i and i+1 (c1 and c2)
!
      use kinds, only: dp
      use cp_interfaces, only: fwfft, invfft
      use fft_base, only: dffts, dfftp
      use gvecp, only: ngm
      use gvecs, only: ngs, nps, nms
      use grid_dimensions, only: nnrx
      use smooth_grid_dimensions, only: nnrsx
      use uspp, only: nkb, vkb
      use uspp_param, only: nhm, nh
      use cvan, only: ish
      use ions_base, only: nsp, na, nat
      use twin_types
      !
      implicit none

      !
      ! input/output vars
      !
      integer, intent(in)  :: i, nbsp, nx, ngw
      real(dp), intent(in)  :: vsic_realspace(nnrx, nx)
      complex(dp), intent(in)  :: vsic_reciprocal(ngm, nx)
      real(dp), intent(in)  :: deeq_sic(nhm, nhm, nat, nx)
      type(twin_matrix), intent(in)  :: bec!(nkb,nbsp) !modified:giovanni
      complex(dp), intent(in)  :: c1(ngw), c2(ngw)
      complex(dp), intent(out) :: vsicpsi(ngw, 2)
      logical, intent(in) :: lgam !added:giovanni

      !
      ! local vars
      !
      character(12) :: subname = 'nksic_eforce'
      real(dp), allocatable      :: vsic_realspace_local(:)
      integer       :: ir, ig, ierr, j
      integer       :: is, iv, jv, isa, ism
      integer       :: ivoff, jvoff, ia, inl, jnl
      real(dp)      :: wfc(2), dd
      complex(dp) :: wfc_c(2)
      complex(dp)   :: fm, fp
      complex(dp), allocatable :: psi(:), psi1(:), psi2(:)
      real(dp), allocatable :: aa(:, :)
      complex(dp), allocatable :: aa_c(:, :)
      complex(dp), parameter :: c_one = CMPLX(1.d0, 0.d0)

      !
      !====================
      ! main body
      !====================
      !
      call start_clock('nk_eforce')
      !
      allocate (vsic_realspace_local(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating vsic_realspace_local', abs(ierr))
      allocate (psi(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi', abs(ierr))
      allocate (psi1(nnrx), stat=ierr)
      if (ierr /= 0) call errore(subname, 'allocating psi1', abs(ierr))
      if (.not. lgam) then
         allocate (psi2(nnrx), stat=ierr)
         if (ierr /= 0) call errore(subname, 'allocating psi2', abs(ierr))
      end if

      !
      ! init
      !
      vsicpsi(:, :) = 0.0d0
      !
      ! take advantage of the smooth and the dense grids
      ! being equal (NCPP case)
      !
      if (nnrsx == nnrx) then !waring:giovanni we are not using ultrasoft
         !
         ! no need to take care of the double grid.
         ! typically, NCPP case

         !
         if (lgam) then
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 2) !warning:giovanni need to change this
         else
            CALL c2psi(psi1, nnrx, c1, c2, ngw, 0) !warning:giovanni need to change this
            CALL c2psi(psi2, nnrx, c2, c1, ngw, 0) !warning:giovanni need to change this
         end if
         !
         CALL invfft('Dense', psi1, dfftp)
         if (.not. lgam) then
            CALL invfft('Dense', psi2, dfftp)
         end if
         !
         ! Transform the ODD potential from reciprocal to real space
         !
         call rho2psi('Dense', psi, dfftp%nnr, vsic_reciprocal(:, i), ngm)
         call invfft('Dense', psi, dfftp)
         vsic_realspace_local = dble(psi)
         !
         call ebl_check(vsic_realspace(:, i), vsic_reciprocal(:, i))
         !
         ! computing the orbital wfcs
         ! and the potentials in real space on the full grid
         !
         if (lgam) then
            do ir = 1, nnrx
               !
               wfc(1) = DBLE(psi1(ir))
               wfc(2) = AIMAG(psi1(ir))
               !
               psi1(ir) = CMPLX(wfc(1)*vsic_realspace(ir, i), wfc(2)*vsic_realspace(ir, i + 1))
               !
            end do
         else
            do ir = 1, nnrx
               !
               wfc_c(1) = psi1(ir)
               wfc_c(2) = psi2(ir)
               !
               psi1(ir) = wfc_c(1)*vsic_realspace(ir, i)
               psi2(ir) = wfc_c(2)*vsic_realspace(ir, i + 1)
               !
            end do
         end if
         !

         CALL fwfft('Dense', psi1, dfftp)
         if (.not. lgam) then
            CALL fwfft('Dense', psi2, dfftp)
         end if
         !
         vsicpsi(:, :) = 0.0_dp
         !
         if (lgam) then
            do ig = 1, ngw
               !
               fp = psi1(nps(ig)) + psi1(nms(ig))
               fm = psi1(nps(ig)) - psi1(nms(ig))
               !
               vsicpsi(ig, 1) = 0.5d0*CMPLX(DBLE(fp), AIMAG(fm))
               vsicpsi(ig, 2) = 0.5d0*CMPLX(AIMAG(fp), -DBLE(fm))
               !
            end do
         else
            do ig = 1, ngw
               !
               fp = psi1(nps(ig))
               fm = psi2(nps(ig))
               !
               vsicpsi(ig, 1) = fp
               vsicpsi(ig, 2) = fm
               !
            end do
         end if

      else
         write (6, *) "WARNING, WE ARE USING USPP"
         !
         ! here we take properly into account the
         ! smooth and the dense grids
         ! typically, USPP case
         !
         CALL nksic_eforce_std(lgam) !warning:giovanni this makes fourier transforms
         !
      end if
      !
      deallocate (psi)
      deallocate (psi1)
      if (.not. lgam) then
         deallocate (psi2)
      end if

      !
      ! add USPP non-local contribution
      ! to the potantial
      ! (this comes from the orbital-dependent piece of
      ! the potential)
      !
      if (nkb > 0) then
!               write(6,*) "WE ARE USING USPP --- WARNING"
         !
         !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
         !
         if (.not. bec%iscmplx) then
            allocate (aa(nkb, 2))
            !
            aa = 0.0d0
            !
            !
            do is = 1, nsp
               !
               do iv = 1, nh(is)
               do jv = 1, nh(is)
                  !
                  isa = 0
                  do ism = 1, is - 1
                     isa = isa + na(ism)
                  end do
                  !
                  ivoff = ish(is) + (iv - 1)*na(is)
                  jvoff = ish(is) + (jv - 1)*na(is)
                  !
                  if (i /= nbsp) then
                     !
                     do ia = 1, na(is)
                        inl = ivoff + ia
                        jnl = jvoff + ia
                        isa = isa + 1
                        !
                        dd = deeq_sic(iv, jv, isa, i)
                        aa(inl, 1) = aa(inl, 1) + dd*bec%rvec(jnl, i)
                        !
                        dd = deeq_sic(iv, jv, isa, i + 1)
                        aa(inl, 2) = aa(inl, 2) + dd*bec%rvec(jnl, i + 1)
                        !
                     end do
                     !
                  else
                     !
                     do ia = 1, na(is)
                        inl = ivoff + ia
                        jnl = jvoff + ia
                        isa = isa + 1
                        !
                        dd = deeq_sic(iv, jv, isa, i)
                        aa(inl, 1) = aa(inl, 1) + dd*bec%rvec(jnl, i)
                        !
                     end do
                     !
                  end if
                  !
               end do
               end do
               !
            end do
            !
!               write(6,*) "deeq_sic"
!               write(6,*) deeq_sic
            !
            call DGEMM('N', 'N', 2*ngw, 2, nkb, 1.0d0, &
                       vkb, 2*ngw, aa, nkb, 1.0d0, vsicpsi(:, :), 2*ngw)
            !
            deallocate (aa)
            !
         else
            allocate (aa_c(nkb, 2))
            !
            aa_c = CMPLX(0.0d0, 0.d0)
            !
            !
            do is = 1, nsp
               !
               do iv = 1, nh(is)
               do jv = 1, nh(is)
                  !
                  isa = 0
                  do ism = 1, is - 1
                     isa = isa + na(ism)
                  end do
                  !
                  ivoff = ish(is) + (iv - 1)*na(is)
                  jvoff = ish(is) + (jv - 1)*na(is)
                  !
                  if (i /= nbsp) then
                     !
                     do ia = 1, na(is)
                        inl = ivoff + ia
                        jnl = jvoff + ia
                        isa = isa + 1
                        !
                        dd = deeq_sic(iv, jv, isa, i)
                        aa_c(inl, 1) = aa_c(inl, 1) + dd*bec%cvec(jnl, i) !warning:giovanni or conjg?
                        !
                        dd = deeq_sic(iv, jv, isa, i + 1)
                        aa_c(inl, 2) = aa_c(inl, 2) + dd*bec%cvec(jnl, i + 1) !warning:giovanni or conjg?
                        !
                     end do
                     !
                  else
                     !
                     do ia = 1, na(is)
                        inl = ivoff + ia
                        jnl = jvoff + ia
                        isa = isa + 1
                        !
                        dd = deeq_sic(iv, jv, isa, i)
                        aa_c(inl, 1) = aa_c(inl, 1) + dd*bec%cvec(jnl, i) !warning:giovanni or conjg?
                        !
                     end do
                     !
                  end if
                  !
               end do
               end do
               !
            end do
            !
            call ZGEMM('N', 'N', ngw, 2, nkb, c_one, &
                       vkb, ngw, aa_c, nkb, c_one, vsicpsi(:, :), ngw)
            !
            deallocate (aa_c)
            !
         end if
      end if

      call stop_clock('nk_eforce')
      return

!
! implementation to deal with both
! the smooth and the dense grids
!
   CONTAINS
      !
      subroutine nksic_eforce_std(lgam)
         !
         use smooth_grid_dimensions, only: nnrsx
         use recvecs_indexes, only: np
         implicit none

         logical, intent(IN) :: lgam
         !
         complex(dp) :: c(ngw, 2)
         complex(dp) :: psis(nnrsx)
         complex(dp) :: vsicg(nnrx)
         complex(dp) :: vsics(nnrsx)
         complex(dp) :: vsicpsis(nnrsx)

         c(:, 1) = c1
         c(:, 2) = c2

         do j = 1, 2
            !
            psis = 0.d0
            if (lgam) then
               do ig = 1, ngw
                  psis(nms(ig)) = CONJG(c(ig, j))
                  psis(nps(ig)) = c(ig, j)
               end do
            else
               do ig = 1, ngw
                  psis(nps(ig)) = c(ig, j)
               end do
            end if
            call invfft('Wave', psis, dffts)
            !
            vsicg(1:nnrx) = vsic_realspace(1:nnrx, i + j - 1)
            call fwfft('Dense', vsicg, dfftp)
            !
            vsics = 0.0_dp
            if (lgam) then
               do ig = 1, ngs
                  vsics(nps(ig)) = vsicg(np(ig))
                  vsics(nms(ig)) = CONJG(vsicg(np(ig)))
               end do
            else
               do ig = 1, ngs
                  vsics(nps(ig)) = vsicg(np(ig))
                  vsics(nms(ig)) = CONJG(vsicg(np(ig)))
               end do
            end if
            !
            call invfft('Smooth', vsics, dffts)
            !
            vsicpsis = 0.0_dp
            if (lgam) then
               do ir = 1, nnrsx
                  vsicpsis(ir) = CMPLX(DBLE(vsics(ir))*DBLE(psis(ir)), 0.0_dp)
               end do
            else
               do ir = 1, nnrsx
                  vsicpsis(ir) = CMPLX(DBLE(vsics(ir))*DBLE(psis(ir)), DBLE(vsics(ir))*AIMAG(psis(ir)))
               end do
            end if
            !
            call fwfft('Wave', vsicpsis, dffts)
            !
            do ig = 1, ngw
               vsicpsi(ig, j) = vsicpsis(nps(ig))
            end do
            !
         end do
         !
      end subroutine nksic_eforce_std
      !
!---------------------------------------------------------------
   end subroutine nksic_eforce
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_dmxc_spin_cp(nnrx, rhoref, f, ispin, rhoele, &
                                 small, wref, wxd)
!---------------------------------------------------------------
!
! the derivative of the xc potential with respect to the local density
! is computed.
! In order to save time, the loop over space coordinates is performed
! inside this routine (inlining).
!
! NOTE: wref and wsic are UPDATED and NOT OVERWRITTEN by this subroutine
!
      USE kinds, ONLY: dp
      USE funct, ONLY: xc_spin, get_iexch, get_icorr
      implicit none
      !
      integer, intent(in)    :: nnrx, ispin
      real(dp), intent(in)    :: rhoref(nnrx, 2), rhoele(nnrx, 2)
      real(dp), intent(in)    :: f, small
      real(dp), intent(inout) :: wref(nnrx), wxd(nnrx, 2)
      !
      character(18) :: subname = 'nksic_dmxc_spin_cp'
      real(dp) :: rhoup, rhodw, rhotot, zeta
      real(dp) :: dmuxc(2, 2)
      real(dp) :: rs, ex, vx, dr, dz, ec, &
                  vcupm, vcdwm, vcupp, vcdwp, dzm, dzp, fact
      !real(dp) :: vxupp, vxdwp, vxupm, vxdwm
      real(dp), external :: dpz, dpz_polarized
      integer :: ir
      !logical :: do_exch, do_corr
      !
      real(dp), parameter :: e2 = 2.0_dp, &
                             pi34 = 0.6203504908994_DP, & ! redefined to pi34=(3/4pi)^(1/3)
                             pi34_old = 0.75_dp/3.141592653589793_dp, third = 1.0_dp/3.0_dp, &
                             p43 = 4.0_dp/3.0_dp, p49 = 4.0_dp/9.0_dp, m23 = -2.0_dp/3.0_dp

      !
      ! mian body
      !
      !CALL start_clock( 'nk_dmxc_spin_cp' )
      !
      ! the current implementation works only on top
      ! of LSD and LDA. Other functionals have to
      ! be implemented explicitly. To do that, we need to
      ! call the proper xc-routine (at the moment we call
      ! slater and pz_corr)
      !
      if (get_iexch() /= 1 .or. get_icorr() /= 1) &
         call errore(subname, 'only LDA/LSD PZ functionals implemented', 10)
      !
      !do_exch = ( get_iexch() == 1 )
      !do_corr = ( get_icorr() == 1 )
      !
      !
      ! main loop
      !
      do ir = 1, nnrx
         !
         dmuxc(:, :) = 0.0_dp
         !
         rhoup = rhoref(ir, 1)
         rhodw = rhoref(ir, 2)
         rhotot = rhoup + rhodw
         !
         if (rhotot < small) cycle
         !
         zeta = (rhoup - rhodw)/rhotot
         if (abs(zeta) > 1.0_dp) zeta = sign(1.0_dp, zeta)

         !
         ! calculate exchange contribution (analytical)
         !
         if (rhoup > small) then
            rs = pi34/(2.0_dp*rhoup)**third
            call slater(rs, ex, vx)
            dmuxc(1, 1) = vx/(3.0_dp*rhoup)
         end if
         !
         if (rhodw > small) then
            rs = pi34/(2.0_dp*rhodw)**third
            call slater(rs, ex, vx)
            dmuxc(2, 2) = vx/(3.0_dp*rhodw)
         end if

         !
         ! calculate correlation contribution (numerical)
         !
         dr = min(1.e-6_dp, 1.e-4_dp*rhotot)
         fact = 0.5d0/dr
         !
         ! the explicit call to the correlation part only
         ! are performed instead of calling xc_spin.
         ! this saves some CPU time.
         ! unfortunately, different functionals have then
         ! to be treated explicitly
         !
         !call xc_spin(rhotot-dr,zeta,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
         !call xc_spin(rhotot+dr,zeta,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
         !
         rs = pi34/(rhotot - dr)**third
         call pz_spin(rs, zeta, ec, vcupm, vcdwm)
         rs = pi34/(rhotot + dr)**third
         call pz_spin(rs, zeta, ec, vcupp, vcdwp)
         !
         dmuxc(1, 1) = dmuxc(1, 1) + (vcupp - vcupm)*fact
         dmuxc(1, 2) = dmuxc(1, 2) + (vcupp - vcupm)*fact
         dmuxc(2, 1) = dmuxc(2, 1) + (vcdwp - vcdwm)*fact
         dmuxc(2, 2) = dmuxc(2, 2) + (vcdwp - vcdwm)*fact

         dz = 1.e-6_dp
         dzp = min(1.0, zeta + dz) - zeta
         dzm = -max(-1.0, zeta - dz) + zeta
         !
         fact = 1.0d0/(rhotot*(dzp + dzm))
         !
         !call xc_spin(rhotot,zeta-dzm,ex,ec,vxupm,vxdwm,vcupm,vcdwm)
         !call xc_spin(rhotot,zeta+dzp,ex,ec,vxupp,vxdwp,vcupp,vcdwp)
         !
         rs = pi34/(rhotot)**third
         call pz_spin(rs, zeta - dzm, ec, vcupm, vcdwm)
         call pz_spin(rs, zeta + dzp, ec, vcupp, vcdwp)

         dmuxc(1, 1) = dmuxc(1, 1) + (vcupp - vcupm)*(1.0_dp - zeta)*fact
         dmuxc(1, 2) = dmuxc(1, 2) - (vcupp - vcupm)*(1.0_dp + zeta)*fact
         dmuxc(2, 1) = dmuxc(2, 1) + (vcdwp - vcdwm)*(1.0_dp - zeta)*fact
         dmuxc(2, 2) = dmuxc(2, 2) - (vcdwp - vcdwm)*(1.0_dp + zeta)*fact

         !
         ! add corrections to the nksic potentials
         !
         wxd(ir, 1) = wxd(ir, 1) + dmuxc(1, ispin)*rhoele(ir, ispin)*f
         wxd(ir, 2) = wxd(ir, 2) + dmuxc(2, ispin)*rhoele(ir, ispin)*f
         !
         wref(ir) = wref(ir) + dmuxc(ispin, ispin)*rhoele(ir, ispin)
         !
      end do

      return
      !
!---------------------------------------------------------------
   end subroutine nksic_dmxc_spin_cp
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rot_emin(nouter, ninner, etot, Omattot, lgam)
!-----------------------------------------------------------------------
!
! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use constants, only: PI
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use io_global, only: stdout, ionode
      use electrons_base, only: nbsp, nbspx, nspin, &
                                iupdwn, nupdwn
      use cp_interfaces, only: invfft
      use nksic, only: vsic_realspace, vsic_reciprocal, pink, &
                       do_nk, do_wref, do_wxd, &
                       innerloop_nmax
      use uspp, only: nkb
      use cp_main_variables, only: bec
      use wavefunctions_module, only: c0, cm
      use control_flags, only: esic_conv_thr
      use cg_module, only: tcg
      use twin_types
      !
      implicit none
      !
      ! in/out vars
      !
      integer                  :: ninner
      integer, intent(in)  :: nouter
      real(dp), intent(in)  :: etot
      complex(dp)                 :: Omattot(nbspx, nbspx)
      logical :: lgam

      !
      ! local variables
      !
      real(dp)    :: esic, esic_old
      integer     :: nbnd1, nbnd2
      integer     :: npassofailmax
      real(dp)    :: dtmp, dalpha
      integer     :: isp
      real(dp)    :: vsicah2sum, deigrms, dmaxeig
      logical     :: do_nonvar, lstopinner
      !
      complex(dp), allocatable :: Omat1tot(:, :)
      complex(dp), allocatable :: Umatbig(:, :)
      real(dp), allocatable :: Heigbig(:)
      complex(dp), allocatable :: wfc_ctmp(:, :)
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)
      complex(dp), allocatable :: vsicah(:, :)
      real(dp), allocatable :: vsic1_realspace(:, :)
      complex(dp), allocatable :: vsic_reciprocal1(:, :)
      type(twin_matrix) :: bec1
!       real(dp),    allocatable :: bec1(:,:)
      real(dp), allocatable :: pink1(:)
      !
      integer, save :: npassofail = 0
      real(dp), save :: passoprod = 0.3d0

      !
      ! variables for test calculations - along gradient line direction
      !
      logical :: ldotest

      !
      ! main body
      !
      CALL start_clock('nk_rot_emin')

      !
      npassofailmax = 5 ! when to stop dividing passoprod by 2
      esic_old = 0.d0

      allocate (Omat1tot(nbspx, nbspx))
      allocate (Umatbig(nbspx, nbspx))
      allocate (Heigbig(nbspx))
      allocate (wfc_ctmp(ngw, nbspx))
      allocate (vsic1_realspace(nnrx, nbspx))
      allocate (vsic_reciprocal1(ngm, nbspx))
      allocate (pink1(nbspx))

      call init_twin(bec1, lgam)
      call allocate_twin(bec1, nkb, nbsp, lgam)
!       allocate( bec1(nkb,nbsp) )
      !
      Umatbig(:, :) = (0.d0, 0.d0)
      Heigbig(:) = 0.d0
      deigrms = 0.d0

      Omattot(:, :) = 0.d0
      do nbnd1 = 1, nbspx
         Omattot(nbnd1, nbnd1) = 1.d0
      end do

      ninner = 0
      ldotest = .false.

      !
      ! init IO
      if (ionode) write (stdout, "(14x,'# iter',6x,'etot',17x,'esic',&
                                  & 17x,'deigrms')")

      !
      ! main loop
      !
      inner_loop: &
         do while (.true.)

         call start_clock("nk_innerloop")
         !
         ninner = ninner + 1

         if (ninner > innerloop_nmax) then
            !
#ifdef __DEBUG
            if (ionode) write (1031, *) '# innerloop_nmax reached.'
            if (ionode) write (1031, *)
#endif
            if (ionode) then
               write (stdout, "(14x,'# innerloop_nmax reached.',/)")
            end if
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if

#ifdef __DEBUG
         !
!$$     ! Now do the test
         !
         if (mod(ninner, 10) == 1 .or. ninner <= 5) ldotest = .true.
         !ldotest=.true.
         if (ldotest) then
            !
            dtmp = 4.d0*PI
            !call nksic_rot_test(dtmp,201,nouter,ninner,etot)
            ldotest = .false.
            !
         end if
#endif

         !
         ! This part calculates the anti-hermitian part of the hamiltonian
         ! vsicah and see whether a convergence has been achieved
         !
         wfc_ctmp(:, :) = (0.d0, 0.d0)
         deigrms = 0.d0

         spin_loop: &
            do isp = 1, nspin
            !
            allocate (Umat(nupdwn(isp), nupdwn(isp)))
            allocate (Heig(nupdwn(isp)))
            allocate (vsicah(nupdwn(isp), nupdwn(isp)))
            !
            call nksic_getvsicah_new2(isp, vsicah, vsicah2sum, lgam)
            !
            call nksic_getHeigU(isp, vsicah, Heig, Umat)

            Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                    iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Umat(:, :)
            Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Heig(:)

            !!
            !! CHP: The following file prints out
            !! the eigenvalues of the force matrix for debugging
            !
            !if (ionode) then
            !     nfile=10000+isp
            !     write(nfile,'(2I10,100F10.6)') ninner,nouter,sum(Heig(:)**2),Heig(:)
            !endif
            !
            deigrms = deigrms + sum(Heig(:)**2)

            deallocate (Umat)
            deallocate (Heig)
            deallocate (vsicah)
            !
         end do spin_loop

         dmaxeig = max(abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
         do isp = 2, nspin
            !
            dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp))))
            dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
            !
         end do

         ! how severe the transform is
         deigrms = sqrt(deigrms/nbsp)

         !
         ! print out ESIC part & other total energy
         !
         esic = sum(pink(:))
         !
#ifdef __DEBUG
         if (ionode) write (1031, '(2I10,3F24.13)') ninner, nouter, etot, esic, deigrms
#endif
         if (ionode) write (stdout, '(10x,2i5,3F21.13)') ninner, nouter, etot, esic, deigrms

         dalpha = passoprod/dmaxeig
         !
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp, Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, dtmp, lgam)

         !
         ! deal with non-variational functionals,
         ! such as NK0
         !
         do_nonvar = (do_nk .and. (.not. do_wref .or. .not. do_wxd))
         !
         if (do_nonvar) then
            lstopinner = (ninner >= 2 .and. &
                          ((esic - dtmp)*(esic - esic_old) > 0.d0))
         else
            lstopinner = (dtmp >= esic)
         end if
         !
         lstopinner = (lstopinner .or. (abs(esic - dtmp) < esic_conv_thr))

         if (lstopinner) then
            !
            npassofail = npassofail + 1
            !
#ifdef __DEBUG
            if (ionode) write (1031, '("# procedure  ",I4," / ",I4, &
                                    & " is finished.",/)') npassofail, npassofailmax
#endif
            if (ionode) write (stdout, '(14x, "# procedure  ",I4," / ",I4, &
                                    & " is finished.",/)') npassofail, npassofailmax
            !
            ! if we reach at the maximum allowed npassofail number,
            ! we exit without further update
            !
            if (npassofail >= npassofailmax) then
               !
               ninner = ninner + 1
               call stop_clock("nk_innerloop")
               exit inner_loop
               !
            end if
            !
            passoprod = passoprod*0.5d0
            ! ldotest=.true.
            cycle
            !
         end if

         !
         ! we keep track of all the rotations to rotate cm later
         !
         Omattot = MATMUL(Omattot, Omat1tot)
         !
         pink(:) = pink1(:)
         vsic_realspace(:, :) = vsic1_realspace(:, :)
         vsic_reciprocal(:, :) = vsic_reciprocal1(:, :)
         call copy_twin(bec, bec1)
!         bec%rvec(:,:)   = bec1(:,:)
         c0(:, :) = wfc_ctmp(:, :)
         esic_old = esic

         call stop_clock("nk_innerloop")
         !
      end do inner_loop

      !
      ! Wavefunction cm rotation according to Omattot
      ! cm is relevant only for damped dynamics
      !
      call start_clock("nk_rot_cm")
      if (.not. tcg) then
         !
         if (ninner >= 2) then
            !
            wfc_ctmp(:, :) = (0.d0, 0.d0)
            !
            do nbnd1 = 1, nbspx
            do nbnd2 = 1, nbspx
               wfc_ctmp(:, nbnd1) = wfc_ctmp(:, nbnd1) + cm(:, nbnd2)*Omattot(nbnd2, nbnd1)
               ! XXX (we can think to use a blas, here, and split over spins)
            end do
            end do
            !
            cm(:, 1:nbspx) = wfc_ctmp(:, 1:nbspx)
            !
         end if
         !
      end if
      !
      deallocate (Omat1tot)
      deallocate (Umatbig)
      deallocate (Heigbig)
      deallocate (wfc_ctmp)
      deallocate (vsic1_realspace)
      deallocate (vsic_reciprocal1)
      call deallocate_twin(bec1)
      deallocate (pink1)
      !
      call stop_clock("nk_rot_cm")
      call stop_clock('nk_rot_emin')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rot_emin
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rot_test(passoprod, nsteps, nouter, ninner, etot)
!-----------------------------------------------------------------------
!
! ... prints out esic by varying the wavefunction along a search direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use io_global, only: ionode
      use electrons_base, only: nbsp, nbspx, nspin, &
                                iupdwn, nupdwn
      use uspp, only: nkb
      use wavefunctions_module, only: c0
      !
      implicit none
      !
      ! in/out vars
      !
      real(dp), intent(in)  :: passoprod
      integer, intent(in)  :: nsteps, ninner, nouter
      real(dp), intent(in)  :: etot

      !
      ! local variables
      !
      real(dp) :: esic
      real(dp)                 :: bec1(nkb, nbsp)
      real(dp)                 :: Omat1tot(nbspx, nbspx)
      real(dp)                 :: vsic1_realspace(nnrx, nbspx)
      complex(dp)              :: vsic_reciprocal1(ngm, nbspx)
      complex(dp), allocatable :: Umat(:, :)
      complex(dp)              :: Umatbig(nbspx, nbspx)
      real(dp), allocatable    :: Heig(:)
      real(dp)                 :: Heigbig(nbspx)
      complex(dp)              :: wfc_ctmp(ngw, nbspx)
      real(dp)                 :: dalpha, dmaxeig
      real(dp)                 :: pink1(nbspx)
      integer                  :: isp, istep
      real(dp), allocatable    :: vsicah(:, :)
      real(dp)                 :: vsicah2sum, deigrms
      integer                  :: nfile

      !
      ! variables for test calculations - along gradient line direction
      !

      !
      ! main body
      !
      CALL start_clock('nk_rot_test')

      Umatbig(:, :) = (0.d0, 0.d0)
      Heigbig(:) = 0.d0
      deigrms = 0.d0

      do isp = 1, nspin

         allocate (Umat(nupdwn(isp), nupdwn(isp)))
         allocate (Heig(nupdwn(isp)))
         allocate (vsicah(nupdwn(isp), nupdwn(isp)))

         call nksic_getvsicah(isp, vsicah, vsicah2sum)
         call nksic_getHeigU(isp, vsicah, Heig, Umat)

         Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Umat(:, :)
         Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Heig(:)

         deigrms = deigrms + sum(Heig(:)**2)

         deallocate (Umat)
         deallocate (Heig)
         deallocate (vsicah)

      end do ! do isp=1,nspin

      ! how severe the transform is
      deigrms = sqrt(deigrms/nbsp)

      dmaxeig = max(abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
      do isp = 2, nspin
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp))))
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
      end do

      nfile = 10000 + 100*nouter + ninner
      if (ionode) write (nfile, *) '# passoprod', passoprod

      do istep = 1, nsteps
         if (nsteps .ne. 1) then
            dalpha = passoprod*(2.d0*istep - nsteps - 1.d0)/(nsteps - 1.d0)/dmaxeig
         else
            dalpha = 0.d0
         end if

         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp, Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, esic)

         if (ionode) write (nfile, '(5F24.13,2I10)') dalpha/3.141592*dmaxeig, dmaxeig, etot, esic, deigrms, ninner, nouter

      end do  !$$ do istep=1,nsteps

      if (ionode) write (nfile, *)

      CALL stop_clock('nk_rot_test')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rot_test
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rot_emin_cg_new(c0, cm, vsic_realspace, vsic_reciprocal, ngw, nnrx, bec, &
                                    nouter, init_n, ninner, etot, Omattot, &
                                    rot_threshold, nbsp, nbspx, nudx, nspin, iupdwn, &
                                    nupdwn, pink, wfc_centers, wfc_spreads, lgam)
!-----------------------------------------------------------------------
!
! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly using cg minimization.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use io_global, only: stdout, ionode
      use gvecp, only: ngm
      use cp_interfaces, only: invfft
      use nksic, only: innerloop_cg_nsd, &
                       innerloop_cg_nreset, &
                       innerloop_nmax, &
                       innerloop_atleast
      use uspp, only: nkb
      use control_flags, only: esic_conv_thr
      use cg_module, only: tcg
      use twin_types
      !
      implicit none
      !
      ! in/out vars
      !
      integer                  :: ninner, nbsp, nbspx, nspin, nudx, nnrx
      integer                  :: init_n, ngw, ispin(nbspx)
      integer, intent(in)  :: nouter
      integer, intent(in)  :: iupdwn(nspin), nupdwn(nspin)
      real(dp), intent(in)  :: etot
      complex(dp)          :: Omattot(nbspx, nbspx), c0(ngw, nbsp), cm(ngw, nbsp)
      real(dp), intent(in)  :: rot_threshold
      real(dp), intent(inout) :: pink(nbsp), vsic_realspace(nnrx, nbspx), vsic_reciprocal(ngm, nbspx), &
                                 wfc_centers(4, nudx, nspin), wfc_spreads(nudx, nspin, 2)
      logical               :: lgam
      type(twin_matrix)     :: bec

      !
      ! local variables for cg routine
      !
      integer     :: nbnd1, nbnd2
      integer     :: isp
      logical     :: ldotest
      real(dp)    :: dtmp
      real(dp)    :: ene0, ene1, enesti, enever, dene0
      real(dp)    :: passo, passov, passof, passomax, spasso
      real(dp)    :: vsicah2sum, vsicah2sum_prev
      integer     :: nidx1, nidx2
      real(dp)    :: dPI, dalpha, dmaxeig, deigrms
      real(dp)    :: pinksumprev, passoprod
      !
      complex(dp), allocatable :: Omat1tot(:, :), Omat2tot(:, :)
      real(dp), allocatable :: Heigbig(:)
      complex(dp), allocatable :: Umatbig(:, :)
      complex(dp), allocatable :: wfc_ctmp(:, :), wfc_ctmp2(:, :)
      complex(dp), allocatable :: gi(:, :), hi(:, :)
      !
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)
      complex(dp), allocatable :: vsicah(:, :)
      real(dp), allocatable :: vsic1_realspace(:, :), vsic2_realspace(:, :)
      complex(dp), allocatable :: vsic_reciprocal1(:, :), vsic_reciprocal2(:, :)
      type(twin_matrix)       :: bec1, bec2
      real(dp), allocatable :: pink1(:), pink2(:)
      logical :: restartcg_innerloop, ene_ok_innerloop, ltresh, setpassomax
      integer :: iter3, nfail
      integer :: maxiter3, numok, minsteps
      real(dp) :: signalpha
      character(len=4) :: marker
      real(dp) :: conv_thr

      !
      ! main body
      !
      CALL start_clock('nk_rot_emin')
      !
      !
      marker = "   "
      maxiter3 = 4
      minsteps = 2
      restartcg_innerloop = .true.
      ene_ok_innerloop = .false.
      ltresh = .false.
      setpassomax = .false.
      nfail = 0
      if (nouter < init_n) THEN
         conv_thr = esic_conv_thr
      ELSE
         conv_thr = rot_threshold
      END IF
      !
      pinksumprev = 1.d8
      dPI = 2.0_DP*asin(1.0_DP)
      passoprod = 0.3d0

      !
      ! local workspace
      !
      allocate (Omat1tot(nbspx, nbspx), Omat2tot(nbspx, nbspx))
      allocate (Umatbig(nbspx, nbspx))
      allocate (Heigbig(nbspx))
      allocate (wfc_ctmp(ngw, nbspx), wfc_ctmp2(ngw, nbspx))
      allocate (hi(nbsp, nbsp))
      allocate (gi(nbsp, nbsp))
      allocate (pink1(nbspx), pink2(nbspx))
      allocate (vsic1_realspace(nnrx, nbspx), vsic2_realspace(nnrx, nbspx))
      allocate (vsic_reciprocal1(ngm, nbspx), vsic_reciprocal2(ngm, nbspx))
      call init_twin(bec1, lgam)
      call allocate_twin(bec1, nkb, nbsp, lgam)
      call init_twin(bec2, lgam)
      call allocate_twin(bec2, nkb, nbsp, lgam)
      !
      Umatbig(:, :) = CMPLX(0.d0, 0.d0)
      Heigbig(:) = 0.d0
      deigrms = 0.d0
      hi(:, :) = 0.d0
      gi(:, :) = 0.d0

      Omattot(:, :) = CMPLX(0.d0, 0.d0)
      do nbnd1 = 1, nbspx
         Omattot(nbnd1, nbnd1) = CMPLX(1.d0, 0.d0)
      end do

      ninner = 0
      ldotest = .false.

      if (ionode) write (stdout, "(14x,'# iter',6x,'etot',17x,'esic',&
                                  & 17x,'deigrms')")

      !
      ! main loop
      !
      inner_loop: &
         do while (.true.)

         call start_clock("nk_innerloop")
         !
         ninner = ninner + 1

         if (ninner > innerloop_nmax) then
            !
#ifdef __DEBUG
            if (ionode) write (1031, *) '# innerloop_nmax reached.'
            if (ionode) write (1031, *)
#endif
            if (ionode) then
               write (stdout, "(14x,'# innerloop_nmax reached.',/)")
            end if
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if

#ifdef __DEBUG
         !
         ! call nksic_printoverlap(ninner,nouter)

!        if(mod(ninner,10).eq.1.or.ninner.le.5) ldotest=.true.
         if (ninner .eq. 31 .or. ninner .eq. 61 .or. ninner .eq. 91) ldotest = .true.
!        if(ninner.le.10.and.nouter.eq.1) ldotest=.true.
!         ldotest=.true.
!        if(ninner.ge.25) ldotest=.true.
         ! Now do the test
         if (ldotest) then
!          dtmp = 1.0d0*3.141592d0
            dtmp = 4.d0*3.141592d0
!          call nksic_rot_test(dtmp,201,nouter,ninner,etot)
            ldotest = .false.
         end if
#endif

         !
         !print out ESIC part & other total energy
         !
         ene0 = sum(pink(1:nbsp))

         !
         ! test convergence
         !
         if (abs(ene0 - pinksumprev) < conv_thr) then
            numok = numok + 1
         else
            numok = 0
         end if
         !
         if (numok >= minsteps .and. ninner >= innerloop_atleast) ltresh = .true.
         !
         if (ltresh) then
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, "(a,/)") '# inner-loop converged.'
               write (1031, "(a,/)") '# inner-loop converged.'
            end if
#endif
            if (ionode .and. numok < minsteps) then
               write (stdout, "(14x,'# innerloop nstep reached',/)")
            else
               write (stdout, "(14x,'# innerloop converged',/)")
            end if
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if
         !
         pinksumprev = ene0

         !
         ! This part calculates the anti-hermitian part of the Hamiltonian vsicah
         ! and see whether a convergence has been achieved
         !
         ! For this run, we obtain the gradient
         !
         vsicah2sum = 0.0d0
         !
         do isp = 1, nspin
            !
            allocate (vsicah(nupdwn(isp), nupdwn(isp)))
            !
            call nksic_getvsicah_new3(ngw, nbsp, nbspx, nspin, c0, &
                                      bec, isp, nupdwn, iupdwn, vsicah, dtmp, lgam)
            !
            gi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
               iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = vsicah(:, :)
            !
            vsicah2sum = vsicah2sum + dtmp
            !
            deallocate (vsicah)
            !
         end do
         !
         if (ninner /= 1) dtmp = vsicah2sum/vsicah2sum_prev
         !
         if (ninner <= innerloop_cg_nsd .or. &
             mod(ninner, innerloop_cg_nreset) == 0 .or. &
             restartcg_innerloop) then
            !
            restartcg_innerloop = .false.
            setpassomax = .false.
            !
            hi(:, :) = gi(:, :)
         else
            hi(:, :) = gi(:, :) + dtmp*hi(:, :)
         end if
         !
         spin_loop: &
            do isp = 1, nspin
            !
            IF (nupdwn(isp) .gt. 0) THEN
               allocate (vsicah(nupdwn(isp), nupdwn(isp)))
               allocate (Umat(nupdwn(isp), nupdwn(isp)))
               allocate (Heig(nupdwn(isp)))

               vsicah(:, :) = hi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                                 iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))

               call nksic_getHeigU_new(nspin, isp, nupdwn, vsicah, Heig, Umat)
               !
               deigrms = deigrms + sum(Heig(:)**2)
               !
               Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                       iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Umat(:, :)
               Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Heig(:)
               !
               deallocate (vsicah)
               deallocate (Umat)
               deallocate (Heig)
            ELSE
               Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                       iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 1.d0
               Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 0.d0
            END IF
            !
         end do spin_loop

         ! how severe the transform is
         deigrms = sqrt(deigrms/nbsp)
#ifdef __DEBUG
         if (ionode) write (1031, '(2I10,3F24.13)') ninner, nouter, etot, ene0, deigrms
#endif
         if (ionode) write (stdout, '(10x,A3,2i5,3F21.13)') marker, ninner, nouter, etot, ene0, deigrms
         !
         !
         dmaxeig = max(dabs(Heigbig(iupdwn(1))), dabs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
         !
         do isp = 2, nspin
            dmaxeig = max(dmaxeig, dabs(Heigbig(iupdwn(isp))))
            dmaxeig = max(dmaxeig, dabs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
         end do
         !
         passomax = passoprod/dmaxeig
         !
         if (ninner == 1 .or. setpassomax) then
            passof = passomax
            setpassomax = .false.
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passof set to passomax'
#endif
         end if

!$$$$        if(passof .gt. passomax*2.d0) then
!$$$$          passof = passomax*2.d0
!$$$$          if(ionode) write(1031,*) '# passof > twice passomax'
!$$$$        endif

!        if(ionode) then
!          write(1037,*)'# deigrms = ',deigrms
!          write(1037,*)'# vsicah2sum = ',vsicah2sum
!          if(ninner.ne.1) write(1037,*)'# vsicah2sum/vsicah2sum_prev = ',dtmp
!        endif

         vsicah2sum_prev = vsicah2sum
         !
         dene0 = 0.d0
         !
         do isp = 1, nspin
            !
            do nbnd1 = 1, nupdwn(isp)
            do nbnd2 = 1, nupdwn(isp)
               !
               nidx1 = nbnd1 - 1 + iupdwn(isp)
               nidx2 = nbnd2 - 1 + iupdwn(isp)
               IF (nidx1 .ne. nidx2) THEN
                  dene0 = dene0 - DBLE(CONJG(gi(nidx1, nidx2))*hi(nidx1, nidx2))
               ELSE  !warning:giovanni: do we need this condition
                  !dene0 = dene0 -DBLE(CONJG(gi(nidx1,nidx2))*hi(nidx1,nidx2))
               END IF
               !
            end do
            end do
            !
         end do

         !$$
         !$$ dene0 = dene0 * 2.d0/nspin
         !
         ! Be careful, the following is correct because A_ji = - A_ij, i.e., the number of
         ! linearly independent variables is half the number of total variables!
         !
         dene0 = dene0*1.d0/nspin
         !
         spasso = 1.d0
         if (dene0 > 0.d0) spasso = -1.d0
         !
         dalpha = spasso*passof
         !
         call nksic_getOmattot_new(nbsp, nbspx, nudx, nspin, ispin, &
                                   iupdwn, nupdwn, wfc_centers, wfc_spreads, &
                                   dalpha, Heigbig, Umatbig, c0, wfc_ctmp, &
                                   Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, ene1, lgam)
         call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)
         !
         ! We neglect this step for paper writing purposes
         !
         if (passo > passomax) then
            passo = passomax
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passo > passomax'
#endif
            !
         end if

         passov = passof
         passof = 2.d0*passo

         dalpha = spasso*passo
         !
!$$ The following line is for dene0 test
!        if(ninner.ge.15) dalpha = spasso*passo*0.00001
!$$
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp2, &
                               Omat2tot, bec2, vsic2_realspace, vsic_reciprocal2, pink2, enever, lgam)
#ifdef __DEBUG
         if (ionode) then
            !
            write (1037, *) ninner, nouter
            write (1037, '("ene0,ene1,enesti,enever")')
            write (1037, '(a3,4f20.10)') 'CG1', ene0, ene1, enesti, enever
            write (1037, '("spasso,passov,passo,passomax,dene0,&
                         & (enever-ene0)/passo/dene0")')
            write (1037, '(a3,4f12.7,e20.10,f12.7)') &
               'CG2', spasso, passov, passo, passomax, dene0, (enever - ene0)/passo/dene0
            write (1037, *)
            !
         end if
#endif
         if (ene0 < ene1 .and. ene0 < enever) then !missed minimum case 3
            !write(6,'("# WARNING: innerloop missed minimum, case 3",/)')
            !
            iter3 = 0
            signalpha = 1.d0
            restartcg_innerloop = .true.
            !
            do while (enever .ge. ene0 .and. iter3 .lt. maxiter3)
               !
               iter3 = iter3 + 1
               !
               signalpha = signalpha*(-0.717d0)
               dalpha = spasso*passo*signalpha
               !
        call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp2, Omat2tot, bec2, vsic2_realspace, vsic_reciprocal2, pink2, enever, lgam)
               !
            end do

            IF (enever .lt. ene0) THEN
               !
               pink(:) = pink2(:)
               vsic_realspace(:, :) = vsic2_realspace(:, :)
               vsic_reciprocal(:, :) = vsic_reciprocal2(:, :)
               c0(:, :) = wfc_ctmp2(:, :)
               call copy_twin(bec, bec2)
               !             bec%rvec(:,:)  = bec2(:,:)
               Omattot = MATMUL(Omattot, Omat2tot)
               !write(6,'("# WARNING: innerloop case 3 interations",3I/)') iter3
               write (marker, '(i1)') iter3
               marker = '*'//marker
               passof = passo*abs(signalpha)
               nfail = 0
               !
            ELSE
               !
               marker = '***'
               ninner = ninner + 1
               nfail = nfail + 1
               numok = 0
               passof = passo*abs(signalpha)
               !
               IF (nfail > 2) THEN
                  write (6, '("# WARNING: innerloop not converged, exit",/)')
                  call stop_clock("nk_innerloop")
                  exit
               END IF
!                ELSE
!                   nfail=0
!                ENDIF
               !
            END IF
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# ene0<ene1 and ene0<enever, exit",/)')
               write (1031, '("# innerloop NOT converged, exit",/)')
            end if
#endif

            !
         else if (ene1 >= enever) then !found minimum
            !
            pink(:) = pink2(:)
            vsic_realspace(:, :) = vsic2_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal2(:, :)
            c0(:, :) = wfc_ctmp2(:, :)
            call copy_twin(bec, bec2)
!             bec%rvec(:,:)  = bec2(:,:)
            Omattot = MATMUL(Omattot, Omat2tot)
            marker = "   "
            nfail = 0
            !
         else !missed minimum, case 1 or 2
            !
            pink(:) = pink1(:)
            vsic_realspace(:, :) = vsic1_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal1(:, :)
            c0(:, :) = wfc_ctmp(:, :)
            call copy_twin(bec, bec1)
            Omattot = MATMUL(Omattot, Omat1tot)
            restartcg_innerloop = .true.
            IF (enever < ene0) THEN
               marker = "*  "
               passof = min(1.5d0*passov, passomax)
            ELSE
               marker = "** "
               passof = passov
            END IF
            nfail = 0
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# It happened that ene1 < enever!!")')
               write (1037, *)
            end if
#endif
            !write(6,'("# WARNING: innerloop missed minimum case 1 or 2",/)')
            !
! =======
!           pink(:) = pink1(:)
!           vsic_realspace(:,:) = vsic1_realspace(:,:)
!           c0(:,:) = wfn_ctmp(:,:)
!           bec%rvec(:,:) = bec1(:,:)
!           Omattot = MATMUL(Omattot,Omat1tot)
!           if(ionode) then
!             write(1037,'("# It happened that ene1 < enever!!")')
!             write(1037,*)
!           endif
! 1.28.2.14
         end if
         !
         call stop_clock("nk_innerloop")
         !
      end do inner_loop

      !
      ! Wavefunction cm rotation according to Omattot
      ! We need this because outer loop could be damped dynamics.
      !
      if (.not. tcg) then
         !
         if (ninner >= 2) then
            !
            wfc_ctmp(:, :) = CMPLX(0.d0, 0.d0)
            !
            do nbnd1 = 1, nbspx
               do nbnd2 = 1, nbspx
                  wfc_ctmp(:, nbnd1) = wfc_ctmp(:, nbnd1) + cm(:, nbnd2)*Omattot(nbnd2, nbnd1) !warning:giovanni CONJUGATE?
                  ! XXX (we can think to use a blas, here, and split over spins)
                  !does not seem we need to make it conjugate
               end do
            end do
            !
            cm(:, 1:nbspx) = wfc_ctmp(:, 1:nbspx)
            !
         end if
         !
      end if
      !
      ! clean local workspace
      !
      deallocate (Omat1tot, Omat2tot)
      deallocate (Umatbig)
      deallocate (Heigbig)
      deallocate (wfc_ctmp, wfc_ctmp2)
      deallocate (hi)
      deallocate (gi)
      deallocate (pink1, pink2)
      deallocate (vsic1_realspace, vsic2_realspace)
      deallocate (vsic_reciprocal1, vsic_reciprocal2)
      call deallocate_twin(bec1)
      call deallocate_twin(bec2)

      CALL stop_clock('nk_rot_emin')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rot_emin_cg_new
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rot_emin_cg(nouter, init_n, ninner, etot, Omattot, &
                                rot_threshold, lgam)
!-----------------------------------------------------------------------
!
! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly using cg minimization.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use io_global, only: stdout, ionode
      use electrons_base, only: nbsp, nbspx, nspin, &
                                iupdwn, nupdwn
      use cp_interfaces, only: invfft
      use nksic, only: vsic_realspace, vsic_reciprocal, pink, &
                       innerloop_cg_nsd, innerloop_cg_nreset, &
                       innerloop_nmax, &
                       innerloop_atleast
      use uspp, only: nkb
      use cp_main_variables, only: bec
      use wavefunctions_module, only: c0, cm
      use control_flags, only: esic_conv_thr
      use cg_module, only: tcg
      use twin_types
      !
      implicit none
      !
      ! in/out vars
      !
      integer                  :: ninner
      integer                  :: init_n
      integer, intent(in)  :: nouter
      real(dp), intent(in)  :: etot
      complex(dp)          :: Omattot(nbspx, nbspx)
      real(dp), intent(in)  :: rot_threshold
      logical               :: lgam

      !
      ! local variables for cg routine
      !
      integer     :: nbnd1, nbnd2
      integer     :: isp
      logical     :: ldotest
      real(dp)    :: dtmp
      real(dp)    :: ene0, ene1, enesti, enever, dene0
      real(dp)    :: passo, passov, passof, passomax, spasso
      real(dp)    :: vsicah2sum, vsicah2sum_prev
      integer     :: nidx1, nidx2
      real(dp)    :: dPI, dalpha, dmaxeig, deigrms
      real(dp)    :: pinksumprev, passoprod
      !
      complex(dp), allocatable :: Omat1tot(:, :), Omat2tot(:, :)
      real(dp), allocatable :: Heigbig(:)
      complex(dp), allocatable :: Umatbig(:, :)
      complex(dp), allocatable :: wfc_ctmp(:, :), wfc_ctmp2(:, :)
      complex(dp), allocatable :: gi(:, :), hi(:, :)
      !
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)
      complex(dp), allocatable :: vsicah(:, :)
      real(dp), allocatable :: vsic1_realspace(:, :), vsic2_realspace(:, :)
      complex(dp), allocatable :: vsic_reciprocal1(:, :), vsic_reciprocal2(:, :)
      type(twin_matrix)       :: bec1, bec2
      real(dp), allocatable :: pink1(:), pink2(:)
      logical :: restartcg_innerloop, ene_ok_innerloop, ltresh, setpassomax
      integer :: iter3, nfail
      integer :: maxiter3, numok
      real(dp) :: signalpha
      character(len=4) :: marker
      real(dp) :: conv_thr

      !
      ! main body
      !
      CALL start_clock('nk_rot_emin')
      !
      !
      marker = "   "
      maxiter3 = 4
      restartcg_innerloop = .true.
      ene_ok_innerloop = .false.
      ltresh = .false.
      setpassomax = .false.
      nfail = 0
      if (nouter < init_n) THEN
         conv_thr = esic_conv_thr
      ELSE
         conv_thr = rot_threshold
      END IF
      !
      pinksumprev = 1.d8
      dPI = 2.0_DP*asin(1.0_DP)
      passoprod = 0.3d0

      !
      ! local workspace
      !
      allocate (Omat1tot(nbspx, nbspx), Omat2tot(nbspx, nbspx))
      allocate (Umatbig(nbspx, nbspx))
      allocate (Heigbig(nbspx))
      allocate (wfc_ctmp(ngw, nbspx), wfc_ctmp2(ngw, nbspx))
      allocate (hi(nbsp, nbsp))
      allocate (gi(nbsp, nbsp))
      allocate (pink1(nbspx), pink2(nbspx))
      allocate (vsic1_realspace(nnrx, nbspx), vsic2_realspace(nnrx, nbspx))
      allocate (vsic_reciprocal1(ngm, nbspx), vsic_reciprocal2(ngm, nbspx))
      call init_twin(bec1, lgam)
      call allocate_twin(bec1, nkb, nbsp, lgam)
      call init_twin(bec2, lgam)
      call allocate_twin(bec2, nkb, nbsp, lgam)
      !
      Umatbig(:, :) = CMPLX(0.d0, 0.d0)
      Heigbig(:) = 0.d0
      deigrms = 0.d0
      hi(:, :) = 0.d0
      gi(:, :) = 0.d0

      Omattot(:, :) = CMPLX(0.d0, 0.d0)
      do nbnd1 = 1, nbspx
         Omattot(nbnd1, nbnd1) = CMPLX(1.d0, 0.d0)
      end do

      ninner = 0
      ldotest = .false.

      if (ionode) write (stdout, "(14x,'# iter',6x,'etot',17x,'esic',&
                                  & 17x,'deigrms')")

      !
      ! main loop
      !
      inner_loop: &
         do while (.true.)

         call start_clock("nk_innerloop")
         !
         ninner = ninner + 1

         if (ninner > innerloop_nmax) then
            !
#ifdef __DEBUG
            if (ionode) write (1031, *) '# innerloop_nmax reached.'
            if (ionode) write (1031, *)
#endif
            if (ionode) then
               write (stdout, "(14x,'# innerloop_nmax reached.',/)")
            end if
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if

#ifdef __DEBUG
         !
         ! call nksic_printoverlap(ninner,nouter)

!        if(mod(ninner,10).eq.1.or.ninner.le.5) ldotest=.true.
         if (ninner .eq. 31 .or. ninner .eq. 61 .or. ninner .eq. 91) ldotest = .true.
!        if(ninner.le.10.and.nouter.eq.1) ldotest=.true.
!         ldotest=.true.
!        if(ninner.ge.25) ldotest=.true.
         ! Now do the test
         if (ldotest) then
!          dtmp = 1.0d0*3.141592d0
            dtmp = 4.d0*3.141592d0
!          call nksic_rot_test(dtmp,201,nouter,ninner,etot)
            ldotest = .false.
         end if
#endif

         !
         !print out ESIC part & other total energy
         !
         ene0 = sum(pink(1:nbsp))

         !
         ! test convergence
         !
         if (abs(ene0 - pinksumprev) < conv_thr) then
            numok = numok + 1
         else
            numok = 0
         end if
         !
         if (numok >= 2 .and. ninner >= innerloop_atleast) ltresh = .true.
         !
         if (ltresh) then
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, "(a,/)") '# inner-loop converged.'
               write (1031, "(a,/)") '# inner-loop converged.'
            end if
#endif
            if (ionode) write (stdout, "(14x,'# innerloop converged',/)")
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if
         !
         pinksumprev = ene0

         !
         ! This part calculates the anti-hermitian part of the Hamiltonian vsicah
         ! and see whether a convergence has been achieved
         !
         ! For this run, we obtain the gradient
         !
         vsicah2sum = 0.0d0
         !
         do isp = 1, nspin
            !
            allocate (vsicah(nupdwn(isp), nupdwn(isp)))
            !
            call nksic_getvsicah_new2(isp, vsicah, dtmp, lgam)
            !
            gi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
               iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = vsicah(:, :)
            !
            vsicah2sum = vsicah2sum + dtmp
            !
            deallocate (vsicah)
            !
         end do
         !
         if (ninner /= 1) dtmp = vsicah2sum/vsicah2sum_prev
         !
         if (ninner <= innerloop_cg_nsd .or. &
             mod(ninner, innerloop_cg_nreset) == 0 .or. &
             restartcg_innerloop) then
            !
            restartcg_innerloop = .false.
            setpassomax = .false.
            !
            hi(:, :) = gi(:, :)
         else
            hi(:, :) = gi(:, :) + dtmp*hi(:, :)
         end if
         !
         spin_loop: &
            do isp = 1, nspin
            !
            IF (nupdwn(isp) .gt. 0) THEN
               allocate (vsicah(nupdwn(isp), nupdwn(isp)))
               allocate (Umat(nupdwn(isp), nupdwn(isp)))
               allocate (Heig(nupdwn(isp)))

               vsicah(:, :) = hi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                                 iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))

               call nksic_getHeigU(isp, vsicah, Heig, Umat)
               !
               deigrms = deigrms + sum(Heig(:)**2)
               !
               Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                       iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Umat(:, :)
               Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Heig(:)
               !
               deallocate (vsicah)
               deallocate (Umat)
               deallocate (Heig)
            ELSE
               Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                       iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 1.d0
               Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 0.d0
            END IF
            !
         end do spin_loop

         ! how severe the transform is
         deigrms = sqrt(deigrms/nbsp)
#ifdef __DEBUG
         if (ionode) write (1031, '(2I10,3F24.13)') ninner, nouter, etot, ene0, deigrms
#endif
         if (ionode) write (stdout, '(10x,A3,2i5,3F21.13)') marker, ninner, nouter, etot, ene0, deigrms
         !
         !
         dmaxeig = max(dabs(Heigbig(iupdwn(1))), dabs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
         !
         do isp = 2, nspin
            dmaxeig = max(dmaxeig, dabs(Heigbig(iupdwn(isp))))
            dmaxeig = max(dmaxeig, dabs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
         end do
         !
         passomax = passoprod/dmaxeig
         !
         if (ninner == 1 .or. setpassomax) then
            passof = passomax
            setpassomax = .false.
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passof set to passomax'
#endif
         end if

!$$$$        if(passof .gt. passomax*2.d0) then
!$$$$          passof = passomax*2.d0
!$$$$          if(ionode) write(1031,*) '# passof > twice passomax'
!$$$$        endif

!        if(ionode) then
!          write(1037,*)'# deigrms = ',deigrms
!          write(1037,*)'# vsicah2sum = ',vsicah2sum
!          if(ninner.ne.1) write(1037,*)'# vsicah2sum/vsicah2sum_prev = ',dtmp
!        endif

         vsicah2sum_prev = vsicah2sum
         !
         dene0 = 0.d0
         !
         do isp = 1, nspin
            !
            do nbnd1 = 1, nupdwn(isp)
            do nbnd2 = 1, nupdwn(isp)
               !
               nidx1 = nbnd1 - 1 + iupdwn(isp)
               nidx2 = nbnd2 - 1 + iupdwn(isp)
               IF (nidx1 .ne. nidx2) THEN
                  dene0 = dene0 - DBLE(CONJG(gi(nidx1, nidx2))*hi(nidx1, nidx2))
               ELSE  !warning:giovanni: do we need this condition
                  !dene0 = dene0 -DBLE(CONJG(gi(nidx1,nidx2))*hi(nidx1,nidx2))
               END IF
               !
            end do
            end do
            !
         end do

         !$$
         !$$ dene0 = dene0 * 2.d0/nspin
         !
         ! Be careful, the following is correct because A_ji = - A_ij, i.e., the number of
         ! linearly independent variables is half the number of total variables!
         !
         dene0 = dene0*1.d0/nspin
         !
         spasso = 1.d0
         if (dene0 > 0.d0) spasso = -1.d0
         !
         dalpha = spasso*passof
         !
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp, &
                               Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, ene1, lgam)
         call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)

         !
         ! We neglect this step for paper writing purposes
         !
         if (passo > passomax) then
            passo = passomax
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passo > passomax'
#endif
            !
         end if

         passov = passof
         passof = 2.d0*passo

         dalpha = spasso*passo
         !
!$$ The following line is for dene0 test
!        if(ninner.ge.15) dalpha = spasso*passo*0.00001
!$$
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp2, &
                               Omat2tot, bec2, vsic2_realspace, vsic_reciprocal2, pink2, enever, lgam)

#ifdef __DEBUG
         if (ionode) then
            !
            write (1037, *) ninner, nouter
            write (1037, '("ene0,ene1,enesti,enever")')
            write (1037, '(a3,4f20.10)') 'CG1', ene0, ene1, enesti, enever
            write (1037, '("spasso,passov,passo,passomax,dene0,&
                         & (enever-ene0)/passo/dene0")')
            write (1037, '(a3,4f12.7,e20.10,f12.7)') &
               'CG2', spasso, passov, passo, passomax, dene0, (enever - ene0)/passo/dene0
            write (1037, *)
            !
         end if
#endif

         if (ene0 < ene1 .and. ene0 < enever) then !missed minimum case 3
            !write(6,'("# WARNING: innerloop missed minimum, case 3",/)')
            !
            iter3 = 0
            signalpha = 1.d0
            restartcg_innerloop = .true.
            !
            do while (enever .ge. ene0 .and. iter3 .lt. maxiter3)
               !
               iter3 = iter3 + 1
               !
               signalpha = signalpha*(-0.717d0)
               dalpha = spasso*passo*signalpha
               !
        call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp2, Omat2tot, bec2, vsic2_realspace, vsic_reciprocal2, pink2, enever, lgam)
               !
            end do

            IF (enever .lt. ene0) THEN
               !
               pink(:) = pink2(:)
               vsic_realspace(:, :) = vsic2_realspace(:, :)
               vsic_reciprocal(:, :) = vsic_reciprocal2(:, :)
               c0(:, :) = wfc_ctmp2(:, :)
               call copy_twin(bec, bec2)
               !             bec%rvec(:,:)  = bec2(:,:)
               Omattot = MATMUL(Omattot, Omat2tot)
               !write(6,'("# WARNING: innerloop case 3 interations",3I/)') iter3
               write (marker, '(i1)') iter3
               marker = '*'//marker
               passof = passo*abs(signalpha)
               nfail = 0
               !
            ELSE
               !
               marker = '***'
               ninner = ninner + 1
               nfail = nfail + 1
               numok = 0
               passof = passo*abs(signalpha)
               !
               IF (nfail > 2) THEN
                  write (6, '("# WARNING: innerloop not converged, exit",/)')
                  call stop_clock("nk_innerloop")
                  exit
               END IF
!                ELSE
!                   nfail=0
!                ENDIF
               !
            END IF
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# ene0<ene1 and ene0<enever, exit",/)')
               write (1031, '("# innerloop NOT converged, exit",/)')
            end if
#endif

            !
         else if (ene1 >= enever) then !found minimum
            !
            pink(:) = pink2(:)
            vsic_realspace(:, :) = vsic2_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal2(:, :)
            c0(:, :) = wfc_ctmp2(:, :)
            call copy_twin(bec, bec2)
!             bec%rvec(:,:)  = bec2(:,:)
            Omattot = MATMUL(Omattot, Omat2tot)
            marker = "   "
            nfail = 0
            !
         else !missed minimum, case 1 or 2
            !
            pink(:) = pink1(:)
            vsic_realspace(:, :) = vsic1_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal1(:, :)
            c0(:, :) = wfc_ctmp(:, :)
            call copy_twin(bec, bec1)
            Omattot = MATMUL(Omattot, Omat1tot)
            restartcg_innerloop = .true.
            IF (enever < ene0) THEN
               marker = "*  "
               passof = min(1.5d0*passov, passomax)
            ELSE
               marker = "** "
               passof = passov
            END IF
            nfail = 0
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# It happened that ene1 < enever!!")')
               write (1037, *)
            end if
#endif
            !write(6,'("# WARNING: innerloop missed minimum case 1 or 2",/)')
            !
! =======
!           pink(:) = pink1(:)
!           vsic_realspace(:,:) = vsic1_realspace(:,:)
!           c0(:,:) = wfn_ctmp(:,:)
!           bec%rvec(:,:) = bec1(:,:)
!           Omattot = MATMUL(Omattot,Omat1tot)
!           if(ionode) then
!             write(1037,'("# It happened that ene1 < enever!!")')
!             write(1037,*)
!           endif
! 1.28.2.14
         end if
         !
         call stop_clock("nk_innerloop")
         !
      end do inner_loop

      !
      ! Wavefunction cm rotation according to Omattot
      ! We need this because outer loop could be damped dynamics.
      !
      if (.not. tcg) then
         !
         if (ninner >= 2) then
            !
            wfc_ctmp(:, :) = CMPLX(0.d0, 0.d0)
            !
            do nbnd1 = 1, nbspx
               do nbnd2 = 1, nbspx
                  wfc_ctmp(:, nbnd1) = wfc_ctmp(:, nbnd1) + cm(:, nbnd2)*Omattot(nbnd2, nbnd1) !warning:giovanni CONJUGATE?
                  ! XXX (we can think to use a blas, here, and split over spins)
                  !does not seem we need to make it conjugate
               end do
            end do
            !
            cm(:, 1:nbspx) = wfc_ctmp(:, 1:nbspx)
            !
         end if
         !
      end if

      !
      ! clean local workspace
      !
      deallocate (Omat1tot, Omat2tot)
      deallocate (Umatbig)
      deallocate (Heigbig)
      deallocate (wfc_ctmp, wfc_ctmp2)
      deallocate (hi)
      deallocate (gi)
      deallocate (pink1, pink2)
      deallocate (vsic1_realspace, vsic2_realspace)
      deallocate (vsic_reciprocal1, vsic_reciprocal2)
      call deallocate_twin(bec1)
      call deallocate_twin(bec2)
!       deallocate( bec1, bec2 )

      CALL stop_clock('nk_rot_emin')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rot_emin_cg
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rot_emin_cg_descla(nouter, ninner, etot, Omattot, lgam)
!-----------------------------------------------------------------------
! !warning:giovanni why not passing wavefunctions as variables???

! ... Finds the orthogonal rotation matrix Omattot that minimizes
!     the orbital-dependent and hence the total energy, and then
!     rotate the wavefunction c0 accordingly using cg minimization.
!     We may need Omattot for further rotation of the gradient for outer loop CG.
!     Right now we do not do that because we set resetcg=.true. after inner loop
!     minimization routine, i.e., setting the search direction to be gradient direction.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use io_global, only: stdout, ionode
      use electrons_base, only: nbsp, nbspx, nspin, &
                                iupdwn, nupdwn
      use cp_interfaces, only: invfft
      use nksic, only: vsic_realspace, vsic_reciprocal, pink, &
                       innerloop_cg_nsd, innerloop_cg_nreset, &
                       innerloop_nmax
      use uspp, only: nkb
      use cp_main_variables, only: bec
      use wavefunctions_module, only: c0, cm
      use control_flags, only: esic_conv_thr
      use cg_module, only: tcg
      use twin_types
      !
      implicit none
      !
      ! in/out vars
      !
      integer                  :: ninner
      integer, intent(in)  :: nouter
      real(dp), intent(in)  :: etot
      complex(dp)          :: Omattot(nbspx, nbspx)
      logical               :: lgam

      !
      ! local variables for cg routine
      !
      integer     :: nbnd1, nbnd2
      integer     :: isp
      logical     :: ldotest
      real(dp)    :: dtmp
      real(dp)    :: ene0, ene1, enesti, enever, dene0
      real(dp)    :: passo, passov, passof, passomax, spasso
      real(dp)    :: vsicah2sum, vsicah2sum_prev
      integer     :: nidx1, nidx2
      real(dp)    :: dPI, dalpha, dmaxeig, deigrms
      real(dp)    :: pinksumprev, passoprod
      !
      complex(dp), allocatable :: Omat1tot(:, :), Omat2tot(:, :)
      real(dp), allocatable :: Heigbig(:)
      complex(dp), allocatable :: Umatbig(:, :)
      complex(dp), allocatable :: wfc_ctmp(:, :), wfc_ctmp2(:, :)
      complex(dp), allocatable :: gi(:, :), hi(:, :)
      !
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)
      complex(dp), allocatable :: vsicah(:, :)
      real(dp), allocatable :: vsic1_realspace(:, :), vsic2_realspace(:, :)
      complex(dp), allocatable :: vsic_reciprocal1(:, :), vsic_reciprocal2(:, :)
      type(twin_matrix)       :: bec1, bec2
      real(dp), allocatable :: pink1(:), pink2(:)

      !
      ! main body
      !
      CALL start_clock('nk_rot_emin')
      !
      !
      pinksumprev = 1.d8
      dPI = 2.0_DP*asin(1.0_DP)
      passoprod = (0.3d0/dPI)*dPI

      !
      ! local workspace
      !
      allocate (Omat1tot(nbspx, nbspx), Omat2tot(nbspx, nbspx))
      allocate (Umatbig(nbspx, nbspx))
      allocate (Heigbig(nbspx))
      allocate (wfc_ctmp(ngw, nbspx), wfc_ctmp2(ngw, nbspx))
      allocate (hi(nbsp, nbsp))
      allocate (gi(nbsp, nbsp))
      allocate (pink1(nbspx), pink2(nbspx))
      allocate (vsic1_realspace(nnrx, nbspx), vsic2_realspace(nnrx, nbspx))
      allocate (vsic_reciprocal1(ngm, nbspx), vsic_reciprocal2(ngm, nbspx))
      call init_twin(bec1, lgam)
      call allocate_twin(bec1, nkb, nbsp, lgam)
      call init_twin(bec2, lgam)
      call allocate_twin(bec2, nkb, nbsp, lgam)
      !
      Umatbig(:, :) = (0.d0, 0.d0)
      Heigbig(:) = 0.d0
      deigrms = 0.d0
      hi(:, :) = 0.d0
      gi(:, :) = 0.d0

      Omattot(:, :) = 0.d0
      do nbnd1 = 1, nbspx
         Omattot(nbnd1, nbnd1) = CMPLX(1.d0, 0.d0)
      end do

      ninner = 0
      ldotest = .false.

      if (ionode) write (stdout, "(14x,'# iter',6x,'etot',17x,'esic',&
                                  & 17x,'deigrms')")

      !
      ! main loop
      !
      inner_loop: &
         do while (.true.)

         call start_clock("nk_innerloop")
         !
         ninner = ninner + 1

         if (ninner > innerloop_nmax) then
            !
#ifdef __DEBUG
            if (ionode) write (1031, *) '# innerloop_nmax reached.'
            if (ionode) write (1031, *)
#endif
            if (ionode) then
               write (stdout, "(14x,'# innerloop_nmax reached.',/)")
            end if
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if

#ifdef __DEBUG
         !
         ! call nksic_printoverlap(ninner,nouter)

!        if(mod(ninner,10).eq.1.or.ninner.le.5) ldotest=.true.
         if (ninner .eq. 31 .or. ninner .eq. 61 .or. ninner .eq. 91) ldotest = .true.
!        if(ninner.le.10.and.nouter.eq.1) ldotest=.true.
!         ldotest=.true.
!        if(ninner.ge.25) ldotest=.true.
         ! Now do the test
         if (ldotest) then
!          dtmp = 1.0d0*3.141592d0
            dtmp = 4.d0*3.141592d0
!          call nksic_rot_test(dtmp,201,nouter,ninner,etot)
            ldotest = .false.
         end if
#endif

         !
         !print out ESIC part & other total energy
         !
         ene0 = sum(pink(1:nbsp))

         !
         ! test convergence
         !
         if (abs(ene0 - pinksumprev) < esic_conv_thr) then
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, "(a,/)") '# inner-loop converged.'
               write (1031, "(a,/)") '# inner-loop converged.'
            end if
#endif
            if (ionode) write (stdout, "(14x,'# innerloop converged',/)")
            !
            call stop_clock("nk_innerloop")
            exit inner_loop
            !
         end if
         !
         pinksumprev = ene0

         !
         ! This part calculates the anti-hermitian part of the Hamiltonian vsicah
         ! and see whether a convergence has been achieved
         !
         ! For this run, we obtain the gradient
         !
         vsicah2sum = 0.0d0
         !
         do isp = 1, nspin
            !
            allocate (vsicah(nupdwn(isp), nupdwn(isp)))
            !
            call nksic_getvsicah_new2(isp, vsicah, dtmp, lgam)
            !
            gi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
               iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = vsicah(:, :)
            !
            vsicah2sum = vsicah2sum + dtmp
            !
            deallocate (vsicah)
            !
         end do
         !
         if (ninner /= 1) dtmp = vsicah2sum/vsicah2sum_prev
         !
         if (ninner <= innerloop_cg_nsd .or. &
             mod(ninner, innerloop_cg_nreset) == 0) then
            !
            hi(:, :) = gi(:, :)
         else
            hi(:, :) = gi(:, :) + dtmp*hi(:, :)
         end if
         !
         spin_loop: &
            do isp = 1, nspin
            !
            allocate (vsicah(nupdwn(isp), nupdwn(isp)))
            allocate (Umat(nupdwn(isp), nupdwn(isp)))
            allocate (Heig(nupdwn(isp)))

            vsicah(:, :) = hi(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                              iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))

            call nksic_getHeigU(isp, vsicah, Heig, Umat)
            !
            deigrms = deigrms + sum(Heig(:)**2)
            !
            Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), &
                    iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Umat(:, :)
            Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Heig(:)
            !
            deallocate (vsicah)
            deallocate (Umat)
            deallocate (Heig)
            !
         end do spin_loop

         ! how severe the transform is
         deigrms = sqrt(deigrms/nbsp)
#ifdef __DEBUG
         if (ionode) write (1031, '(2I10,3F24.13)') ninner, nouter, etot, ene0, deigrms
#endif
         if (ionode) write (stdout, '(10x,2i5,3F21.13)') ninner, nouter, etot, ene0, deigrms
         !
         !
         dmaxeig = max(abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
         !
         do isp = 2, nspin
            dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp))))
            dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
         end do
         !
         passomax = passoprod/dmaxeig
         !
         if (ninner == 1) then
            passof = passomax
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passof set to passomax'
#endif
         end if

!$$$$        if(passof .gt. passomax*2.d0) then
!$$$$          passof = passomax*2.d0
!$$$$          if(ionode) write(1031,*) '# passof > twice passomax'
!$$$$        endif

!        if(ionode) then
!          write(1037,*)'# deigrms = ',deigrms
!          write(1037,*)'# vsicah2sum = ',vsicah2sum
!          if(ninner.ne.1) write(1037,*)'# vsicah2sum/vsicah2sum_prev = ',dtmp
!        endif

         vsicah2sum_prev = vsicah2sum
         !
         dene0 = 0.d0
         !
         do isp = 1, nspin
            !
            do nbnd1 = 1, nupdwn(isp)
            do nbnd2 = 1, nupdwn(isp)
               !
               nidx1 = nbnd1 - 1 + iupdwn(isp)
               nidx2 = nbnd2 - 1 + iupdwn(isp)
               IF (nidx1 .eq. nidx2) THEN
                  dene0 = dene0 - DBLE(CONJG(gi(nidx1, nidx2))*hi(nidx1, nidx2))
               ELSE
                  dene0 = dene0 - 0.5d0*DBLE(CONJG(gi(nidx1, nidx2))*hi(nidx1, nidx2))
               END IF
               !
            end do
            end do
            !
         end do

         !$$
         !$$ dene0 = dene0 * 2.d0/nspin
         !
         ! Be careful, the following is correct because A_ji = - A_ij, i.e., the number of
         ! linearly independent variables is half the number of total variables!
         !
         dene0 = dene0*2.d0/nspin
         !
         spasso = 1.d0
         if (dene0 > 0.d0) spasso = -1.d0
         !
         dalpha = spasso*passof
         !
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp, &
                               Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, ene1, lgam)
         call minparabola(ene0, spasso*dene0, ene1, passof, passo, enesti)

         !
         ! We neglect this step for paper writing purposes
         !
         if (passo > passomax) then
            passo = passomax
#ifdef __DEBUG
            if (ionode) write (1031, *) '# passo > passomax'
#endif
            !
         end if

         passov = passof
         passof = 2.d0*passo

         dalpha = spasso*passo
         !
!$$ The following line is for dene0 test
!        if(ninner.ge.15) dalpha = spasso*passo*0.00001
!$$
         call nksic_getOmattot(dalpha, Heigbig, Umatbig, c0, wfc_ctmp2, &
                               Omat2tot, bec2, vsic2_realspace, vsic_reciprocal2, pink2, enever, lgam)

#ifdef __DEBUG
         if (ionode) then
            !
            write (1037, *) ninner, nouter
            write (1037, '("ene0,ene1,enesti,enever")')
            write (1037, '(a3,4f20.10)') 'CG1', ene0, ene1, enesti, enever
            write (1037, '("spasso,passov,passo,passomax,dene0,&
                         & (enever-ene0)/passo/dene0")')
            write (1037, '(a3,4f12.7,e20.10,f12.7)') &
               'CG2', spasso, passov, passo, passomax, dene0, (enever - ene0)/passo/dene0
            write (1037, *)
            !
         end if
#endif

         if (ene0 < ene1 .and. ene0 < enever) then
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# ene0<ene1 and ene0<enever, exit",/)')
               write (1031, '("# innerloop NOT converged, exit",/)')
            end if
#endif
            !
            ninner = ninner + 1
            call stop_clock("nk_innerloop")
            !
            exit
            !
         end if

         if (ene1 >= enever) then
            !
            pink(:) = pink2(:)
            vsic_realspace(:, :) = vsic2_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal2(:, :)
            c0(:, :) = wfc_ctmp2(:, :)
            call copy_twin(bec, bec2)
!             bec%rvec(:,:)  = bec2(:,:)
            Omattot = MATMUL(Omattot, Omat2tot)
            !
         else
            !
            pink(:) = pink1(:)
            vsic_realspace(:, :) = vsic1_realspace(:, :)
            vsic_reciprocal(:, :) = vsic_reciprocal1(:, :)
            c0(:, :) = wfc_ctmp(:, :)
            call copy_twin(bec, bec1)
            Omattot = MATMUL(Omattot, Omat1tot)
            !
#ifdef __DEBUG
            if (ionode) then
               write (1037, '("# It happened that ene1 < enever!!")')
               write (1037, *)
            end if
#endif
            !
! =======
!           pink(:) = pink1(:)
!           vsic_realspace(:,:) = vsic1_realspace(:,:)
!           c0(:,:) = wfn_ctmp(:,:)
!           bec%rvec(:,:) = bec1(:,:)
!           Omattot = MATMUL(Omattot,Omat1tot)
!           if(ionode) then
!             write(1037,'("# It happened that ene1 < enever!!")')
!             write(1037,*)
!           endif
! 1.28.2.14
         end if
         !
         call stop_clock("nk_innerloop")
         !
      end do inner_loop

      !
      ! Wavefunction cm rotation according to Omattot
      ! We need this because outer loop could be damped dynamics.
      !
      if (.not. tcg) then
         !
         if (ninner >= 2) then
            !
            wfc_ctmp(:, :) = (0.d0, 0.d0)
            !
            do nbnd1 = 1, nbspx
            do nbnd2 = 1, nbspx
               wfc_ctmp(:, nbnd1) = wfc_ctmp(:, nbnd1) + cm(:, nbnd2)*Omattot(nbnd2, nbnd1) !warning:giovanni CONJUGATE?
               ! XXX (we can think to use a blas, here, and split over spins)
            end do
            end do
            !
            cm(:, 1:nbspx) = wfc_ctmp(:, 1:nbspx)
            !
         end if
         !
      end if

      !
      ! clean local workspace
      !
      deallocate (Omat1tot, Omat2tot)
      deallocate (Umatbig)
      deallocate (Heigbig)
      deallocate (wfc_ctmp, wfc_ctmp2)
      deallocate (hi)
      deallocate (gi)
      deallocate (pink1, pink2)
      deallocate (vsic1_realspace, vsic2_realspace)
      deallocate (vsic_reciprocal1, vsic_reciprocal2)
      call deallocate_twin(bec1)
      call deallocate_twin(bec2)
!       deallocate( bec1, bec2 )

      CALL stop_clock('nk_rot_emin')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rot_emin_cg_descla
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_getOmattot_new(nbsp, nbspx, nudx, nspin, ispin, &
                                   iupdwn, nupdwn, wfc_centers, wfc_spreads, &
                                   dalpha, Heigbig, Umatbig, wfc0, &
                                   wfc1, Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, ene1, lgam, is_empty)!warning:giovanni bec1 here needs to be a twin!
!---------------------------------------------------------------
!
! ... This routine rotates the wavefunction wfc0 into wfc1 according to
!     the force matrix (Heigbig, Umatbig) and the step of size dalpha.
!     Other quantities such as bec, vsic_realspace, pink are also calculated for wfc1.

      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use ions_base, only: nsp
      use uspp, only: becsum
      use cp_main_variables, only: eigr, rhor
      use nksic, only: deeq_sic, wtot_realspace, wtot_reciprocal, fsic
      use control_flags, only: gamma_only, do_wf_cmplx
      use twin_types
      use electrons_module, only: icompute_spread
      use core, only: rhoc
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in) :: nbsp, nbspx, nudx, nspin
      integer, intent(in) :: ispin(nbspx), nupdwn(nspin), &
                             iupdwn(nspin)
      real(dp), intent(in) :: dalpha
      complex(dp), intent(in) :: Umatbig(nbspx, nbspx)
      real(dp), intent(in) :: Heigbig(nbspx), wfc_centers(4, nudx, nspin), &
                              wfc_spreads(nudx, nspin, 2)
      complex(dp), intent(in) :: wfc0(ngw, nbspx)
      !
      complex(dp)                    :: wfc1(ngw, nbspx)
      complex(dp)                       :: Omat1tot(nbspx, nbspx)
      type(twin_matrix)     :: bec1 !(nkb,nbsp) !modified:giovanni
      real(dp)                       :: vsic1_realspace(nnrx, nbspx)
      complex(dp)                    :: vsic_reciprocal1(ngm, nbspx)
      real(dp)                       :: pink1(nbspx)
      real(dp)                       :: ene1
      logical :: lgam, is_empty

      !
      ! local variables for cg routine
      !
      integer    :: isp, nbnd1
      real(dp)   :: dmaxeig
      complex(dp), allocatable :: Omat1(:, :)
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)

      !
      call start_clock("nk_getOmattot")
      !

!       call init_twin(bec1,lgam)
!       call allocate_twin(bec1,nkb,nbsp, lgam)

      Omat1tot(:, :) = 0.d0
      do nbnd1 = 1, nbspx
         Omat1tot(nbnd1, nbnd1) = 1.d0
      end do

      wfc1(:, :) = CMPLX(0.d0, 0.d0)

      dmaxeig = max(abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
      do isp = 2, nspin
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp))))
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
      end do

      spin_loop: &
         do isp = 1, nspin
         !
         IF (nupdwn(isp) .gt. 0) THEN
            !
            allocate (Umat(nupdwn(isp), nupdwn(isp)))
            allocate (Heig(nupdwn(isp)))
            allocate (Omat1(nupdwn(isp), nupdwn(isp)))

            Umat(:, :) = Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))
            Heig(:) = Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))

            call nksic_getOmat1(isp, Heig, Umat, dalpha, Omat1, lgam)

!$$ Wavefunction wfc0 is rotated into wfc0 using Omat1
            call nksic_rotwfn(isp, Omat1, wfc0, wfc1)

! Assigning the rotation matrix for a specific spin isp
            Omat1tot(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Omat1(:, :)

            deallocate (Umat)
            deallocate (Heig)
            deallocate (Omat1)
            !
         ELSE
            Omat1tot(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 1.d0
         END IF
         !
      end do spin_loop

      !
      ! recalculate bec & vsic_realspace according to the new wavefunction
      !
      call calbec(1, nsp, eigr, wfc1, bec1)

      vsic1_realspace(:, :) = 0.d0
      vsic_reciprocal1(:, :) = 0.d0
      pink1(:) = 0.d0
      !
      call nksic_potential(nbsp, nbspx, wfc1, fsic, bec1, becsum, deeq_sic, &
                      ispin, iupdwn, nupdwn, rhor, rhoc, wtot_realspace, wtot_reciprocal, vsic1_realspace, vsic_reciprocal1, pink1, nudx, wfc_centers, &
                           wfc_spreads, icompute_spread, is_empty)
      !
      ene1 = sum(pink1(:))

!       call deallocate_twin(bec1)

      !
      call stop_clock("nk_getOmattot")
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getOmattot_new
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_getOmattot(dalpha, Heigbig, Umatbig, wfc0, wfc1, Omat1tot, bec1, vsic1_realspace, vsic_reciprocal1, pink1, ene1, lgam)!warning:giovanni bec1 here needs to be a twin!
!---------------------------------------------------------------
!
! ... This routine rotates the wavefunction wfc0 into wfc1 according to
!     the force matrix (Heigbig, Umatbig) and the step of size dalpha.
!     Other quantities such as bec, vsic_realspace, pink are also calculated for wfc1.
!

      use kinds, only: dp
      use grid_dimensions, only: nnrx
      use gvecp, only: ngm
      use gvecw, only: ngw
      use electrons_base, only: nbsp, nbspx, nspin, ispin, &
                                iupdwn, nupdwn, nudx
      use ions_base, only: nsp
      use uspp, only: becsum
      use cp_main_variables, only: eigr, rhor
      use nksic, only: deeq_sic, wtot_realspace, wtot_reciprocal, fsic, do_wxd, &
                       valpsi, odd_alpha
      use control_flags, only: gamma_only, do_wf_cmplx
      use twin_types
      use electrons_module, only: wfc_centers, wfc_spreads, &
                                  icompute_spread
      use core, only: rhoc
      use input_parameters, only: odd_nkscalfact
      !
      implicit none
      !
      ! in/out vars
      !
      real(dp), intent(in) :: dalpha
      complex(dp), intent(in) :: Umatbig(nbspx, nbspx)
      real(dp), intent(in) :: Heigbig(nbspx)
      complex(dp), intent(in) :: wfc0(ngw, nbspx)
      !
      complex(dp)                    :: wfc1(ngw, nbspx)
      complex(dp)                       :: Omat1tot(nbspx, nbspx)
      type(twin_matrix)     :: bec1 !(nkb,nbsp) !modified:giovanni
      real(dp)                       :: vsic1_realspace(nnrx, nbspx)
      complex(dp)                    :: vsic_reciprocal1(ngm, nbspx)
      real(dp)                       :: pink1(nbspx)
      real(dp)                       :: ene1
      logical :: lgam

      !
      ! local variables for cg routine
      !
      integer    :: isp, nbnd1
      real(dp)   :: dmaxeig
      complex(dp), allocatable :: Omat1(:, :)
      complex(dp), allocatable :: Umat(:, :)
      real(dp), allocatable :: Heig(:)

      !
      call start_clock("nk_getOmattot")
      !

!       call init_twin(bec1,lgam)
!       call allocate_twin(bec1,nkb,nbsp, lgam)

      Omat1tot(:, :) = 0.d0
      do nbnd1 = 1, nbspx
         Omat1tot(nbnd1, nbnd1) = 1.d0
      end do

      wfc1(:, :) = CMPLX(0.d0, 0.d0)

      dmaxeig = max(abs(Heigbig(iupdwn(1))), abs(Heigbig(iupdwn(1) + nupdwn(1) - 1)))
      do isp = 2, nspin
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp))))
         dmaxeig = max(dmaxeig, abs(Heigbig(iupdwn(isp) + nupdwn(isp) - 1)))
      end do

      spin_loop: &
         do isp = 1, nspin
         !
         IF (nupdwn(isp) .gt. 0) THEN
            !
            allocate (Umat(nupdwn(isp), nupdwn(isp)))
            allocate (Heig(nupdwn(isp)))
            allocate (Omat1(nupdwn(isp), nupdwn(isp)))

            Umat(:, :) = Umatbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))
            Heig(:) = Heigbig(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp))

            call nksic_getOmat1(isp, Heig, Umat, dalpha, Omat1, lgam)

!$$ Wavefunction wfc0 is rotated into wfc0 using Omat1
            call nksic_rotwfn(isp, Omat1, wfc0, wfc1)

! Assigning the rotation matrix for a specific spin isp
            Omat1tot(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = Omat1(:, :)

            deallocate (Umat)
            deallocate (Heig)
            deallocate (Omat1)
            !
         ELSE
            Omat1tot(iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp), iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = 1.d0
         END IF
         !
      end do spin_loop
      !
      ! recalculate bec & vsic_realspace according to the new wavefunction
      !
      call calbec(1, nsp, eigr, wfc1, bec1)
      !
      if (odd_nkscalfact) then
         !
         valpsi(:, :) = (0.0_DP, 0.0_DP)
         odd_alpha(:) = 0.0_DP
         !
         call odd_alpha_routine(nbspx, .false.)
         !
      end if
      !
      vsic1_realspace(:, :) = 0.d0
      vsic_reciprocal1(:, :) = 0.d0
      pink1(:) = 0.d0
      !
      !
      call nksic_potential(nbsp, nbspx, wfc1, fsic, bec1, becsum, deeq_sic, &
              ispin, iupdwn, nupdwn, rhor, rhoc, wtot_realspace, wtot_reciprocal, vsic1_realspace, vsic_reciprocal1, do_wxd, pink1, nudx, wfc_centers, &
                           wfc_spreads, icompute_spread, .false.)
      !
      ene1 = sum(pink1(:))

!       call deallocate_twin(bec1)

      !
      call stop_clock("nk_getOmattot")
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getOmattot
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_rotwfn(isp, Omat1, wfc1, wfc2)
!-----------------------------------------------------------------------
!
! ... Simple rotation of wfc1 into wfc2 by Omat1.
!     wfc2(n) = sum_m wfc1(m) Omat1(m,n)
!
      use electrons_base, only: iupdwn, nupdwn, nbspx
      use gvecw, only: ngw
      use kinds, only: dp
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp
      complex(dp), intent(in)  :: Omat1(nupdwn(isp), nupdwn(isp))
      complex(dp), intent(in)  :: wfc1(ngw, nbspx)
      complex(dp)              :: wfc2(ngw, nbspx)

      !
      ! local variables for cg routine
      !
      integer                  :: nbnd1, nbnd2

      CALL start_clock('nk_rotwfn')
      !
      wfc2(:, iupdwn(isp):iupdwn(isp) - 1 + nupdwn(isp)) = CMPLX(0.d0, 0.d0)

      !
      ! a blas could be used here XXX
      !
      do nbnd1 = 1, nupdwn(isp)
      do nbnd2 = 1, nupdwn(isp)
         !
         wfc2(:, iupdwn(isp) - 1 + nbnd1) = wfc2(:, iupdwn(isp) - 1 + nbnd1) &
                                            + wfc1(:, iupdwn(isp) - 1 + nbnd2)*Omat1(nbnd2, nbnd1)
         !
      end do
      end do

      CALL stop_clock('nk_rotwfn')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_rotwfn
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getHeigU_new(nspin, isp, nupdwn, vsicah, Heig, Umat)
!-----------------------------------------------------------------------
!
! ... solves the eigenvalues (Heig) and eigenvectors (Umat) of the force
!     matrix vsicah.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use mp, only: mp_bcast
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp, nspin, nupdwn(nspin)
      real(dp)     :: Heig(nupdwn(isp))
      complex(dp)  :: Umat(nupdwn(isp), nupdwn(isp))
      complex(dp)     :: vsicah(nupdwn(isp), nupdwn(isp))

      !
      ! local variables
      !
      complex(dp)              :: Hmat(nupdwn(isp), nupdwn(isp))
      complex(dp)              :: ci

      ci = CMPLX(0.d0, 1.d0)

!$$ Now this part diagonalizes Hmat = iWmat
      Hmat(:, :) = ci*vsicah(:, :)
!$$ diagonalize Hmat
!      if(ionode) then
      CALL zdiag(nupdwn(isp), nupdwn(isp), Hmat(1, 1), Heig(1), Umat(1, 1), 1)
!      endif

!      CALL mp_bcast(Umat, ionode_id, intra_image_comm)
!      CALL mp_bcast(Heig, ionode_id, intra_image_comm)

      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getHeigU_new
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getHeigU(isp, vsicah, Heig, Umat)
!-----------------------------------------------------------------------
!
! ... solves the eigenvalues (Heig) and eigenvectors (Umat) of the force
!     matrix vsicah.
!     (Ultrasoft pseudopotential case is not implemented.)
!
      use kinds, only: dp
      use mp, only: mp_bcast
      use electrons_base, only: nupdwn
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp
      real(dp)     :: Heig(nupdwn(isp))
      complex(dp)  :: Umat(nupdwn(isp), nupdwn(isp))
      complex(dp)     :: vsicah(nupdwn(isp), nupdwn(isp))

      !
      ! local variables
      !
      complex(dp)              :: Hmat(nupdwn(isp), nupdwn(isp))
      complex(dp)              :: ci

      ci = CMPLX(0.d0, 1.d0)

!$$ Now this part diagonalizes Hmat = iWmat
      Hmat(:, :) = ci*vsicah(:, :)
!$$ diagonalize Hmat
!      if(ionode) then
      CALL zdiag(nupdwn(isp), nupdwn(isp), Hmat(1, 1), Heig(1), Umat(1, 1), 1)
!      endif

!      CALL mp_bcast(Umat, ionode_id, intra_image_comm)
!      CALL mp_bcast(Heig, ionode_id, intra_image_comm)

      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getHeigU
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_printoverlap(ninner, nouter)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!
      use kinds, only: dp
      use grid_dimensions, only: nr1x, nr2x, nr3x, nnrx
      use gvecw, only: ngw
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use io_global, only: ionode
      use electrons_base, only: nbspx
      use cp_interfaces, only: invfft
      use fft_base, only: dfftp
      use nksic, only: vsic_realspace
      use wavefunctions_module, only: c0
      !
      implicit none
      !
      ! in/out vars
      !
      integer   :: ninner, nouter
      real(dp)  :: overlap(nbspx, nbspx), vsicah(nbspx, nbspx)

      !
      ! local variables
      !
      complex(dp)              :: psi1(nnrx), psi2(nnrx)
      real(dp)                 :: overlaptmp, vsicahtmp
      integer                  :: i, nbnd1, nbnd2
      real(dp)                 :: dwfnnorm

      dwfnnorm = 1.0/(DBLE(nr1x)*DBLE(nr2x)*DBLE(nr3x))

      vsicah(:, :) = 0.d0
      overlap(:, :) = 0.d0

      do nbnd1 = 1, nbspx
         CALL c2psi(psi1, nnrx, c0(:, nbnd1), (0.d0, 0.d0), ngw, 1)
         CALL invfft('Dense', psi1, dfftp)

         do nbnd2 = 1, nbspx
            if (nbnd2 .lt. nbnd1) then
               vsicahtmp = -vsicah(nbnd2, nbnd1)
               overlaptmp = overlap(nbnd2, nbnd1)
            else
               CALL c2psi(psi2, nnrx, c0(:, nbnd2), (0.d0, 0.d0), ngw, 1)
               CALL invfft('Dense', psi2, dfftp)

               vsicahtmp = 0.d0
               overlaptmp = 0.d0

               do i = 1, nnrx
!$$ Imposing Pederson condition
                  vsicahtmp = vsicahtmp &
                              + 2.d0*DBLE(CONJG(psi1(i))*(vsic_realspace(i, nbnd2) &
                                                          - vsic_realspace(i, nbnd1))*psi2(i))*dwfnnorm
!$$ The following two lines give exactly the same results: checked
                  overlaptmp = overlaptmp + DBLE(CONJG(psi1(i))*psi2(i))*dwfnnorm
!              overlaptmp = overlaptmp + dble(psi1(i)) * dble(psi2(i)) * dwfnnorm
               end do

               CALL mp_sum(vsicahtmp, intra_image_comm)
               CALL mp_sum(overlaptmp, intra_image_comm)
            end if ! if(nbnd2.lt.nbnd1)

            vsicah(nbnd1, nbnd2) = vsicahtmp
            overlap(nbnd1, nbnd2) = overlaptmp

         end do ! nbnd2=1,nbspx

      end do ! nbspx

      if (ionode) then
         write (1021, *) ninner, nouter
         write (1022, *) ninner, nouter
         do nbnd1 = 1, nbspx
            write (1021, '(100F12.7)') (overlap(nbnd1, nbnd2), nbnd2=1, nbspx)
            write (1022, '(100F12.7)') (vsicah(nbnd1, nbnd2), nbnd2=1, nbspx)
         end do
         write (1021, *)
         write (1022, *)
      end if

      return
      !
!---------------------------------------------------------------
   end subroutine nksic_printoverlap
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getvsicah(isp, vsicah, vsicah2sum)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!
      use kinds, only: dp
      use grid_dimensions, only: nr1x, nr2x, nr3x, nnrx
      use gvecw, only: ngw
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use electrons_base, only: nspin, iupdwn, nupdwn
      use cp_interfaces, only: invfft
      use fft_base, only: dfftp
      use nksic, only: vsic_realspace, fsic
      use wavefunctions_module, only: c0
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp
      real(dp)                 :: vsicah(nupdwn(isp), nupdwn(isp))
      real(dp)                 :: vsicah2sum

      !
      ! local variables
      !
      complex(dp)     :: psi1(nnrx), psi2(nnrx)
      real(dp)        :: vsicahtmp, cost
      real(dp)        :: dwfnnorm
      integer         :: nbnd1, nbnd2
      integer         :: i, j1, j2

      CALL start_clock('nk_get_vsicah')
      !
      dwfnnorm = 1.0d0/(DBLE(nr1x)*DBLE(nr2x)*DBLE(nr3x))
      cost = 2.0d0*DBLE(nspin)*0.5d0*dwfnnorm
      !
      vsicah(:, :) = 0.d0
      vsicah2sum = 0.d0

      !
      ! Imposing Pederson condition
      !
      do nbnd1 = 1, nupdwn(isp)
         !
         j1 = iupdwn(isp) - 1 + nbnd1
         !
         CALL c2psi(psi1, nnrx, c0(:, j1), (0.d0, 0.d0), ngw, 1)
         CALL invfft('Dense', psi1, dfftp)

         do nbnd2 = 1, nbnd1 - 1
            !
            j2 = iupdwn(isp) - 1 + nbnd2
            !
            CALL c2psi(psi2, nnrx, c0(:, j2), (0.0d0, 0.0d0), ngw, 1)
            CALL invfft('Dense', psi2, dfftp)
            !
            vsicahtmp = 0.d0
            !
            do i = 1, nnrx
               !
               vsicahtmp = vsicahtmp + &
                           DBLE(CONJG(psi1(i))*psi2(i) &
                                *(vsic_realspace(i, j2)*fsic(j2) &
                                  - vsic_realspace(i, j1)*fsic(j1)))
               !
            end do
            vsicahtmp = vsicahtmp*cost
            !
            vsicah(nbnd1, nbnd2) = vsicahtmp
            vsicah(nbnd2, nbnd1) = -vsicahtmp
            !
         end do
         !
      end do
      !
      call mp_sum(vsicah, intra_image_comm)
      !
      vsicah2sum = 0.0d0
      do nbnd1 = 1, nupdwn(isp)
      do nbnd2 = 1, nbnd1 - 1
         vsicah2sum = vsicah2sum + 2.0d0*vsicah(nbnd2, nbnd1)*vsicah(nbnd2, nbnd1)
      end do
      end do
      !
      call stop_clock('nk_get_vsicah')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getvsicah
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getvsicah_new1(isp, vsicah, vsicah2sum)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!     Exploit fft of wfc pairs.
!
      use kinds, only: dp
      use grid_dimensions, only: nr1x, nr2x, nr3x, nnrx
      use gvecw, only: ngw
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use electrons_base, only: nspin, iupdwn, nupdwn
      use cp_interfaces, only: invfft
      use fft_base, only: dfftp
      use nksic, only: vsic_realspace, fsic
      use wavefunctions_module, only: c0
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp
      real(dp)                 :: vsicah(nupdwn(isp), nupdwn(isp))
      real(dp)                 :: vsicah2sum

      !
      ! local variables
      !
      real(dp)        :: vsicahtmp, cost
      real(dp)        :: dwfnnorm
      integer         :: nbnd1, nbnd2
      integer         :: i, j1, jj1, j2, jj2
      !
      complex(dp), allocatable :: psi1(:), psi2(:)
      real(dp), allocatable :: wfc1(:, :), wfc2(:, :)

      CALL start_clock('nk_get_vsicah')
      !
      dwfnnorm = 1.0d0/(DBLE(nr1x)*DBLE(nr2x)*DBLE(nr3x))
      cost = 2.0d0*DBLE(nspin)*0.5d0*dwfnnorm
      !
      allocate (wfc1(nnrx, 2))
      allocate (wfc2(nnrx, 2))
      allocate (psi1(nnrx))
      allocate (psi2(nnrx))

      !
      ! Imposing Pederson condition
      !
      vsicah(:, :) = 0.d0
      !
      do nbnd1 = 1, nupdwn(isp), 2
         !
         j1 = iupdwn(isp) - 1 + nbnd1
         !
         CALL c2psi(psi1, nnrx, c0(:, j1), c0(:, j1 + 1), ngw, 2)
         CALL invfft('Dense', psi1, dfftp)
         !
         wfc1(:, 1) = DBLE(psi1(:))
         wfc1(:, 2) = AIMAG(psi1(:))
         !
         do jj1 = 1, 2
            !
            if (nbnd1 + jj1 - 1 > nupdwn(isp)) cycle
            !
            !
            do nbnd2 = 1, nbnd1 - 1 + jj1 - 1, 2
               !
               j2 = iupdwn(isp) - 1 + nbnd2
               !
               CALL c2psi(psi2, nnrx, c0(:, j2), c0(:, j2 + 1), ngw, 2)
               CALL invfft('Dense', psi2, dfftp)
               !
               wfc2(:, 1) = DBLE(psi2(:))
               wfc2(:, 2) = AIMAG(psi2(:))
               !
               do jj2 = 1, 2
                  !
                  if (nbnd2 + jj2 - 1 > nbnd1 - 1 + jj1 - 1) cycle
                  !
                  vsicahtmp = 0.d0
                  !
                  do i = 1, nnrx
                     !
                     vsicahtmp = vsicahtmp + &
                                 cost*DBLE(wfc1(i, jj1)*wfc2(i, jj2) &
                                           *(vsic_realspace(i, j2 + jj2 - 1)*fsic(j2 + jj2 - 1) &
                                             - vsic_realspace(i, j1 + jj1 - 1)*fsic(j1 + jj1 - 1)))
                     !
                  end do
                  !
                  vsicah(nbnd1 + jj1 - 1, nbnd2 + jj2 - 1) = vsicahtmp
                  vsicah(nbnd2 + jj2 - 1, nbnd1 + jj1 - 1) = -vsicahtmp
                  !
               end do
            end do
            !
         end do
      end do
      !
      call mp_sum(vsicah, intra_image_comm)
      !
      vsicah2sum = 0.0d0
      !
      do nbnd1 = 1, nupdwn(isp)
      do nbnd2 = 1, nbnd1 - 1
         vsicah2sum = vsicah2sum + 2.0d0*vsicah(nbnd2, nbnd1)*vsicah(nbnd2, nbnd1)
      end do
      end do
      !
      !
      deallocate (wfc1, wfc2)
      deallocate (psi1, psi2)
      !
      call stop_clock('nk_get_vsicah')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getvsicah_new1
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getvsicah_new3(ngw, nbsp, nbspx, nspin, c0, bec, &
                                   isp, nupdwn, iupdwn, vsicah, vsicah2sum, lgam)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!     makes use of nksic_eforce to compute   h_i | phi_i >
!     and then computes   < phi_j | h_i | phi_i >  in reciprocal space.
!
      use kinds, only: dp
      use reciprocal_vectors, only: gstart
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use cp_interfaces, only: invfft
      use nksic, only: vsic_realspace, vsic_reciprocal, vsicpsi, &
                       deeq_sic  ! to be passed directly
      use twin_types
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp, nspin, ngw, nbsp, nbspx, &
                              nupdwn(nspin), iupdwn(nspin)
      complex(dp)              :: vsicah(nupdwn(isp), nupdwn(isp)), c0(ngw, nbsp)
      real(dp)                 :: vsicah2sum
      logical                  :: lgam
      type(twin_matrix)        :: bec

      !
      ! local variables
      !
      real(dp)        :: cost
      integer         :: nbnd1, nbnd2
      integer         :: j1, jj1, j2
      !
      !complex(dp), allocatable :: vsicpsi(:,:)
      complex(dp), allocatable :: hmat(:, :)

      CALL start_clock('nk_get_vsicah')
      !
      cost = DBLE(nspin)*2.0d0
      !
      !allocate( vsicpsi(npw,2) )
      allocate (hmat(nupdwn(isp), nupdwn(isp)))

      !
      ! compute < phi_j | Delta h_i | phi_i >
      !
      do nbnd1 = 1, nupdwn(isp), 2
         !
         ! NOTE: USPP not implemented
         !
         j1 = nbnd1 + iupdwn(isp) - 1
         CALL nksic_eforce(j1, nbsp, nbspx, vsic_realspace, vsic_reciprocal, &
                           deeq_sic, bec, ngw, c0(:, j1), c0(:, j1 + 1), vsicpsi, lgam)
         !
         do jj1 = 1, 2
            !
            if (nbnd1 + jj1 - 1 > nupdwn(isp)) cycle
            !
            do nbnd2 = 1, nupdwn(isp)
               !
               j2 = nbnd2 + iupdwn(isp) - 1
               IF (lgam) THEN
                  hmat(nbnd2, nbnd1 + jj1 - 1) = 2.d0*DBLE(DOT_PRODUCT(c0(:, j2), vsicpsi(:, jj1)))
                  !
                  if (gstart == 2) then
                     hmat(nbnd2, nbnd1 + jj1 - 1) = hmat(nbnd2, nbnd1 + jj1 - 1) - &
                                                    DBLE(c0(1, j2)*vsicpsi(1, jj1))
                  end if
               ELSE
                  hmat(nbnd2, nbnd1 + jj1 - 1) = DOT_PRODUCT(c0(:, j2), vsicpsi(:, jj1))
               END IF
               !
            end do
            !
         end do
      end do
      !
      call mp_sum(hmat, intra_image_comm)
      hmat = hmat*cost

      !
      ! Imposing Pederson condition
      !
      vsicah(:, :) = 0.d0
      vsicah2sum = 0.0d0
      !
      do nbnd1 = 1, nupdwn(isp)
      do nbnd2 = 1, nbnd1 - 1
         !
         IF (lgam) THEN
            vsicah(nbnd2, nbnd1) = DBLE(hmat(nbnd2, nbnd1) - CONJG(hmat(nbnd1, nbnd2)))
            vsicah(nbnd1, nbnd2) = DBLE(hmat(nbnd1, nbnd2) - CONJG(hmat(nbnd2, nbnd1)))
         ELSE
            vsicah(nbnd2, nbnd1) = hmat(nbnd2, nbnd1) - CONJG(hmat(nbnd1, nbnd2))
            vsicah(nbnd1, nbnd2) = hmat(nbnd1, nbnd2) - CONJG(hmat(nbnd2, nbnd1))
         END IF
         vsicah2sum = vsicah2sum + DBLE(CONJG(vsicah(nbnd2, nbnd1))*vsicah(nbnd2, nbnd1))
         !
      end do
      !IF(.not.lgam) THEN
      !  vsicah( nbnd1, nbnd1) = hmat(nbnd1,nbnd1) -CONJG(hmat(nbnd1,nbnd1))
      !  vsicah2sum =  vsicah2sum + 2.d0*DBLE(CONJG(vsicah(nbnd1,nbnd1))*vsicah(nbnd1,nbnd1))
      !ENDIF
      end do
      !
      deallocate (hmat)
      !
      call stop_clock('nk_get_vsicah')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getvsicah_new3
!---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getvsicah_new2(isp, vsicah, vsicah2sum, lgam)
!-----------------------------------------------------------------------
!
! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
!     makes use of nksic_eforce to compute   h_i | phi_i >
!     and then computes   < phi_j | h_i | phi_i >  in reciprocal space.
!
      use kinds, only: dp
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use electrons_base, only: nspin, iupdwn, nupdwn, nbsp, nbspx
      use cp_interfaces, only: invfft
      use nksic, only: vsic_realspace, vsic_reciprocal, vsicpsi, &
                       valpsi, deeq_sic  ! to be passed directly
      use wavefunctions_module, only: c0
      use cp_main_variables, only: bec  ! to be passed directly
      use input_parameters, only: odd_nkscalfact
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: isp
      complex(dp)                 :: vsicah(nupdwn(isp), nupdwn(isp))
      real(dp)                 :: vsicah2sum
      logical                  :: lgam

      !
      ! local variables
      !
      real(dp)        :: cost
      integer         :: nbnd1, nbnd2
      integer         :: j1, jj1, j2
      !
      !complex(dp), allocatable :: vsicpsi(:,:)
      complex(dp), allocatable :: hmat(:, :)

      CALL start_clock('nk_get_vsicah')
      !
      cost = DBLE(nspin)*2.0d0
      !
      !allocate( vsicpsi(npw,2) )
      allocate (hmat(nupdwn(isp), nupdwn(isp)))

      !
      ! compute < phi_j | Delta h_i | phi_i >
      !
      do nbnd1 = 1, nupdwn(isp), 2
         !
         ! NOTE: USPP not implemented
         !
         j1 = nbnd1 + iupdwn(isp) - 1
         CALL nksic_eforce(j1, nbsp, nbspx, vsic_realspace, vsic_reciprocal, &
                           deeq_sic, bec, ngw, c0(:, j1), c0(:, j1 + 1), vsicpsi, lgam)
         !
         do jj1 = 1, 2
            !
            if (nbnd1 + jj1 - 1 > nupdwn(isp)) cycle
            !
            do nbnd2 = 1, nupdwn(isp)
               !
               j2 = nbnd2 + iupdwn(isp) - 1
               !
               IF (odd_nkscalfact) THEN
                  !
                  vsicpsi(:, jj1) = vsicpsi(:, jj1) + valpsi(nbnd1 + jj1 - 1, :)
                  !
               END IF
               !
               IF (lgam) THEN
                  hmat(nbnd2, nbnd1 + jj1 - 1) = 2.d0*DBLE(DOT_PRODUCT(c0(:, j2), vsicpsi(:, jj1)))
                  !
                  if (gstart == 2) then
                     hmat(nbnd2, nbnd1 + jj1 - 1) = hmat(nbnd2, nbnd1 + jj1 - 1) - &
                                                    DBLE(c0(1, j2)*vsicpsi(1, jj1))
                  end if
               ELSE
                  hmat(nbnd2, nbnd1 + jj1 - 1) = DOT_PRODUCT(c0(:, j2), vsicpsi(:, jj1))
               END IF
               !
            end do
            !
         end do
      end do
      !
      call mp_sum(hmat, intra_image_comm)
      hmat = hmat*cost

      !
      ! Imposing Pederson condition
      !
      vsicah(:, :) = 0.d0
      vsicah2sum = 0.0d0
      !
      do nbnd1 = 1, nupdwn(isp)
      do nbnd2 = 1, nbnd1 - 1
         !
         IF (lgam) THEN
            vsicah(nbnd2, nbnd1) = DBLE(hmat(nbnd2, nbnd1) - CONJG(hmat(nbnd1, nbnd2)))
            vsicah(nbnd1, nbnd2) = DBLE(hmat(nbnd1, nbnd2) - CONJG(hmat(nbnd2, nbnd1)))
         ELSE
            vsicah(nbnd2, nbnd1) = hmat(nbnd2, nbnd1) - CONJG(hmat(nbnd1, nbnd2))
            vsicah(nbnd1, nbnd2) = hmat(nbnd1, nbnd2) - CONJG(hmat(nbnd2, nbnd1))
         END IF
         vsicah2sum = vsicah2sum + DBLE(CONJG(vsicah(nbnd2, nbnd1))*vsicah(nbnd2, nbnd1))
         !
      end do
      !IF(.not.lgam) THEN
      !  vsicah( nbnd1, nbnd1) = hmat(nbnd1,nbnd1) -CONJG(hmat(nbnd1,nbnd1))
      !  vsicah2sum =  vsicah2sum + 2.d0*DBLE(CONJG(vsicah(nbnd1,nbnd1))*vsicah(nbnd1,nbnd1))
      !ENDIF
      end do
      !
      deallocate (hmat)
      !
      call stop_clock('nk_get_vsicah')
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_getvsicah_new2
!---------------------------------------------------------------

! !-----------------------------------------------------------------------
!       subroutine nksic_getvsicah_twin( vsicah, vsicah2sum, nlam, descla, lgam)
! !-----------------------------------------------------------------------
! ! warning:giovanni IMPLEMENT without spin, call spin-by-spin initialize vsicah outside
! !IT IS JUST LIKE LAMBDA MATRIX, NEED NO FURTHER DESCLA INITIALIZATION!!.. DO AS
! ! IN ORTHO_GAMMA... PASS DESCLA MATRIX
! ! ... Calculates the anti-hermitian part of the SIC hamiltonian, vsicah.
! !     makes use of nksic_eforce to compute   h_i | phi_i >
! !     and then computes   < phi_j | h_i | phi_i >  in reciprocal space.
! !
!       use kinds,                      only : dp
!       use grid_dimensions,            only : nr1x, nr2x, nr3x, nnrx
!       use gvecw,                      only : ngw
!       use reciprocal_vectors,         only : gstart
!       USE mp,             ONLY: mp_sum,mp_bcast, mp_root_sum
!       use mp_global,                  only : intra_image_comm, leg_ortho
!       use electrons_base,             only : nspin, iupdwn, nupdwn, nbsp,nbspx
!       use cp_interfaces,              only : invfft
!       use fft_base,                   only : dfftp
!       use nksic,                      only : vsic_realspace, fsic, vsicpsi, &
!                                              deeq_sic  ! to be passed directly
!       use wavefunctions_module,       only : c0
!       use cp_main_variables,          only : bec  ! to be passed directly
!       use twin_types
! !       USE cp_main_variables,        ONLY : collect_lambda, distribute_lambda, descla, nrlx
!       USE descriptors,       ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
!                                    descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
!                                    nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , la_me_ , la_nrl_, nlax_
!       !
!       implicit none
!       !
!       ! in/out vars
!       !
! !       integer,     intent(in)  :: nspin
!       type(twin_matrix), dimension(nspin) :: vsicah!( nupdwn(isp),nupdwn(isp))
!       real(dp)                 :: vsicah2sum
!       logical                  :: lgam
!       INTEGER     :: descla( descla_siz_ )
!       INTEGER :: np_rot, me_rot, comm_rot, nrl
!       !
!       ! local variables
!       !
!       real(dp)        :: cost
!       integer         :: nbnd1,nbnd2,isp
!       integer         :: i, j1, jj1, j2, jj2, nss, istart, is
!       INTEGER :: np(2), coor_ip(2), ipr, ipc, nr, nc, ir, ic, ii, jj, root, j, nlam, nlax
!       INTEGER :: desc_ip( descla_siz_ )
!       LOGICAL :: la_proc
!       !
!       !complex(dp), allocatable :: vsicpsi(:,:)
!       real(dp), allocatable :: mtmp(:,:)
!       complex(dp),    allocatable :: h0c0(:,:), mtmp_c(:,:)
! !       type(twin_matrix) :: c0hc0(nspin)!modified:giovanni
!
!       CALL start_clock('nk_get_vsicah')
!
!       nlax    = descla( nlax_ )
!       la_proc = ( descla( lambda_node_ ) > 0 )
!       nlam    = 1
!       if ( la_proc ) nlam = nlax_
!       !
!       !
!       ! warning:giovanni:put a check on dimensions here?? (like in ortho_base/ortho)
!       ! this check should be on dimensionality of vsicah
!       !
!       cost     = dble( nspin ) * 2.0d0
!       !
!       !allocate( vsicpsi(npw,2) )
! !       allocate(c0hc0(nspin))
!       allocate(h0c0(ngw,nbspx))
!
! !       do is=1,nspin
! !         call init_twin(c0hc0(is),lgam)
! !         call allocate_twin(c0hc0(is),nlam,nlam,lgam)
! !       enddo
!
!       !
!       ! compute < phi_j | Delta h_i | phi_i >
!       !
! !
!       do j1 = 1, nbsp, 2
!           !
!           ! NOTE: USPP not implemented
!           !
!           CALL nksic_eforce( j1, nbsp, nbspx, vsic_realspace, &
!                              deeq_sic, bec, ngw, c0(:,j1), c0(:,j1+1), h0c0(:,j1:j1+1), lgam )
!           !
!       enddo
!
!       DO is= 1, nspin
!
!         nss= nupdwn( is )
!         istart= iupdwn( is )
!
!         np(1) = descla( la_npr_ , is )
!         np(2) = descla( la_npc_ , is )
!
!         DO ipc = 1, np(2)
!             DO ipr = 1, np(1)
!
!               coor_ip(1) = ipr - 1
!               coor_ip(2) = ipc - 1
!               CALL descla_init( desc_ip, descla( la_n_ , is ), descla( la_nx_ , is ), np, coor_ip, descla( la_comm_ , is ), 1 )
!
!               nr = desc_ip( nlar_ )
!               nc = desc_ip( nlac_ )
!               ir = desc_ip( ilar_ )
!               ic = desc_ip( ilac_ )
!
!               CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
!                                 desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )
!               !
!               root = root * leg_ortho
!
!               IF(.not.c0hc0(is)%iscmplx) THEN
!                 ALLOCATE( mtmp( nr, nc ) )
!                 mtmp = 0.0d0
!                 CALL DGEMM( 'T', 'N', nr, nc, 2*ngw, - 2.0d0, c0( 1, istart + ir - 1 ), 2*ngw, &
!                           h0c0( 1, istart + ic - 1 ), 2*ngw, 0.0d0, mtmp, nr )
!                 IF (gstart == 2) THEN
!                   DO jj = 1, nc
!                       DO ii = 1, nr
!                         i = ii + ir - 1
!                         j = jj + ic - 1
!                         mtmp(ii,jj) = mtmp(ii,jj) + DBLE( c0( 1, i + istart - 1 ) ) * DBLE( h0c0( 1, j + istart - 1 ) )
!                       END DO
!                   END DO
!                 END IF
!                 mtmp=mtmp*cost
!               ELSE
!                 ALLOCATE( mtmp_c( nr, nc ) )
!                 mtmp_c = CMPLX(0.0d0,0.d0)
!                 CALL ZGEMM( 'C', 'N', nr, nc, ngw, CMPLX(- 1.0d0,0.d0), c0( 1, istart + ir - 1 ),ngw, &
!                           h0c0( 1, istart + ic - 1 ), ngw, CMPLX(0.0d0,0.d0), mtmp_c, nr )
!               ENDIF
!               mtmp_c=mtmp_c*cost
!               IF(.not.c0hc0(is)%iscmplx) THEN
!                 CALL mp_root_sum( mtmp, vsicah(is)%rvec(1:nr,1:nc), root, intra_image_comm )
!                 DEALLOCATE( mtmp )
!               ELSE
!                 CALL mp_root_sum( mtmp_c, vsicah(is)%cvec(1:nr,1:nc), root, intra_image_comm )
!                 DEALLOCATE( mtmp_c )
!               ENDIF
! !                  IF( coor_ip(1) == descla( la_myr_ , is ) .AND. &
! !                      coor_ip(2) == descla( la_myc_ , is ) .AND. descla( lambda_node_ , is ) > 0 ) THEN
! !                     c0hc0(1:nr,1:nc,is) = mtmp
! !                  END IF
!             END DO
!         END DO
! !
! ! fill mtmp or mtmp_c with hermitian conjugate of vsicah
! ! and
! ! antisymmetrize vsicah
!         IF(lgam) THEN
!           allocate(mtmp(nlam,nlam))
!           mtmp=0.d0
!           CALL sqr_tr_cannon( nupdw(is), vsicah(is)%rvec, nlam, mtmp, nlam, descla )
!           DO i=1,nr
!               DO j=1,nc
!                 vsicah(is)%rvec(i,j) = vsicah(is)%rvec(i,j)-mtmp(i,j)
!               END DO
!           END DO
!          deallocate(mtmp)
!         ELSE
!           allocate(mtmp_c(nlam,nlam))
!           mtmp_c=0.d0
!           CALL sqr_tr_cannon( nupdw(is), vsicah(is)%cvec, nlam, mtmp_c, nlam, descla )
!           DO i=1,nr
!               DO j=1,nc
!                 vsicah(is)%cvec(i,j) = vsicah(is)%cvec(i,j)-mtmp(i,j)
!               END DO
!           END DO
!           deallocate(mtmp_c)
!         ENDIF
!
!       END DO
!
!       !
!       ! Imposing Pederson condition
!       !
!
! !       vsicah(:,:) = 0.d0
! !       vsicah2sum = 0.0d0
! !       !
! !       do nbnd1 = 1, nupdwn(isp)
! !       do nbnd2 = 1, nbnd1-1
! !           !
! !           IF(lgam) THEN
! !             vsicah( nbnd2, nbnd1) = DBLE(hmat(nbnd2,nbnd1) -CONJG(hmat(nbnd1,nbnd2)))
! !             vsicah( nbnd1, nbnd2) = DBLE(hmat(nbnd1,nbnd2) -CONJG(hmat(nbnd2,nbnd1)))
! !           ELSE
! !             vsicah( nbnd2, nbnd1) = hmat(nbnd2,nbnd1) -CONJG(hmat(nbnd1,nbnd2))
! !             vsicah( nbnd1, nbnd2) = hmat(nbnd1,nbnd2) -CONJG(hmat(nbnd2,nbnd1))
! !           ENDIF
! !           vsicah2sum            =  vsicah2sum + 2.0d0*CONJG(vsicah(nbnd2,nbnd1))*vsicah(nbnd2,nbnd1)
! !           !
! !       enddo
! !       enddo
!       !
!       deallocate( h0c0 )
!
!       !
!       call stop_clock('nk_get_vsicah')
!       !
!       return
!       !
! !---------------------------------------------------------------
! end subroutine nksic_getvsicah_twin
! !---------------------------------------------------------------

!-----------------------------------------------------------------------
   subroutine nksic_getOmat1(isp, Heig, Umat, passof, Omat1, lgam)
!-----------------------------------------------------------------------
!
! ... Obtains the rotation matrix from the force-related matrices Heig and Umat
!     and also from the step size (passof).
!
      use kinds, only: dp
      use constants, only: ci
      use electrons_base, only: nupdwn
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in) :: isp
      real(dp), intent(in) :: Heig(nupdwn(isp))
      complex(dp), intent(in) :: Umat(nupdwn(isp), nupdwn(isp))
      real(dp), intent(in) :: passof
      complex(dp)                 :: Omat1(nupdwn(isp), nupdwn(isp))
      logical :: lgam
      !
      ! local variables
      !
      complex(dp) :: Cmattmp(nupdwn(isp), nupdwn(isp))
      complex(dp) :: exp_iHeig(nupdwn(isp))

      integer     :: nbnd1
      real(dp)    :: dtmp

      call start_clock("nk_getOmat1")

!$$ We set the step size in such a way that the phase change
!$$ of the wavevector with largest eigenvalue upon rotation is fixed
!          passof = passoprod/max(abs(Heig(1)),abs(Heig(nupdwn(isp))))
!$$ Now the above step is done outside.

      do nbnd1 = 1, nupdwn(isp)
         dtmp = passof*Heig(nbnd1)
         exp_iHeig(nbnd1) = DCOS(dtmp) + ci*DSIN(dtmp)
      end do

!$$ Cmattmp = exp(i * passof * Heig) * Umat^dagger   ; Omat = Umat * Cmattmp
      do nbnd1 = 1, nupdwn(isp)
         Cmattmp(nbnd1, :) = exp_iHeig(nbnd1)*CONJG(Umat(:, nbnd1))
      end do

!           Omat1 = MATMUL( CONJG(TRANSPOSE(Umat)), Cmattmp) !modified:giovanni
      IF (lgam) THEN
         Omat1 = DBLE(MATMUL(Umat, Cmattmp)) !modified:giovanni
      ELSE
         Omat1 = MATMUL(Umat, Cmattmp) !modified:giovanni
      END IF

      call stop_clock("nk_getOmat1")

      return
!---------------------------------------------------------------
   end subroutine nksic_getOmat1
!---------------------------------------------------------------

   SUBROUTINE compute_nksic_centers(nnrx, nx, nudx, nbsp, nspin, iupdwn, &
                                    nupdwn, ispin, orb_rhor, wfc_centers, wfc_spreads, j, k)

      USE kinds, ONLY: DP
      USE ions_positions, ONLY: taus
      USE ions_base, ONLY: ions_cofmass, pmass, na, nsp
      USE cell_base, ONLY: s_to_r
      USE mp, ONLY: mp_bcast
      !
      !   INPUT VARIABLES
      !
      INTEGER, INTENT(IN)      :: nx, nnrx, nudx
      INTEGER, INTENT(IN)      :: ispin(nx), j, k, nspin, nbsp, &
                                  nupdwn(nspin), iupdwn(nspin)
      !ispin is 1 or 2 for each band (listed as in c0),
      !nx is nudx, j and k the two bands involved in the
      !spread calculation
      REAL(DP), INTENT(in)  :: orb_rhor(nnrx, 2)
      REAL(DP) :: wfc_centers(4, nudx, nspin)
      REAL(DP) :: wfc_spreads(nudx, nspin, 2)
      !orbital spreads: both wannier(1) and self-hartree(2)
      !self-hartree is stored separately, within nksic subroutines
      !
      !INTERNAL VARIABLES
      !
      INTEGER :: myspin1, myspin2, mybnd1, mybnd2
      REAL(DP):: r0(3)
      REAL(DP), external :: ddot

      !
      myspin1 = ispin(j)
      !
      mybnd1 = j - iupdwn(myspin1) + 1
      !
      ! compute ionic center of mass
      !
      CALL ions_cofmass(taus, pmass, na, nsp, r0)
      ! and use it as reference position
      !
      call compute_dipole(nnrx, 1, orb_rhor(1, 1), r0, wfc_centers(1:4, mybnd1, myspin1), &
                          wfc_spreads(mybnd1, myspin1, 1))
!!! NB: NLN: I modify the wfc_spread to quadrupole form, it does not equal to wannier definition
      wfc_spreads(mybnd1, myspin1, 1) = wfc_spreads(mybnd1, myspin1, 1) - &
                                        ddot(3, wfc_centers(2:4, mybnd1, myspin1), 1, wfc_centers(2:4, mybnd1, myspin1), 1)
      !
      ! now shift wavefunction centers by r0
      !
      wfc_centers(2:4, mybnd1, myspin1) = wfc_centers(2:4, mybnd1, myspin1) + r0(1:3)
      !
      IF (k .le. nbsp) THEN
         !
         myspin2 = ispin(k)
         mybnd2 = k - iupdwn(myspin2) + 1
         !
         call compute_dipole(nnrx, 1, orb_rhor(1, 2), r0, wfc_centers(1:4, mybnd2, myspin2), &
                             wfc_spreads(mybnd2, myspin2, 1))
!!! NB: NLN: I modify the wfc_spread to quadrupole form, it does not equal to wannier definition
         wfc_spreads(mybnd2, myspin2, 1) = wfc_spreads(mybnd2, myspin2, 1) - &
                                           ddot(3, wfc_centers(2:4, mybnd2, myspin2), 1, wfc_centers(2:4, mybnd2, myspin2), 1)
         !
         ! now shift wavefunction centers by r0
         !
         wfc_centers(2:4, mybnd2, myspin2) = wfc_centers(2:4, mybnd2, myspin2) + r0(1:3)
         !
      END IF
      !
      RETURN

   END SUBROUTINE compute_nksic_centers
!
   SUBROUTINE spread_sort(ngw, nspin, nbsp, nudx, nupdwn, iupdwn, tempspreads, wfc_centers, sort_spreads)

      USE kinds, ONLY: DP
      USE input_parameters, only: draw_pot, sortwfc_spread !added:linh draw vsic_realspace potentials
      USE wavefunctions_module, only: c0, cm
      USE mp, only: mp_bcast

      IMPLICIT NONE

      !COMPLEX(DP) :: c0(ngw, nbsp), cm(ngw,nbsp)
      INTEGER :: ngw, nspin, nbsp, nudx, nupdwn(nspin), iupdwn(nspin)
      REAL(DP) :: tempspreads(nudx, nspin, 2)
      REAL(DP) :: wfc_centers(4, nudx, nspin)
      INTEGER :: sort_spreads(nudx, nspin)
      !
      INTEGER :: isp, j, k, refnum, i, ig
      INTEGER, ALLOCATABLE :: aidarray(:, :)
      !REAL(DP), ALLOCATABLE :: tempspreads(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: tempwfc(:, :)

      ! do nothing if one is drawing the potential: to avoid mismatch between potential and orbital
      IF (draw_pot) THEN
         return
      END IF
      !
      !allocate(tempspreads(nudx,nspin,2))
      allocate (aidarray(nudx, 2), tempwfc(ngw, 2))
      !
      !tempspreads(:,:,:) = wfc_spreads(:,:,:)
      !
      !write(*,*) mpime, "spreads", tempspreads(:,2,2)
      !
      do isp = 1, nspin
         !
         !if(ionode) then
         do j = 1, nupdwn(isp) !initialize sort-decodification array
            !
            aidarray(j, 1) = j
            aidarray(j, 2) = 0
            !
         end do
         !
         do j = 1, nupdwn(isp) - 1 !bubble-sort the decodification array
            !
            do k = nupdwn(isp), j + 1, -1
               !
               IF (tempspreads(k, isp, 2) .lt. tempspreads(k - 1, isp, 2)) THEN
                  !
                  call swap_real(tempspreads(k, isp, 2), tempspreads(k - 1, isp, 2))
                  call swap_real(tempspreads(k, isp, 1), tempspreads(k - 1, isp, 1))
                  do i = 1, 4
                     call swap_real(wfc_centers(i, k, isp), wfc_centers(i, k - 1, isp))
                  end do
                  call swap_integer(aidarray(k, 1), aidarray(k - 1, 1))
                  !
               END IF
               !
            end do
            !
         end do
         !write(*,*) mpime, "aidarray", aidarray(:,1)
         j = 1
         k = 1
         refnum = 0
         !write(*,*) mpime, "before", c0(2,:)
         !
         if (sortwfc_spread) then
            !
            do while (k .le. nupdwn(isp))
               !
               write (6, *) j, aidarray(j, 2), aidarray(j, 1), refnum
               IF (aidarray(j, 2) == 0 .and. j /= aidarray(j, 1)) THEN
                  !
                  IF (aidarray(j, 1) /= refnum) THEN
                     !
                     IF (refnum == 0) THEN
                        !
                        do ig = 1, ngw
                           !
                           tempwfc(ig, 1) = c0(ig, iupdwn(isp) + j - 1)
                           tempwfc(ig, 2) = cm(ig, iupdwn(isp) + j - 1)
                           !
                        end do
                        refnum = j
                        !
                     END IF
                     !
                     do ig = 1, ngw
                        !
                        c0(ig, iupdwn(isp) + j - 1) = c0(ig, iupdwn(isp) + aidarray(j, 1) - 1)
                        cm(ig, iupdwn(isp) + j - 1) = cm(ig, iupdwn(isp) + aidarray(j, 1) - 1)
                        !
                     end do
                     !
                     aidarray(j, 2) = 1
                     j = aidarray(j, 1)
                     !
                  ELSE
                     !
                     do ig = 1, ngw
                        !
                        c0(ig, iupdwn(isp) + j - 1) = tempwfc(ig, 1)
                        cm(ig, iupdwn(isp) + j - 1) = tempwfc(ig, 2)
                        !
                     end do
                     !
                     aidarray(j, 2) = 1
                     j = refnum + 1
                     refnum = 0
                     !
                  END IF
                  k = k + 1
                  !
               ELSE
                  !
                  IF (j == aidarray(j, 1)) THEN
                     !
                     k = k + 1
                     !
                  END IF
                  !
                  j = j + 1
                  !
                  if (j .gt. nupdwn(isp)) THEN
                     exit
                  ELSE
                     cycle
                  END IF
                  !
               END IF
               !
            end do
         end if
         !
         sort_spreads(:, isp) = aidarray(:, 1)
         !
      end do

      !
      if (allocated(tempwfc)) deallocate (tempwfc)
      deallocate (aidarray)
      !
      return

   contains

      subroutine swap_integer(a, b)

         use kinds, ONLY: DP

         implicit none

         INTEGER :: a, b
         INTEGER :: c

         c = a
         a = b
         b = c

         return

      end subroutine swap_integer

      subroutine swap_real(a, b)

         use kinds, ONLY: DP

         implicit none

         REAL(DP) :: a, b
         REAL(DP) :: c

         c = a
         a = b
         b = c

         return

      end subroutine swap_real

   END SUBROUTINE spread_sort

   SUBROUTINE compute_complexification_index(ngw, nnrx, nnrsx, nbsp, nbspx, nspin, ispin, iupdwn, nupdwn, c0, bec, &
                                             complexification_index)
      !
      ! Here the overlap between the wavefunction manifold and its conjugate is calculated
      !
      ! As it is now, this routine works only with Norm Conserving Pseudopotentials
      !
      USE kinds, ONLY: DP
      USE twin_types
      USE mp, ONLY: mp_sum
      USE mp_global, ONLY: intra_image_comm
      use cell_base, only: omega
      use cp_interfaces, only: fwfft, invfft
      use fft_base, only: dfftp

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ngw, nnrx, nnrsx, nbsp, nbspx, nspin, &
                             iupdwn(nspin), nupdwn(nspin), ispin(nbspx)
      type(twin_matrix) :: bec
      COMPLEX(DP) :: c0(ngw, nbspx), complexification_index

      INTEGER :: i, j, ir
      COMPLEX(DP), allocatable :: temp_array(:, :), psi1(:), psi2(:)
      REAL(DP) :: sa1

      sa1 = 1.0d0/omega
      !
      allocate (temp_array(nbsp, nbsp))
      !
      temp_array = CMPLX(0.d0, 0.d0)
      !
      if (nnrsx == nnrx) then ! warning this is a bad way to say we are using ultrasoft
         !
         allocate (psi1(nnrx), psi2(nnrx))
         !
         do i = 1, nbsp
            !
            do j = 1, i
               !
               IF (ispin(i) == ispin(j)) THEN
                  !
                  call c2psi(psi1, nnrx, c0(:, i), c0(:, j), ngw, 0)
                  call c2psi(psi2, nnrx, c0(:, j), c0(:, i), ngw, 0)
                  !
                  CALL invfft('Dense', psi1, dfftp)
                  CALL invfft('Dense', psi2, dfftp)
                  !
                  do ir = 1, nnrx
                     !
                     temp_array(i, j) = temp_array(i, j) + psi1(ir)*psi2(ir)
                     !
                  end do
                  !
               END IF
               !
            end do
            !
         end do
         !
      else !if using uspp
         !
         allocate (psi1(nnrx), psi2(nnrx))
         !
         ! for the moment: do nothing
         !
      end if
      !
      call mp_sum(temp_array, intra_image_comm)
      !
      temp_array = temp_array/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
      complexification_index = 0.d0
      !
      do i = 1, nbsp
         !
         do j = 1, i - 1
            !
            IF (ispin(j) == ispin(i)) THEN
               !
               complexification_index = complexification_index + 2.d0*abs(temp_array(i, j))**2
               !
            END IF
            !
         end do
         !
         complexification_index = complexification_index + abs(temp_array(i, i))**2
         !
      end do
      !
      complexification_index = (1.d0 - complexification_index/nbsp)*100.d0 ! the index is in percentage
      !
      deallocate (temp_array)
      deallocate (psi1)
      deallocate (psi2)
      !
      return

   END subroutine compute_complexification_index

!-----------------------------------------------------------------------
   subroutine nksic_potential_non_ortho(nbsp, nx, c, cdual, f_diag, &
                                        bec, becdual, becsum, &
                                        deeq_sic, ispin, iupdwn, nupdwn, &
                                        rhor, rhoc, wtot_realspace, wtot_reciprocal, vsic_realspace, vsic_reciprocal, do_wxd_, pink, nudx, &
                                        wfc_centers, wfc_spreads, &
                                        icompute_spread)
!-----------------------------------------------------------------------
!
! ....calculate orbital dependent potentials,
!     following the Non-Koopmans' (NK) scheme,
!     but also Perdew-Zunger (PZ),
!     Non-Koopmans' integral definition (NKI),
!     Non-Joopmans on Perdew Zunger (PZNK)
!
!     subroutine writte for non-orthogonal functions
!     note that non-linear core correction is not working
!     in this particular subroutine
!
!
      use kinds, only: dp
      use gvecp, only: ngm
      use gvecw, only: ngw
      use grid_dimensions, only: nnrx
      use electrons_base, only: nspin
      use funct, only: dft_is_gradient
      use nksic, only: orb_rhor, wxdsic_realspace, wxdsic_reciprocal, &
                       wrefsic, rhoref, rhobar, &
                       do_nk, do_nki, do_pz, do_nkpz, &
                       do_nkipz, grhobar, fion_sic, &
                       pzalpha => odd_alpha, do_pz_renorm, edens, &
                       tauw, taukin, kfact
      use nksic_corrections, only: nksic_correction_nk, nksic_correction_nki, &
                       nksic_correction_nkpz, nksic_correction_nkipz, &
                       nksic_correction_pz
      use ions_base, only: nat
      use control_flags, only: gamma_only, do_wf_cmplx !added:giovanni
      use uspp_param, only: nhm
      use cp_interfaces, only: nksic_get_orbitalrho !added:giovanni
      use twin_types !added:giovanni
      use input_parameters, only: draw_pot, pot_number  !added:linh draw vsic_realspace potentials
      use io_pot_sic_xml, only: write_pot_sic  !added:linh draw vsic_realspace potentials
      USE io_global, ONLY: stdout
      !
      implicit none
      !
      ! in/out vars
      !
      integer, intent(in)  :: nbsp, nx, nudx
      complex(dp), intent(in)  :: c(ngw, nx), cdual(ngw, nx)
      type(twin_matrix), intent(in)  :: bec, becdual!(nkb,nbsp) !modified:giovanni
      real(dp), intent(in)  :: becsum(nhm*(nhm + 1)/2, nat, nspin)
      integer, intent(in)  :: ispin(nx)
      integer, intent(in)  :: iupdwn(nspin), nupdwn(nspin)
      real(dp), intent(in)  :: f_diag(nx)
      real(dp), intent(in)  :: rhor(nnrx, nspin)
      real(dp), intent(in)  :: rhoc(nnrx)
      real(dp), intent(out) :: vsic_realspace(nnrx, nx), wtot_realspace(nnrx, 2)
      complex(dp), intent(out) :: vsic_reciprocal(ngm, nx), wtot_reciprocal(ngm, 2)
      real(dp), intent(out) :: deeq_sic(nhm, nhm, nat, nx)
      logical, intent(in)  :: do_wxd_
      real(dp), intent(out) :: pink(nx)
      logical                  :: icompute_spread
      real(DP) :: wfc_centers(4, nudx, nspin)
      real(DP) :: wfc_spreads(nudx, nspin, 2)

      !
      ! local variables
      !
      integer  :: i, j, jj, ibnd, isp
      real(dp) :: focc, pinkpz, shart
      real(dp), allocatable :: vsicpz_realspace(:)
      complex(dp), allocatable :: vsicpz_reciprocal(:)
      complex(dp), allocatable :: rhobarg(:, :)
      logical :: lgam
      !
      ! main body
      !
      CALL start_clock('nksic_drv')
      lgam = gamma_only .and. .not. do_wf_cmplx
      !
      ! compute potentials
      !
      if (dft_is_gradient()) then
         allocate (rhobarg(ngm, 2))
         !write(6,*) "allocated rhobarg"
      else
         allocate (rhobarg(1, 1))
      end if

      if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
         wtot_realspace = 0.0_dp
         wtot_reciprocal = 0.0_dp
      end if
      !
      if (do_nkpz .or. do_nkipz) then
         allocate (vsicpz_realspace(nnrx))
         allocate (vsicpz_reciprocal(ngm))
         vsicpz_realspace = 0.0_dp
         vsicpz_reciprocal = 0.0_dp
      end if
      !
      pink = 0.0_dp
      vsic_realspace = 0.0_dp
      vsic_reciprocal = 0.0_dp
      !
      ! if using pz_renorm factors, compute here tauw and upsilonw
      !
      if (do_pz_renorm) THEN
         !
         edens = 0.d0
         !
         do isp = 1, nspin
            !
            call nksic_get_taukin_pz(1.d0, nspin, isp, rhor(1, isp), tauw, 1)
            !
         end do
         !
      END IF
      !
      ! loop over bands (2 ffts at the same time)
      !
!       call compute_overlap(c, ngw, nbsp, overlap_)
      !
      do j = 1, nbsp, 2
         !
         ! compute orbital densities
         ! n odd => c(:,n+1) is already set to zero
         !
         call nksic_get_orbitalrho(ngw, nnrx, bec, becdual, ispin, nbsp, &
                                   c(:, j), c(:, j + 1), cdual(:, j), cdual(:, j + 1), orb_rhor, &
                                   j, j + 1, lgam) !note:giovanni change here for non-orthogonality flavour
!begin_added:giovanni
!            orb_rhor(:,1) = orb_rhor(:,1)/overlap_(j+1-iupdwn(ispin(j)),j+1-iupdwn(ispin(j)),ispin(j))
!            orb_rhor(:,2) = orb_rhor(:,2)/overlap_(j+2-iupdwn(ispin(j+1)),j+2-iupdwn(ispin(j+1)),ispin(j+1))
         !compute centers and spreads of nksic or pz minimizing orbitals
         IF (icompute_spread) THEN
            !
            call compute_nksic_centers(nnrx, nx, nudx, nbsp, nspin, iupdwn, &
                                       nupdwn, ispin, orb_rhor, wfc_centers, wfc_spreads, j, j + 1)
            !
         END IF
         !
!end_added:giovanni
         !
         shart = 0.d0
         !
         ! compute orbital potentials
         !
         inner_loop: do jj = 1, 2
            !
            i = j + jj - 1
            !
            ! this condition is important when n is odd
            !
            if (i > nbsp) exit inner_loop
            !
            ibnd = i
            if (nspin == 2) then
               if (i >= iupdwn(2)) ibnd = i - iupdwn(2) + 1
            end if
            !
            ! note: iupdwn(2) is set to zero if nspin = 1
            !
            focc = f_diag(i)*DBLE(nspin)/2.0d0
            !
            ! compute parameters needed for PZ-renormalization
            !
            IF (do_pz_renorm) THEN
               !
               call nksic_get_taukin_pz(focc, nspin, ispin(i), orb_rhor(:, jj), &
                                        taukin, ibnd)
               !
            END IF
            !
            !
            ! define rhoref and rhobar
            !
            call nksic_get_rhoref(i, nnrx, ispin(i), nspin, &
                                  focc, rhor, orb_rhor(:, jj), &
                                  rhoref, rhobar, rhobarg, grhobar)

            !
            ! compute nk pieces to build the potentials and the energy
            !
            if (do_nk .or. do_nkpz) then
               !
               call nksic_correction_nk(focc, ispin(i), orb_rhor(:, jj), &
                                        rhor, rhoref, rhobar, rhobarg, grhobar, &
                                        vsic_realspace(:, i), wxdsic_realspace, wrefsic, do_wxd_, &
                                        pink(i), ibnd, shart)
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
               ! here information is accumulated over states
               ! (wtot_realspace is added in the next loop)
               !
               wtot_realspace(1:nnrx, 1:2) = wtot_realspace(1:nnrx, 1:2) + wxdsic_realspace(1:nnrx, 1:2)
               !
               ! ths sic potential is partly updated here to save some memory
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + wrefsic(1:nnrx) &
                                 - wxdsic_realspace(1:nnrx, ispin(i))
               !
            end if

            !
            ! compute nkpz pieces to build the potential and the energy
            !
            if (do_nkpz) then
               !
               call nksic_correction_nkpz(focc, orb_rhor(:, jj), vsicpz_realspace, &
                                          wrefsic, pinkpz, ibnd, ispin(i))
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + vsicpz_realspace(1:nnrx) &
                                 + wrefsic(1:nnrx)
               !
               pink(i) = pink(i) + pinkpz
               !
            end if

            !
            ! compute pz potentials and energy
            !
            if (do_pz) then
               !
               call nksic_correction_pz(focc, ispin(i), orb_rhor(:, jj), &
                                        vsic_realspace(:, i), vsic_reciprocal(:, i), pink(i), pzalpha(i), ibnd, shart)
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
               if (do_pz_renorm) then
                  !
                  edens(:, ispin(i)) = edens(:, ispin(i)) + pink(i)*orb_rhor(:, jj)
                  !
               end if
               !
            end if

            !
            ! compute nki pieces to build the potentials and the energy
            !
            if (do_nki .or. do_nkipz) then
               !
               call nksic_correction_nki(focc, ispin(i), orb_rhor(:, jj), &
                                         rhor, rhoref, rhobar, rhobarg, grhobar, &
                                         vsic_realspace(:, i), vsic_reciprocal(:, i), wxdsic_realspace, &
                                         wxdsic_reciprocal, do_wxd_, pink(i), ibnd, shart)
               !
               ! here information is accumulated over states
               ! (wtot_realspace is added in the next loop)
               !
               wtot_realspace(1:nnrx, 1:2) = wtot_realspace(1:nnrx, 1:2) + wxdsic_realspace(1:nnrx, 1:2)
               wtot_reciprocal(1:ngm, 1:2) = wtot_reciprocal(1:ngm, 1:2) + wxdsic_reciprocal(1:ngm, 1:2)
               !
               ! ths sic potential is partly updated here to save some memory
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) - wxdsic_realspace(1:nnrx, ispin(i))
               vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i) - wxdsic_reciprocal(1:ngm, ispin(i))
               !
            end if

            if (do_nkipz) then
               !
               call nksic_correction_nkipz(focc, ispin(i), orb_rhor(:, jj), vsicpz_realspace, vsicpz_reciprocal, &
                                           pinkpz, ibnd, shart)
               !
               vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + vsicpz_realspace(1:nnrx)
               vsic_reciprocal(:, i) = vsic_reciprocal(:, i) + vsicpz_reciprocal(:)
               !
               pink(i) = pink(i) + pinkpz
               !
               wfc_spreads(ibnd, ispin(i), 2) = shart
               !
            end if
            !
            ! take care of spin symmetry
            !
            pink(i) = f_diag(i)*pink(i)
            !
            if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
               !
               if (nspin == 1) then
                  !
                  wtot_realspace(1:nnrx, 1) = wtot_realspace(1:nnrx, 1) + wxdsic_realspace(1:nnrx, 2)
                  wtot_realspace(1:nnrx, 2) = wtot_realspace(1:nnrx, 2) + wxdsic_realspace(1:nnrx, 1)
                  !
                  wtot_reciprocal(1:ngm, 1) = wtot_reciprocal(1:ngm, 1) + wxdsic_reciprocal(1:ngm, 2)
                  wtot_reciprocal(1:ngm, 2) = wtot_reciprocal(1:ngm, 2) + wxdsic_reciprocal(1:ngm, 1)
               end if
               !
            end if
            !
         end do inner_loop
         !
      end do

      !
      ! Switch off the icompute_spread flag if present
      !
      IF (icompute_spread) THEN
         !
         icompute_spread = .false.
         !
      END IF
      !
      ! now wtot_realspace is completely built and can be added to vsic
      !
      if (do_nk .or. do_nkpz .or. do_nki .or. do_nkipz) then
         !
         do i = 1, nbsp
            !
            vsic_realspace(1:nnrx, i) = vsic_realspace(1:nnrx, i) + wtot_realspace(1:nnrx, ispin(i))
            vsic_reciprocal(1:ngm, i) = vsic_reciprocal(1:ngm, i) + wtot_reciprocal(1:ngm, ispin(i))
            !
         end do
         !
      end if
      !
      ! if pz is renormalized, here we compute the potential, and multiply here pink by renormalization factor
      !
      if (do_pz_renorm) then
         !
         do j = 1, nbsp, 2
            !
            call nksic_get_orbitalrho(ngw, nnrx, bec, ispin, nbsp, &
                                      c(:, j), c(:, j + 1), orb_rhor, j, j + 1, lgam)
            !
            inner_loop_renorm: do jj = 1, 2
               !
               i = j + jj - 1
               if (i > nbsp) exit inner_loop_renorm
               !
               ibnd = i
               focc = f_diag(i)*DBLE(nspin)/2.0d0
               !
               if (nspin == 2) then
                  if (i >= iupdwn(2)) ibnd = i - iupdwn(2) + 1
               end if
               !
               call nksic_get_pz_factor(nspin, ispin(i), orb_rhor(:, jj), rhor, &
                                        taukin, tauw, pzalpha(i), ibnd, kfact)
               !
               ! update vsic_realspace with factor here: it works for pz, will it work for
               ! nk-type functionals?
               !
!                 vsic_realspace(:,i) = vsic_realspace(:,i)*pzalpha(i)
!                 pink(i) = pink(i)*pzalpha(i)
               !
               !
!                 call nksic_get_pzfactor_potential(focc, nspin, ispin(i), rhor, orb_rhor(:,jj), &
!                                       pink(i), taukin, tauw, edens, upsilonkin, upsilonw, vsic_realspace(:,i), pzalpha(i), ibnd, kfact)
               !
            end do inner_loop_renorm
            !
         end do
         !
      end if
      !
      !
      if (draw_pot) then !added:linh draw vsic_realspace potentials
         !
         write (stdout, *) "I am writing out vsic", nbsp
         do i = 1, nbsp
            !
            if (i == pot_number) call write_pot_sic(vsic_realspace(:, i))
            !
         end do
         !
      end if !added:linh draw vsic_realspace potentials
      !
      if (allocated(vsicpz_realspace)) deallocate (vsicpz_realspace)
      !
      if (allocated(vsicpz_reciprocal)) deallocate (vsicpz_reciprocal)
      !
      call ebl_check(vsic_realspace(:, 10), vsic_reciprocal(:, 10))

      !
      ! USPP:
      ! compute corrections to the D coefficients of the pseudopots
      ! due to vsic_realspace(r, i) in the case of orbital dependent functionals.
      ! The corresponding contributions to the forces are computed.
      !
      ! IMPORTANT: the following call makes use of newd.
      !            It must be done before we call newd for the
      !            total potentials, because deeq is overwritten at every call
      !
      fion_sic(:, :) = 0.0d0
      !
      IF (nhm > 0) then
         !
         deeq_sic(:, :, :, :) = 0.0d0
         !
         DO i = 1, nbsp
            !
            CALL nksic_newd(i, nnrx, ispin(i), nspin, vsic_realspace(:, i), vsic_reciprocal(:, i), nat, nhm, &
                            becsum, fion_sic, deeq_sic(:, :, :, i))
            !
         END DO
         !
      END IF
      !
      deallocate (rhobarg)
      !
      CALL stop_clock('nksic_drv')
      return
      !
!-----------------------------------------------------------------------
   end subroutine nksic_potential_non_ortho
!-----------------------------------------------------------------------

   subroutine ebl_check(vsic_realspace, vsic_reciprocal)

      use kinds, only: dp
      use cp_interfaces, only: invfft
      use fft_base, only: dfftp
      use gvecp, only: ngm
      use grid_dimensions, only: nnrx
      ! use ifcore, only: tracebackqq
      use mp, only: mp_sum
      use mp_global, only: intra_pool_comm

      implicit none

      real(dp), intent(in) :: vsic_realspace(nnrx)
      complex(dp), intent(in) :: vsic_reciprocal(ngm)

      real(dp), allocatable :: vsic_realspace_local(:)
      real(dp) :: sumdiff
      complex(dp), allocatable :: psi(:)
      integer :: return_code

      allocate (vsic_realspace_local(size(vsic_realspace)))
      allocate (psi(size(vsic_realspace)))

      call rho2psi('Dense', psi, dfftp%nnr, vsic_reciprocal, ngm)
      call invfft('Dense', psi, dfftp)

      vsic_realspace_local = dble(psi)

      sumdiff = sum(vsic_realspace_local - vsic_realspace)
      call mp_sum(sumdiff, intra_pool_comm)
      write (*, *) 'ebl_check:', sumdiff

      return_code = -1
      if (abs(sumdiff) > 1d-8) return_code = 0

      ! call tracebackqq(user_exit_code=return_code)

      deallocate (vsic_realspace_local, psi)
   end subroutine ebl_check
