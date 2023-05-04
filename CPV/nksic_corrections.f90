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

module nksic_corrections

   implicit none

contains
!---------------------------------------------------------------
   subroutine nksic_correction_pz(f, ispin, orb_rhor, &
                                  vsic, pink, pzalpha, ibnd, shart)
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds, only: dp
      use constants, only: e2, fpi, hartree_si, electronvolt_si
      use cell_base, only: tpiba2, omega
      use nksic, only: etxc => etxc_sic, vxc => vxc_sic, nknmax, &
                       nkscalfact, do_pz_renorm
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use gvecp, only: ngm
      use recvecs_indexes, only: np, nm
      use reciprocal_vectors, only: gstart, g
      use eecp_mod, only: do_comp
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dft_is_gradient
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use control_flags, only: gamma_only, do_wf_cmplx
      use control_flags, only: hartree_only_sic
      !
      implicit none
      integer, intent(in)  :: ispin, ibnd
      real(dp), intent(in)  :: f, orb_rhor(nnrx), pzalpha
      complex(dp), intent(out) :: vsic(ngm)
      real(dp), intent(out) :: pink, shart
      !
      !character(19) :: subname='nksic_correction_pz'
      integer       :: ig
      real(dp)      :: ehele, fact
      !
      complex(dp), allocatable :: rhogaux(:, :)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: vtmp(:)
      real(dp), allocatable :: vsic_realspace(:)
      !
      real(dp), allocatable :: grhoraux(:, :, :)
      real(dp), allocatable :: haux(:, :, :)
      logical :: lgam
      real(dp) :: dexc_dummy(3, 3)
      !
      !==================
      ! main body
      !==================
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      vsic = 0.0_dp
      pink = 0.0_dp
      !
      if (ibnd > nknmax .and. nknmax .ge. 0) return
      if (f < 1.0d-6) return
      !
      CALL start_clock('nk_corr')
      CALL start_clock('nk_corr_h')
      !
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      !allocate(rhoelef(nnrx,2))
      allocate (rhogaux(ngm, 2))
      allocate (vtmp(ngm))
      allocate (vcorr(ngm))
      allocate (vhaux(nnrx))
      allocate(vsic_realspace(nnrx), source=0.0_dp)
      !
      !rhoelef=0.0d0
      !rhoelef(:,ispin) = f * orb_rhor(:)
      !
      ! Compute self-hartree contributions
      !
      rhogaux = 0.0_dp
      !
      ! rhoelef contains occupations
      !
      !vhaux(:) = rhoelef(:,ispin)
      vhaux(:) = f*orb_rhor(:)
      !
      call fwfft('Dense', vhaux, dfftp)
      !
      do ig = 1, ngm
         rhogaux(ig, ispin) = vhaux(np(ig))
      end do

      !
      ! compute hartree-like potential
      !
      if (gstart == 2) vtmp(1) = (0.d0, 0.d0)
      do ig = gstart, ngm
         vtmp(ig) = rhogaux(ig, ispin)*fpi/(tpiba2*g(ig))
      end do
      !
      ! compute periodic corrections
      !
      if (do_comp) then
         !
         call calc_compensation_potential(vcorr, rhogaux(:, ispin), .true.)
         vtmp(:) = vtmp(:) + vcorr(:)
         !
      end if
      !
      vhaux = 0.0_dp
!       if(lgam) then  !!!### uncomment for k points
      do ig = 1, ngm
         !
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG(vtmp(ig))
         !
      end do
!       else !!!### uncomment for k points
!           do ig=1,ngm !!!### uncomment for k points
      !
!               vhaux(np(ig)) = vtmp(ig) !!!### uncomment for k points
!               vhaux(nm(ig)) = conjg(vtmp(ig))
      !
!           enddo !!!### uncomment for k points
!       endif !!!### uncomment for k points
      call invfft('Dense', vhaux, dfftp)
      !
      ! init vsic
      !
      vsic_realspace(1:nnrx) = -DBLE(vhaux(1:nnrx))
      ehele = 0.5_dp*sum(DBLE(vhaux(1:nnrx)) &
                         *orb_rhor(1:nnrx))
      !
      ! set ehele as measure of spread
      !
      !IF(icompute_spread) THEN
      shart = abs(ehele)*fact*hartree_si/electronvolt_si
      call mp_sum(shart, intra_image_comm)
      !ENDIF
      !
      ehele = ehele*f !this is to make ehele quadratic in f (check this)
      !
      ! partial cleanup
      !
      deallocate (vtmp)
      deallocate (vcorr)
      deallocate (vhaux)
      !
      CALL stop_clock('nk_corr_h')
      !
      ! Compute xc-contributions
      !
      if (.not. hartree_only_sic) then
         !
         if (dft_is_gradient()) then
            !
            allocate (grhoraux(nnrx, 3, 2))
            !
            allocate (haux(nnrx, 2, 2))
            !
            ! note: rhogaux contains the occupation
            !
            grhoraux = 0.0_dp
            call fillgrad(1, rhogaux(:, ispin:ispin), grhoraux(:, :, ispin:ispin), lgam)
            !
            !
         else
            allocate (grhoraux(1, 1, 1))
            allocate (haux(1, 1, 1))
            !
            grhoraux = 0.0_dp
         end if
         !
         !
         vxc = 0.0_dp
         haux = 0.0_dp
         etxc = 0.0_dp
         !
         vxc(:, ispin) = f*orb_rhor(:)
         ! call exch_corr_wrapper(nnrx,2,grhoraux,rhoelef,etxc,vxc,haux)
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxc, etxc) !proposed:giovanni fixing PBE, warning, check array dimensions
         !
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
            call gradh(2, grhoraux, rhogaux, vxc, dexc_dummy, lgam)
            !  grhoraux(nnr,3,nspin)?yes; rhogaux(ng,nspin)? rhoref(nnr, nspin)
            !
         end if
!$$
         vsic_realspace(1:nnrx) = vsic_realspace(1:nnrx) - vxc(1:nnrx, ispin)

      else
         !
         etxc = 0.
         !
      end if
!      vsic_realspace(1:nnrx) = -vxc(1:nnrx,ispin)
!$$
      !
      ! energy correction terms
      !
!$$
      pink = fact*(-etxc - ehele)
!$$
!      pink = fact * ( -ehele )
!      pink = fact * ( -etxc )
!$$

!$$ This is for screened pz functional; apparently, I should have used a different variable name.
      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      IF (.not. do_pz_renorm) THEN
         !
         pink = pink*nkscalfact
         vsic_realspace = vsic_realspace*nkscalfact
         !
      ELSE
         !I do not renormalize here, I will do it outside the subroutine
         !
      END IF
      !
      call mp_sum(pink, intra_image_comm)
      !
      if (.not. hartree_only_sic) then
         deallocate (grhoraux)
         deallocate (haux)
      end if
      deallocate (rhogaux)
      !
      CALL stop_clock('nk_corr')
      !
      ! Store the vsic potential in reciprocal space
      vtmp = vsic_realspace(:)
      call fwfft('Dense', vtmp, dfftp)
      call psi2rho('Dense', vtmp, dfftp%nnr, vsic, ngm)
      deallocate(vsic_realspace, vtmp)
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_pz
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_nkipz(f, ispin, orb_rhor, &
                                     vsic, pink, ibnd, shart, is_empty)
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds, only: dp
      use constants, only: e2, fpi, hartree_si, electronvolt_si
      use cell_base, only: tpiba2, omega
      use nksic, only: nknmax, nkscalfact
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use gvecp, only: ngm
      use recvecs_indexes, only: np, nm
      use reciprocal_vectors, only: gstart, g
      use eecp_mod, only: do_comp
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dft_is_gradient
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use control_flags, only: gamma_only, do_wf_cmplx, hartree_only_sic
      !
      implicit none
      integer, intent(in)  :: ispin, ibnd
      real(dp), intent(in)  :: f, orb_rhor(nnrx)
      complex(dp), intent(out) :: vsic(ngm)
      real(dp), intent(out) :: pink, shart
      logical, optional, intent(in) :: is_empty
      !
      !character(19) :: subname='nksic_correction_pz'
      integer       :: ig
      real(dp)      :: ehele, fact, w2cst, etmp, etxc_
      !
      real(dp), allocatable :: vxc_(:, :)
      complex(dp), allocatable :: rhogaux(:, :)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: vtmp(:)
      complex(dp), allocatable :: psi(:)
      real(dp), allocatable    :: vsic_realspace(:)
      !
      real(dp), allocatable :: grhoraux(:, :, :)
      logical :: lgam
      real(dp) :: icoeff
      real(dp) :: dexc_dummy(3, 3)
      logical :: is_empty_
      !
      !==================
      ! main body
      !==================
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      if (lgam) then
         icoeff = 2.d0
      else
         icoeff = 1.d0
      end if
      !
      IF (present(is_empty)) THEN
         !
         is_empty_ = is_empty
         !
      ELSE
         !
         is_empty_ = .false.
         !
      END IF
      !
      vsic = 0.0_dp
      pink = 0.0_dp
      !
      if (ibnd > nknmax .and. nknmax .ge. 0) return
      !
      CALL start_clock('nk_corr')
      CALL start_clock('nk_corr_h')
      !
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      allocate (rhogaux(ngm, 2))
      allocate (vtmp(ngm))
      allocate (vcorr(ngm))
      allocate (vxc_(nnrx, 2))
      allocate (vhaux(nnrx))
      !
      ! Compute self-hartree contributions
      !
      rhogaux = 0.0_dp
      !
      ! vhaux does not contain occupations
      !
      vhaux(:) = orb_rhor(:)
      !
      call fwfft('Dense', vhaux, dfftp)
      !
      do ig = 1, ngm
         rhogaux(ig, ispin) = vhaux(np(ig))
      end do
      !
      ! compute hartree-like potential
      !
      if (gstart == 2) vtmp(1) = (0.d0, 0.d0)
      do ig = gstart, ngm
         vtmp(ig) = rhogaux(ig, ispin)*fpi/(tpiba2*g(ig))
      end do
      !
      ! compute periodic corrections
      !
      if (do_comp) then
         !
         call calc_compensation_potential(vcorr, rhogaux(:, ispin), .true.)
         vtmp(:) = vtmp(:) + vcorr(:)
         !
      end if
      !
      vhaux = 0.0_dp
      do ig = 1, ngm
         !
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG(vtmp(ig))
         !
      end do
      !
      call invfft('Dense', vhaux, dfftp)
      !
      ! init vsic
      !
      allocate(vsic_realspace(nnrx), source=0.0_dp)
      vsic_realspace(1:nnrx) = -DBLE(vhaux(1:nnrx))   ! -v_hartree[n_i](r)
      !
      ehele = icoeff*DBLE(DOT_PRODUCT(vtmp(1:ngm), rhogaux(1:ngm, ispin)))
      if (gstart == 2) ehele = ehele + (1.d0 - icoeff)*DBLE(CONJG(vtmp(1))*rhogaux(1, ispin))
      !
      w2cst = 0.0_dp
      !
      w2cst = 0.5_dp*ehele*omega  ! -E_H[n_i] + \int( v_H[n_i](r) n_i(r) )dr --> -E_H[n_i] + 2E_H[n_i] = E_H[n_i]
      !
      call mp_sum(w2cst, intra_image_comm)
      !
      vsic_realspace = vsic_realspace + w2cst    ! Hartree part of first and third terms in eq. (A15) Borghi PRB
      !
      ehele = 0.5d0*ehele*omega/fact
      !
      shart = abs(ehele)*fact*hartree_si/electronvolt_si
      !
      call mp_sum(shart, intra_image_comm)
      !
      ! NsC >>>
      etmp = 0.D0
      etmp = sum(vsic_realspace(1:nnrx)*orb_rhor(1:nnrx))
      etmp = etmp*fact*hartree_si/electronvolt_si
      call mp_sum(etmp, intra_image_comm)
      ! NsC <<<
      ! partial cleanup
      !
      deallocate (vtmp)
      deallocate (vcorr)
      deallocate (vhaux)
      !
      CALL stop_clock('nk_corr_h')
      !
      IF (.NOT. hartree_only_sic) THEN
         ! Compute xc-contributions
         !
         if (dft_is_gradient()) then
            allocate (grhoraux(nnrx, 3, 2))
            !
            ! note: rhogaux does not contain the occupation
            !
            grhoraux = 0.0_dp
            call fillgrad(1, rhogaux(:, ispin:ispin), grhoraux(:, :, ispin:ispin), lgam)
         else
            allocate (grhoraux(1, 1, 1))
            !
            grhoraux = 0.0_dp
         end if
         !
         !
         vxc_ = 0.0_dp
         etxc_ = 0.0_dp
         !
         vxc_(:, ispin) = orb_rhor(:)
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxc_, etxc_)
         !proposed:giovanni fixing PBE, warning, check array dimensions
         !
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc
            !  should be dexc(3,3), is lgam a variable here?
            call gradh(2, grhoraux, rhogaux, vxc_, dexc_dummy, lgam)
            !
         end if
         !
      ELSE
         !
         vxc_ = 0.D0
         etxc_ = 0.D0
         !
      END IF
      !
      IF (.not. is_empty_) THEN
         !
         etmp = sum(vxc_(1:nnrx, ispin)*orb_rhor(1:nnrx))
         !
         w2cst = -etxc_ + etmp
         w2cst = w2cst*fact
         !
         call mp_sum(w2cst, intra_image_comm)
         !
         pink = -f*(etxc_ + ehele)
         !
      ELSE
         !
         etmp = sum(vxc_(1:nnrx, ispin)*orb_rhor(1:nnrx))
         !
         w2cst = -etxc_ + etmp
         w2cst = w2cst*fact
         !
         call mp_sum(w2cst, intra_image_comm)
         !
         pink = -(etxc_ + ehele)
         !
      END IF
      !
      pink = pink*fact
      !
      call mp_sum(pink, intra_image_comm)
      !
      vsic_realspace(1:nnrx) = vsic_realspace(1:nnrx) - vxc_(1:nnrx, ispin) + w2cst
      !
      ! NsC >>>
      etmp = 0.D0
      etmp = sum(vsic_realspace(1:nnrx)*orb_rhor(1:nnrx))
      etmp = etmp*fact*hartree_si/electronvolt_si
      call mp_sum(etmp, intra_image_comm)
      ! NsC <<<
      !
      pink = pink*nkscalfact
      vsic_realspace = vsic_realspace*nkscalfact
      !
      IF (.not. hartree_only_sic) deallocate (grhoraux)
      deallocate (rhogaux)
      deallocate (vxc_)
      !
      CALL stop_clock('nk_corr')
      !
      ! Store the vsic potential in reciprocal space
      allocate(psi(nnrx))
      psi = vsic_realspace(:)
      call fwfft('Dense', psi, dfftp)
      call psi2rho('Dense', psi, dfftp%nnr, vsic, ngm)
      deallocate(vsic_realspace, psi)
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_nkipz
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_nki(f, ispin, orb_rhor, rhor, rhoref, rhobar, rhobarg, grhobar, &
                                   vsic, wxdsic, do_wxd_, pink, ibnd, shart, is_empty)
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans (integrated, NKI)
!     potential from the orbital density
!
!     note that fref=1.0 when performing NKI (i.e. it has a diff
!     meaning)
!     then  rho_ref = rho - rho_i + n_i
!           rho_bar = rho - rho_i
!
      use kinds, only: dp
      use constants, only: e2, fpi, hartree_si, electronvolt_si
      use cell_base, only: tpiba2, omega
      use nksic, only: fref, rhobarfact, nknmax, &
                       nkscalfact, &
                       etxc => etxc_sic, vxc => vxc_sic
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use gvecp, only: ngm
      use recvecs_indexes, only: np, nm
      use reciprocal_vectors, only: gstart, g
      use eecp_mod, only: do_comp
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dmxc_spin, dft_is_gradient
      use mp, only: mp_sum
      use mp_global, only: intra_image_comm
      use electrons_base, only: nspin
      use control_flags, only: gamma_only, do_wf_cmplx
      !
      implicit none
      integer, intent(in)  :: ispin, ibnd
      real(dp), intent(in)  :: f, orb_rhor(nnrx)
      real(dp), intent(in)  :: rhor(nnrx, nspin)
      real(dp), intent(in)  :: rhoref(nnrx, 2)
      real(dp), intent(in)  :: rhobar(nnrx, 2)
      complex(dp), intent(in)  :: rhobarg(ngm, 2)
      real(dp), intent(in)  :: grhobar(nnrx, 3, 2)
      complex(dp), intent(out) :: vsic(ngm)
      complex(dp), intent(out) :: wxdsic(ngm, 2)
      logical, intent(in)  :: do_wxd_
      real(dp), intent(out) :: pink, shart
      logical, optional, intent(in) :: is_empty
      !
      integer       :: ig, i
      real(dp)      :: fact, ehele, etmp
      real(dp)      :: etxcref, etxc0, w2cst
      !
      real(dp), allocatable :: rhoele(:, :)
      real(dp), allocatable :: rhoraux(:, :)
      real(dp), allocatable :: vxc0(:, :)
      real(dp), allocatable :: vxcref(:, :)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:, :)
      complex(dp), allocatable :: vtmp(:)
      !
      real(dp), allocatable :: grhoraux(:, :, :)
      real(dp), allocatable :: orb_grhor(:, :, :)
      complex(dp), allocatable :: orb_rhog(:, :)
      real(dp), allocatable :: haux(:, :, :)
      logical :: lgam, is_empty_
      real(dp) :: icoeff
      real(dp) :: dexc_dummy(3, 3)
      real(dp), allocatable :: vsic_realspace(:)
      real(dp), allocatable :: wxdsic_realspace(:, :)
      complex(dp), allocatable :: psi(:)
      !
      !==================
      ! main body
      !==================
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      !
      if (lgam) then
         icoeff = 2.d0
      else
         icoeff = 1.d0
      end if
      !
      IF (present(is_empty)) THEN
         !
         is_empty_ = is_empty
         !
      ELSE
         !
         is_empty_ = .false.
         !
      END IF
      !
      if (ibnd > nknmax .and. nknmax .ge. 0) return
      !
      CALL start_clock('nk_corr')
      CALL start_clock('nk_corr_h')
      !
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      allocate (rhoele(nnrx, 2))
      allocate (rhogaux(ngm, 2))
      allocate (vtmp(ngm))
      allocate (orb_rhog(ngm, 1))
      allocate (vcorr(ngm))
      allocate (vhaux(nnrx))
      !
      rhoele = 0.0d0
      rhoele(:, ispin) = orb_rhor(:)
      !
      allocate(vsic_realspace(nnrx), source=0.0_dp)
      vsic = 0.0_dp
      allocate(wxdsic_realspace(nnrx, 2), source=0.0_dp)
      wxdsic = 0.0_dp
      pink = 0.0_dp
      !
      ! Compute self-hartree contributions
      !
      orb_rhog = 0.0_dp
      !
      ! rhoele has no occupation
      !
      vhaux(:) = rhoele(:, ispin)
      !
      call fwfft('Dense', vhaux, dfftp)
      !
      do ig = 1, ngm
         orb_rhog(ig, 1) = vhaux(np(ig))
      end do
      !
      ! compute hartree-like potential
      !
      if (gstart == 2) vtmp(1) = (0.d0, 0.d0)
      do ig = gstart, ngm
         vtmp(ig) = orb_rhog(ig, 1)*fpi/(tpiba2*g(ig))
      end do
      !
      ! compute periodic corrections
      !
      if (do_comp) then
         !
         call calc_compensation_potential(vcorr, orb_rhog(:, 1), .true.)
         vtmp(:) = vtmp(:) + vcorr(:)
         !
      end if
      !
      vhaux = 0.0_dp
      do ig = 1, ngm
         !
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG(vtmp(ig))
         !
      end do
      !
      call invfft('Dense', vhaux, dfftp)
      !
      ! init here vsic_realspace to save some memory
      !
      ! this is just the self-hartree potential
      !
      vsic_realspace(1:nnrx) = (1.0_dp - f)*DBLE(vhaux(1:nnrx))
      !
      ! self-hartree contrib to pink
      ! and w2cst for vsic
      !
      ehele = icoeff*DBLE(DOT_PRODUCT(vtmp(1:ngm), orb_rhog(1:ngm, 1)))
      if (gstart == 2) ehele = ehele + (1.d0 - icoeff)*DBLE(CONJG(vtmp(1))*orb_rhog(1, 1))
      !
      shart = abs(ehele)*omega*0.5d0*hartree_si/electronvolt_si
      !
      call mp_sum(shart, intra_image_comm)
      ! -self-hartree energy to be added to the vsic_realspace potential
      !
      ! the scalar Hatree term of both empty and occupied states is
      ! in the same form: -E_H[n_i]
      !
      w2cst = 0.0_dp
      !
      w2cst = -0.5_dp*ehele*omega
      !
      call mp_sum(w2cst, intra_image_comm)
      !
      vsic_realspace = vsic_realspace + w2cst
      !
      ! the f * (1-f) term is added here
      !
      IF (.not. is_empty_) THEN
         !
         ehele = 0.5_dp*f*(1.0_dp - f)*ehele*omega/fact
         !
      ELSE !this is for the fake functional for empty states
         !
         ehele = 0.5_dp*ehele*omega/fact
         !
      END IF
      !
      ! NsC >>>
      etmp = 0.D0
      etmp = sum(vsic_realspace(1:nnrx)*orb_rhor(1:nnrx))
      etmp = etmp*fact*hartree_si/electronvolt_si
      call mp_sum(etmp, intra_image_comm)
      ! NsC <<<
      !
      deallocate (vtmp)
      deallocate (vcorr)
      deallocate (vhaux)
      !
      CALL stop_clock('nk_corr_h')
      !
      CALL start_clock('nk_corr_vxc')
      !
      !
      !   add self-xc contributions
      !
      if (dft_is_gradient()) then
         !
         allocate (grhoraux(nnrx, 3, 2))
         allocate (orb_grhor(nnrx, 3, 1))
         allocate (haux(nnrx, 2, 2))
         !
         ! compute the gradient of n_i(r)
         call fillgrad(1, orb_rhog, orb_grhor(:, :, 1:1), lgam)
         !
      else
         !
         allocate (grhoraux(1, 1, 1))
         allocate (haux(1, 1, 1))
         grhoraux = 0.0_dp
         !
      end if
      !
      allocate (vxc0(nnrx, 2))
      allocate (vxcref(nnrx, 2))
      !
      ! this term is computed for ibnd, ispin == 1 and stored
      ! or if rhobarfact < 1
      !
      if ((ibnd == 1 .and. ispin == 1) .OR. rhobarfact < 1.0_dp) then
         !
         etxc = 0.0_dp
         vxc = 0.0_dp
         !
         ! some meory can be same in the nspin-2 case,
         ! considering that rhobar + f*rhoele is identical to rho
         ! when rhobarfact == 1
         !
         ! call exch_corr_wrapper(nnrx,2,grhoraux,rhor,etxc,vxc,haux)
         !
         if (dft_is_gradient()) then
            !
            grhoraux(:, :, 1:2) = grhobar(:, :, 1:2)
            grhoraux(:, :, ispin) = grhobar(:, :, ispin) &
                                    + f*orb_grhor(:, :, 1)
            !
            rhogaux(:, 1:2) = rhobarg(:, 1:2)
            rhogaux(:, ispin) = rhobarg(:, ispin) + f*orb_rhog(:, 1)
            !
         end if
         !
         allocate (rhoraux(nnrx, 2))
         !
         rhoraux = rhobar + f*rhoele
         vxc = rhoraux
         !
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxc, etxc)
         !proposed:giovanni warning rhoraux is overwritten with vxc,
         !check array dimensions
         !
         !begin_added:giovanni fixing PBE potential
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc
            !  should be dexc(3,3), is lgam a variable here?
            !
            call gradh(2, grhoraux, rhogaux, vxc, dexc_dummy, lgam)
            !
         end if
         !end_added:giovanni fixing PBE potential
         deallocate (rhoraux)
         !
      end if
      !
      etxcref = 0.0_dp
      vxcref = 0.0_dp
      !
      if (f == 1.0_dp) then
         !
         vxcref = vxc
         etxcref = etxc
         !
      else
         !
         if (dft_is_gradient()) then
            !
            grhoraux(:, :, 1:2) = grhobar(:, :, 1:2)
            grhoraux(:, :, ispin) = grhobar(:, :, ispin) &
                                    + fref*orb_grhor(:, :, 1)
            !
            rhogaux(:, 1:2) = rhobarg(:, 1:2)
            rhogaux(:, ispin) = rhobarg(:, ispin) + fref*orb_rhog(:, 1)
            !
         end if
         !
         vxcref = rhoref
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxcref, etxcref)
         !
         !proposed:giovanni warning rhoraux is overwritten with vxc,
         !check array dimensions
         !
         !begin_added:giovanni fixing PBE potential
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc
            !  should be dexc(3,3), is lgam a variable here?
            call gradh(2, grhoraux, rhogaux, vxcref, dexc_dummy, lgam)
            !  grhobar(nnr,3,nspin)? rhogbar(ng,nspin)? rhor(nnr, nspin)
            !
         end if
         !end_added:giovanni fixing PBE potential
         !
      end if
      !
      !rhoraux = rhobar
      !
      etxc0 = 0.0_dp
      vxc0 = 0.0_dp
      !
      vxc0 = rhobar
      CALL exch_corr_cp(nnrx, 2, grhobar, vxc0, etxc0)
      !proposed:giovanni
      !
      !begin_added:giovanni fixing PBE potential
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc
         !  should be dexc(3,3), is lgam a variable here?
         call gradh(2, grhobar, rhobarg, vxc0, dexc_dummy, lgam)
         !  grhobar(nnr,3,nspin)? rhogbar(ng,nspin)? rhor(nnr, nspin)
         !
      end if
      !end_added:giovanni fixing PBE potential
      !
      ! update potential (including other constant terms)
      ! and define pink
      !
      IF (.not. is_empty_) THEN
         !
         etmp = sum(vxcref(1:nnrx, ispin)*rhoele(1:nnrx, ispin))
         w2cst = (etxcref - etxc0) - etmp
         w2cst = w2cst*fact
         !
         call mp_sum(w2cst, intra_image_comm)
         !
         pink = (1.0_dp - f)*etxc0 - etxc + f*etxcref + ehele
         !
      ELSE
         !
         etmp = sum(vxcref(1:nnrx, ispin)*rhoele(1:nnrx, ispin))
         w2cst = (etxcref - etxc0) - etmp
         w2cst = w2cst*fact
         !
         call mp_sum(w2cst, intra_image_comm)
         !
         etmp = sum(vxc(1:nnrx, ispin)*rhoele(1:nnrx, ispin))
         !
         pink = etxcref - etxc0 - etmp + ehele
         !
      END IF
      !
      pink = pink*fact
      !
      call mp_sum(pink, intra_image_comm)
      !
      vsic_realspace(1:nnrx) = vsic_realspace(1:nnrx) &
                               + vxcref(1:nnrx, ispin) - vxc(1:nnrx, ispin) + w2cst
      !
      !   calculate wxd
      !
      if (do_wxd_) then
         !
         wxdsic_realspace(:, 1:2) = (1.0_dp - f)*vxc0(:, 1:2) - vxc(:, 1:2) + f*vxcref(:, 1:2)
         !
         ! Store the wxdsic potential in reciprocal space
         allocate (psi(nnrx))
         do i = 1, 2
            psi = wxdsic_realspace(:, i)
            call fwfft('Dense', psi, dfftp)
            call psi2rho('Dense', psi, dfftp%nnr, wxdsic(:, i), ngm)
         end do
         deallocate (psi)
      else
         wxdsic(:, :) = 0.0d0
      end if
      !
      call stop_clock('nk_corr_vxc')
      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink*nkscalfact
      vsic_realspace = vsic_realspace*nkscalfact
      !
      if (do_wxd_) then
         !
         wxdsic = wxdsic*nkscalfact
         !
      else
         !
         wxdsic = 0.d0
         !
      end if
      !
      ! Store the vsic potential in reciprocal space
      allocate (psi(nnrx))
      psi = vsic_realspace(:)
      call fwfft('Dense', psi, dfftp)
      call psi2rho('Dense', psi, dfftp%nnr, vsic, ngm)
      deallocate(psi, vsic_realspace)
      !
      deallocate (vxc0)
      deallocate (vxcref)
      deallocate (rhoele)
      !
      deallocate (grhoraux)
      deallocate (haux)
      deallocate (rhogaux)
      !
      if (allocated(orb_grhor)) deallocate (orb_grhor)
      !
      CALL stop_clock('nk_corr')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_nki
!---------------------------------------------------------------

end module nksic_corrections
