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
   subroutine nksic_correction_nk(f, ispin, orb_rhor, rhor, &
                                  rhoref, rhobar, rhobarg, grhobar, &
                                  vsic_realspace, wxdsic_realspace, wrefsic, do_wxd_, &
                                  pink, ibnd, shart)
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential from the orbital density
!
      use kinds, only: dp
      use constants, only: e2, fpi, hartree_si, electronvolt_si
      use cell_base, only: tpiba2, omega
      use nksic, only: fref, rhobarfact, nknmax, &
                       vanishing_rho_w, &
                       nkscalfact, do_wref, &
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
      real(dp), intent(out) :: vsic_realspace(nnrx), wrefsic(nnrx)
      real(dp), intent(out) :: wxdsic_realspace(nnrx, 2)
      logical, intent(in)  :: do_wxd_
      real(dp), intent(out) :: pink
      !
      !character(19) :: subname='nksic_correction_nk'
      integer       :: ig, ir
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
      logical :: lgam !!added:giovanni
      real(dp) :: icoeff
      real(dp) :: dexc_dummy(3, 3)
      real(dp) :: shart
      !
      !==================
      ! main body
      !==================
      !
      lgam = gamma_only .and. .not. do_wf_cmplx !added:giovanni
      if (lgam) then
         icoeff = 2.d0
      else
         icoeff = 1.d0
      end if

      if (ibnd > nknmax .and. nknmax .ge. 0) return
      !
      CALL start_clock('nk_corr')
      CALL start_clock('nk_corr_h')

      !
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      allocate (rhoele(nnrx, 2))
      allocate (rhogaux(ngm, 2))
      allocate (orb_rhog(ngm, 1))
      allocate (vtmp(ngm))
      allocate (vcorr(ngm))
      allocate (vhaux(nnrx))
      !
      rhoele = 0.0d0
      rhoele(:, ispin) = orb_rhor(:)
      !
      vsic_realspace = 0.0_dp
      wrefsic = 0.0_dp
      wxdsic_realspace = 0.0_dp
      pink = 0.0_dp

      !
      ! Compute self-hartree contributions
      !
      orb_rhog = 0.0_dp
      !
      ! rhoele has no occupation
      !
      ! f-fref is NOT included here in vhaux
      ! (will be added afterwards)
      !
      vhaux = 0.d0
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

      vhaux = 0.0_dp
!       IF(lgam) THEN !!!### uncomment for k points
      do ig = 1, ngm
         !
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG(vtmp(ig))
         !
      end do
!       ELSE !!!### uncomment for k points
!           do ig=1,ngm !!!### uncomment for k points
      !
!               vhaux(np(ig)) = vtmp(ig) !!!### uncomment for k points
!               vhaux(nm(ig)) = conjg(vtmp(ig))
      !
!           enddo !!!### uncomment for k points
!       ENDIF !!!### uncomment for k points

      call invfft('Dense', vhaux, dfftp)
      !
      ! init here wref sic to save some memory
      !
      ! this is just the self-hartree potential
      ! (to be multiplied by fref later on)
      !
      wrefsic(1:nnrx) = DBLE(vhaux(1:nnrx))
      !
      ! self-hartree contrib to pink
      ! and init vsic
      !
      !ehele=0.5_dp * sum(dble(vhaux(1:nnrx))*rhoele(1:nnrx,ispin))
      !
      ehele = icoeff*DBLE(DOT_PRODUCT(vtmp(1:ngm), orb_rhog(1:ngm, 1)))
      if (gstart == 2) ehele = ehele + (1.d0 - icoeff)*DBLE(CONJG(vtmp(1))*orb_rhog(1, 1))
      !
      shart = 0.5_dp*ehele*omega*hartree_si/electronvolt_si
      call mp_sum(shart, intra_image_comm)

      ! the f * (2.0d0 * fref-f) term is added here
      ehele = 0.5_dp*f*(2.0_dp*fref - f)*ehele*omega/fact
      !shart = 0.5_dp * ehele * omega / fact

      !
      ! fref-f has to be included explicitly in rhoele
      !
      vsic_realspace(1:nnrx) = (fref - f)*DBLE(vhaux(1:nnrx))

      deallocate (vtmp)
      deallocate (vcorr)
      deallocate (vhaux)
      !
      CALL stop_clock('nk_corr_h')

      CALL start_clock('nk_corr_vxc')
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
         allocate (grhoraux(1, 1, 1))
         allocate (haux(1, 1, 1))
         grhoraux = 0.0_dp
         !
      end if
      !
      !
      allocate (vxc0(nnrx, 2))
      allocate (vxcref(nnrx, 2))
      !
      etxcref = 0.0_dp
      vxcref = 0.0_dp
      !
      !rhoraux = rhoref
      !
      if (dft_is_gradient()) then
         !
         grhoraux(:, :, 1:2) = grhobar(:, :, 1:2)
         grhoraux(:, :, ispin) = grhobar(:, :, ispin) &
                                 + fref*orb_grhor(:, :, 1)
         !
         rhogaux(:, 1:2) = rhobarg(:, 1:2)
         rhogaux(:, ispin) = rhobarg(:, ispin) + fref*orb_rhog(:, 1)

      end if
      !
      !call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
      vxcref = rhoref
      CALL exch_corr_cp(nnrx, 2, grhoraux, vxcref, etxcref) !proposed:giovanni fixing PBE, warning, rhoref overwritten with vxcref, check array dimensions
      !NB grhoaux(nnr,3,nspin)? yes; rhoref(nnr,nspin)? yes
!begin_added:giovanni fixing PBE potential
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
         call gradh(2, grhoraux, rhogaux, vxcref, dexc_dummy, lgam)
         !  grhoraux(nnr,3,nspin)?yes; rhogaux(ng,nspin)? rhoref(nnr, nspin)
         !
      end if
      !

!end_added:giovanni fixing PBE potential

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
         allocate (rhoraux(nnrx, 2))
         !
         rhoraux = rhobar + f*rhoele
         !
         if (dft_is_gradient()) then
            !
            grhoraux(:, :, 1:2) = grhobar(:, :, 1:2)
            grhoraux(:, :, ispin) = grhobar(:, :, ispin) &
                                    + f*orb_grhor(:, :, 1)
            !
            rhogaux(:, 1:2) = rhobarg(:, 1:2)
            rhogaux(:, ispin) = rhobarg(:, ispin) + f*orb_rhog(:, 1)

         end if
         !
         !call exch_corr_wrapper(nnrx,2,grhoraux,rhoraux,etxc,vxc,haux)
         vxc = rhoraux
         CALL exch_corr_cp(nnrx, 2, grhoraux, vxc, etxc) !proposed:giovanni warning rhoraix is overwritten with vxc, check array dimensions
         !NB grhoraux(nnr,3,nspin)? rhoraux(nnr,nspin)?
         !begin_added:giovanni fixing PBE potential
         if (dft_is_gradient()) then
            !
            !  Add second part of the xc-potential to rhor
            !  Compute contribution to the stress dexc
            !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
            call gradh(2, grhoraux, rhogaux, vxc, dexc_dummy, lgam)
            !  grhoraux(nnr,3,nspin)? rhogaux(ng,nspin)? rhoraux(nnr, nspin)
            !
         end if
         !end_added:giovanni fixing PBE potential
         !
         deallocate (rhoraux)
         !
      end if
      !
      deallocate (rhogaux)
      deallocate (orb_rhog)
      !
      etxc0 = 0.0_dp
      vxc0 = 0.0_dp
      !
      !rhoraux = rhobar
      !
      !call exch_corr_wrapper(nnrx,2,grhobar,rhobar,etxc0,vxc0,haux)
      vxc0 = rhobar
      CALL exch_corr_cp(nnrx, 2, grhobar, vxc0, etxc0) !proposed:giovanni warning rhobar is overwritten with vxc0, check array dimensions
      !NB grhobar(nnr,3,nspin)? rhobar(nnr,nspin)?
!begin_added:giovanni fixing PBE potential
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
         call gradh(2, grhobar, rhobarg, vxc0, dexc_dummy, lgam)
         !  grhobar(nnr,3,nspin)? rhogbar(ng,nspin)? rhor(nnr, nspin)
         !
      end if
!end_added:giovanni fixing PBE potential
      !
      ! update vsic_realspace pot
      !
      vsic_realspace(1:nnrx) = vsic_realspace(1:nnrx) &
                     + vxcref(1:nnrx, ispin) - vxc(1:nnrx, ispin)
      !
      ! define pink
      !
      etmp = f*sum(vxcref(1:nnrx, ispin)*rhoele(1:nnrx, ispin))
      !
      pink = (etxc0 - etxc) + etmp + ehele
      pink = pink*fact
      !
      call mp_sum(pink, intra_image_comm)
      !
      call stop_clock('nk_corr_vxc')

      !
      !   calculate wref and wxd
      !
      CALL start_clock('nk_corr_fxc')
      !
      wxdsic_realspace(:, :) = 0.0d0
      !
      if (do_wref .or. do_wxd_) then
         !
         ! note that vxd and wref are updated
         ! (and not overwritten) by the next call
         !
         call nksic_dmxc_spin_cp_update(nnrx, rhoref, f, ispin, rhoele, &
                                        vanishing_rho_w, wrefsic, wxdsic_realspace) !modified:linh
         !
         !
         if (do_wref) then
            !
            w2cst = sum(wrefsic(1:nnrx)*rhoele(1:nnrx, ispin))*fact
            !
            call mp_sum(w2cst, intra_image_comm)
            !
            do ir = 1, nnrx
               wrefsic(ir) = fref*(wrefsic(ir) - w2cst)
            end do
            !
         end if
         !
         if (do_wxd_) then
            !
            wxdsic_realspace(:, 1:2) = rhobarfact*(wxdsic_realspace(:, 1:2) &
                                         + vxc0(:, 1:2) - vxc(:, 1:2))
            !
         end if
         !
      end if
      !
      CALL stop_clock('nk_corr_fxc')

      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink*nkscalfact
      vsic_realspace = vsic_realspace*nkscalfact
      !
      if (do_wxd_) then
         wxdsic_realspace = wxdsic_realspace*nkscalfact
      else
         wxdsic_realspace = 0.d0
      end if
      !
      if (do_wref) then
         wrefsic = wrefsic*nkscalfact
      else
         wrefsic = 0.d0
      end if

      !
      deallocate (vxc0)
      deallocate (vxcref)
      deallocate (rhoele)
      !
      deallocate (grhoraux)
      deallocate (haux)
      !
      if (allocated(orb_grhor)) deallocate (orb_grhor)
      !
      CALL stop_clock('nk_corr')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_nk
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_pz(f, ispin, orb_rhor, &
                                  vsic_realspace, vsic_reciprocal, pink, pzalpha, ibnd, shart)
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
      real(dp), intent(out) :: vsic_realspace(nnrx)
      complex(dp), intent(out) :: vsic_reciprocal(ngm)
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
      vsic_realspace = 0.0_dp
      vsic_reciprocal = 0.0_dp
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
         !pink = pink * pzalpha
         !vsic_realspace = vsic_realspace * pzalpha
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
      ! Store the vsic_realspace potential in reciprocal space
      vtmp = vsic_realspace(:)
      call fwfft('Dense', vtmp, dfftp)
      call psi2rho('Dense', vtmp, dfftp%nnr, vsic_reciprocal, ngm)
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_pz
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_nkpz(f, orb_rhor, vsic_realspace, wrefsic, pink, ibnd, ispin)
!---------------------------------------------------------------
!
! ... calculate the non-Koopmans potential on top of Perdew-Zunger,
!     from the orbital densities
!
      use kinds, only: dp
      use constants, only: e2, fpi
      use cell_base, only: tpiba2, omega
      use nksic, only: fref, nkscalfact, &
                       do_wref, vanishing_rho_w
      use grid_dimensions, only: nnrx, nr1, nr2, nr3
      use gvecp, only: ngm
      use recvecs_indexes, only: np, nm
      use reciprocal_vectors, only: gstart, g
      use eecp_mod, only: do_comp
      use cp_interfaces, only: fwfft, invfft, fillgrad
      use fft_base, only: dfftp
      use funct, only: dmxc_spin, dft_is_gradient
      use mp_global, only: intra_image_comm
      use mp, only: mp_sum
      use control_flags, only: gamma_only, do_wf_cmplx

      !
      implicit none
      real(dp), intent(in)  :: f, orb_rhor(nnrx)
      integer, intent(in)  :: ispin, ibnd
      real(dp), intent(out) :: vsic_realspace(nnrx), wrefsic(nnrx)
      real(dp), intent(out) :: pink
      !
      integer     :: ig, ir
      real(dp)    :: fact, etxcref
      real(dp)    :: w2cst
      !
      real(dp), allocatable :: rhoele(:, :)
      real(dp), allocatable :: rhoref(:, :)
      real(dp), allocatable :: vxcref(:, :)
      real(dp), allocatable :: wxdsic_realspace(:, :)
      real(dp), allocatable :: grhoraux(:, :, :)
      real(dp), allocatable :: haux(:, :, :)
      complex(dp), allocatable :: vhaux(:)
      complex(dp), allocatable :: vcorr(:)
      complex(dp), allocatable :: rhogaux(:, :)
      complex(dp), allocatable :: vtmp(:)
      logical :: lgam
      real(dp) :: dexc_dummy(3, 3)
      !
      CALL start_clock('nk_corr')
      CALL start_clock('nk_corr_h')
      !
      lgam = gamma_only .and. .not. do_wf_cmplx
      fact = omega/DBLE(nr1*nr2*nr3)
      !
      allocate (wxdsic_realspace(nnrx, 2))
      allocate (rhoele(nnrx, 2))
      allocate (rhoref(nnrx, 2))
      allocate (rhogaux(ngm, 2))
      allocate (vtmp(ngm))
      allocate (vcorr(ngm))
      allocate (vhaux(nnrx))
      !
      rhoele = 0.0d0
      rhoele(:, ispin) = orb_rhor(:)
      !
      vsic_realspace = 0.0_dp
      wrefsic = 0.0_dp
      wxdsic_realspace = 0.0_dp
      pink = 0.0_dp
      !
      ! compute self-hartree contributions
      !
      rhogaux = 0.0_dp
      !
      ! rhoele has no occupation
      !
      vhaux(:) = rhoele(:, ispin)
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
!       IF(lgam) THEN  !!!### uncomment for k points
      do ig = 1, ngm
         !
         vhaux(np(ig)) = vtmp(ig)
         vhaux(nm(ig)) = CONJG(vtmp(ig))
         !
      end do
!       ELSE !!!### uncomment for k points
!           do ig=1,ngm !!!### uncomment for k points
      !
!             vhaux(np(ig)) = vtmp(ig) !!!### uncomment for k points
!             vhaux(nm(ig)) = conjg(vtmp(ig))
      !
!           enddo !!!### uncomment for k points
!       ENDIF !!!### uncomment for k points
      !
      call invfft('Dense', vhaux, dfftp)
      !
      ! init here wref sic to save some memory
      !
      ! this is just the self-hartree potential
      ! (to be multiplied by fref later on)
      !
      wrefsic(1:nnrx) = DBLE(vhaux(1:nnrx))
      !
      ! the term - fref has to be included explicitly in rhoele
      !
      vsic_realspace(1:nnrx) = -fref*DBLE(vhaux(1:nnrx))
      !
      deallocate (vtmp)
      deallocate (vcorr)
      deallocate (vhaux)
      !
      call stop_clock('nk_corr_h')
      call start_clock('nk_corr_vxc')
      !
      !   add self-xc contributions
      !
      rhoref = fref*rhoele
      !
      if (dft_is_gradient()) then
         allocate (grhoraux(nnrx, 3, 2))
         allocate (haux(nnrx, 2, 2))
         !
         grhoraux = 0.0_dp
         call fillgrad(1, rhogaux, grhoraux(:, :, ispin:ispin), lgam)
         !
         grhoraux(:, :, ispin) = grhoraux(:, :, ispin)*fref
      else
         allocate (grhoraux(1, 1, 1))
         allocate (haux(1, 1, 1))
         grhoraux = 0.0_dp
      end if
      !

      allocate (vxcref(nnrx, 2))
      !
      etxcref = 0.0_dp
      vxcref = 0.0_dp
      !
      vxcref = rhoref
      !
      CALL exch_corr_cp(nnrx, 2, grhoraux, vxcref, etxcref) !proposed:giovanni fixing PBE, warning, rhoref overwritten with vxcref, check array dimensions
      !
      !begin_added:giovanni fixing PBE potential
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !  Need a dummy dexc here, need to cross-check gradh! dexc should be dexc(3,3), is lgam a variable here?
         call gradh(2, grhoraux, rhogaux, vxcref, dexc_dummy, lgam)
         !  grhoraux(nnr,3,nspin)?yes; rhogaux(ng,nspin)? rhoref(nnr, nspin)
         !
      end if
!end_added:giovanni fixing PBE potential
      deallocate (rhogaux)
!       call exch_corr_wrapper(nnrx,2,grhoraux,rhoref,etxcref,vxcref,haux)
      !
      ! update vsic_realspace pot
      !
      vsic_realspace(1:nnrx) = vsic_realspace(1:nnrx) - vxcref(1:nnrx, ispin)
      !
      ! define pink
      !
      pink = f*sum(vsic_realspace(1:nnrx)*rhoele(1:nnrx, ispin))*fact
      call mp_sum(pink, intra_image_comm)
      !
      call stop_clock('nk_corr_vxc')
      !
      !   calculate wref
      !
      CALL start_clock('nk_corr_fxc')
      !
      if (do_wref) then
         !
         ! note that wxd and wref are updated
         ! (and not overwritten) by the next call
         !
         call nksic_dmxc_spin_cp_update(nnrx, rhoref, f, ispin, rhoele, &
                                        vanishing_rho_w, wrefsic, wxdsic_realspace)!modified:linh
         !
         w2cst = sum(wrefsic(1:nnrx)*rhoele(1:nnrx, ispin))*fact
         !
         call mp_sum(w2cst, intra_image_comm)
         !
         do ir = 1, nnrx
            wrefsic(ir) = -fref*(wrefsic(ir) - w2cst)
         end do
         !
      end if
      !
      CALL stop_clock('nk_corr_fxc')
      !
      !   rescale contributions with the nkscalfact parameter
      !   take care of non-variational formulations
      !
      pink = pink*nkscalfact
      vsic_realspace = vsic_realspace*nkscalfact
      !
      if (do_wref) then
         wrefsic = wrefsic*nkscalfact
      else
         wrefsic = 0.d0
      end if
      !
      deallocate (wxdsic_realspace)
      deallocate (vxcref)
      deallocate (rhoele)
      deallocate (rhoref)
      deallocate (grhoraux)
      deallocate (haux)
      !
      CALL stop_clock('nk_corr')
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_nkpz
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_nkipz(f, ispin, orb_rhor, &
                                     vsic_realspace, vsic_reciprocal, pink, ibnd, shart, is_empty)
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
      real(dp), intent(out) :: vsic_realspace(nnrx)
      complex(dp), intent(out) :: vsic_reciprocal(ngm)
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
      vsic_realspace = 0.0_dp
      vsic_reciprocal = 0.0_dp
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
      ! Store the vsic_realspace potential in reciprocal space
      vtmp = vsic_realspace(:)
      call fwfft('Dense', vtmp, dfftp)
      call psi2rho('Dense', vtmp, dfftp%nnr, vsic_reciprocal, ngm)
      !
      return
      !
!---------------------------------------------------------------
   end subroutine nksic_correction_nkipz
!---------------------------------------------------------------

!---------------------------------------------------------------
   subroutine nksic_correction_nki(f, ispin, orb_rhor, rhor, &
                                   rhoref, rhobar, rhobarg, grhobar, &
                                   vsic_realspace, vsic_reciprocal, wxdsic_realspace, &
                                   wxdsic_reciprocal, do_wxd_, pink, ibnd, shart, is_empty)
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
      real(dp), intent(out) :: vsic_realspace(nnrx)
      complex(dp), intent(out) :: vsic_reciprocal(ngm)
      real(dp), intent(out) :: wxdsic_realspace(nnrx, 2)
      complex(dp), intent(out) :: wxdsic_reciprocal(ngm, 2)
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
      complex(dp), allocatable :: wxdtmp(:, :)
      !
      real(dp), allocatable :: grhoraux(:, :, :)
      real(dp), allocatable :: orb_grhor(:, :, :)
      complex(dp), allocatable :: orb_rhog(:, :)
      real(dp), allocatable :: haux(:, :, :)
      logical :: lgam, is_empty_
      real(dp) :: icoeff
      real(dp) :: dexc_dummy(3, 3)
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
      vsic_realspace = 0.0_dp
      vsic_reciprocal = 0.0_dp
      wxdsic_realspace = 0.0_dp
      wxdsic_reciprocal = 0.0_dp
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
      wxdsic_realspace(:, :) = 0.0d0
      wxdsic_reciprocal(:, :) = 0.0d0
      !
      if (do_wxd_) then
         !
         wxdsic_realspace(:, 1:2) = (1.0_dp - f)*vxc0(:, 1:2) - vxc(:, 1:2) + f*vxcref(:, 1:2)
         !
      end if
      !
      ! Transform wxdsic_realspace to reciprocal space
      !
      allocate (wxdtmp(ngm, 2))
      wxdtmp = wxdsic_realspace(:, :)
      do i = 1, 2
         call fwfft('Dense', wxdtmp(:, i), dfftp)
         call psi2rho('Dense', wxdtmp(:, i), dfftp%nnr, wxdsic_reciprocal(:, i), ngm)
      end do
      deallocate (wxdtmp)
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
         wxdsic_realspace = wxdsic_realspace*nkscalfact
         wxdsic_reciprocal = wxdsic_reciprocal*nkscalfact
         !
      else
         !
         wxdsic_realspace = 0.d0
         wxdsic_reciprocal = 0.d0
         !
      end if
      !
      ! Store the vsic_realspace potential in reciprocal space
      vtmp = vsic_realspace(:)
      call fwfft('Dense', vtmp, dfftp)
      call psi2rho('Dense', vtmp, dfftp%nnr, vsic_reciprocal, ngm)
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

!---------------------------------------------------------------
   subroutine nksic_dmxc_spin_cp_update(nnrx, rhoref, f, ispin, rhoele, &
                                        small, wref, wxd)
!---------------------------------------------------------------

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
      real(dp) :: rhoup, rhodw, rhotot, zeta
      real(dp) :: dmuxc(2, 2)
      real(dp) :: rs, ex, vx, dr, dz, ec, &
                  vcupm, vcdwm, vcupp, vcdwp, &
                  vxupm, vxdwm, vxupp, vxdwp, &
                  dzm, dzp, fact
      !
      real(dp), external :: dpz, dpz_polarized
      integer :: ir
      !logical :: do_exch, do_corr
      !
      real(dp), parameter :: e2 = 2.0_dp, &
                             pi34 = 0.6203504908994_DP, & ! redefined to pi34=(3/4pi)^(1/3)
                             pi34_old = 0.75_dp/3.141592653589793_dp, third = 1.0_dp/3.0_dp, &
                             p43 = 4.0_dp/3.0_dp, p49 = 4.0_dp/9.0_dp, m23 = -2.0_dp/3.0_dp
      !
      if (get_iexch() == 1 .and. get_icorr() == 1) THEN
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
         !
      else
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

            dr = min(1.e-6_dp, 1.e-4_dp*rhotot)
            fact = 0.5d0/dr

            call xc_spin(rhotot - dr, zeta, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
            call xc_spin(rhotot + dr, zeta, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
            !
            dmuxc(1, 1) = dmuxc(1, 1) + (vxupp + vcupp - vxupm - vcupm)*fact
            dmuxc(1, 2) = dmuxc(1, 2) + (vxupp + vcupp - vxupm - vcupm)*fact
            dmuxc(2, 1) = dmuxc(2, 1) + (vxdwp + vcdwp - vxdwm - vcdwm)*fact
            dmuxc(2, 2) = dmuxc(2, 2) + (vxdwp + vcdwp - vxdwm - vcdwm)*fact
            !
            dz = 1.E-6_DP
            dzp = min(1.0, zeta + dz) - zeta
            dzm = -max(-1.0, zeta - dz) + zeta
            !
            fact = 1.0d0/(rhotot*(dzp + dzm))
            !
            call xc_spin(rhotot, zeta - dzm, ex, ec, vxupm, vxdwm, vcupm, vcdwm)
            call xc_spin(rhotot, zeta + dzp, ex, ec, vxupp, vxdwp, vcupp, vcdwp)
            !
            dmuxc(1, 1) = dmuxc(1, 1) + (vxupp + vcupp - vxupm - vcupm)*(1.0_DP - zeta)*fact
            dmuxc(1, 2) = dmuxc(1, 2) - (vxupp + vcupp - vxupm - vcupm)*(1.0_DP + zeta)*fact
            dmuxc(2, 1) = dmuxc(2, 1) + (vxdwp + vcdwp - vxdwm - vcdwm)*(1.0_DP - zeta)*fact
            dmuxc(2, 2) = dmuxc(2, 2) - (vxdwp + vcdwp - vxdwm - vcdwm)*(1.0_DP + zeta)*fact
            !
            ! add corrections to the nksic potentials
            !
            wxd(ir, 1) = wxd(ir, 1) + dmuxc(1, ispin)*rhoele(ir, ispin)*f
            wxd(ir, 2) = wxd(ir, 2) + dmuxc(2, ispin)*rhoele(ir, ispin)*f
            !
            wref(ir) = wref(ir) + dmuxc(ispin, ispin)*rhoele(ir, ispin)
            !
         end do
         !
      end if

      return

!---------------------------------------------------------------
   end subroutine nksic_dmxc_spin_cp_update
!---------------------------------------------------------------
end module nksic_corrections