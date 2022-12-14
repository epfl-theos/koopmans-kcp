!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!====================================================================
   SUBROUTINE inner_loop_cold( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter, vpot )
!====================================================================
      !
      ! minimizes the total free energy
      ! using cold smearing,
      !
      !

      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, taurdr, &
                                tprnfor, ndr, ndw, nbeg, nomore, &
                                tsde, tortho, tnosee, tnosep, trane, &
                                tranp, tsdp, tcp, tcap, ampre, &
                                amprp, tnoseh, gamma_only, do_wf_cmplx !added:giovanni gamma_only do_wf_cmplx
      USE core,           ONLY: nlcc_any
      USE energies,       ONLY: eht, epseu, exc, etot, eself, enl, &
                                ekin, atot, entropy, egrand
      USE electrons_base, ONLY: f, nspin, nel, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE ensemble_dft,   ONLY: tens,  ninner, ismear, etemp, &
                                ef, z0t, c0diag, becdiag, &
                                e0, psihpsi, compute_entropy2, &
                                compute_entropy_der, compute_entropy, &
                                niter_cold_restart, lambda_cold
      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: gstart
      USE cvan,           ONLY: nvb, ish
      USE ions_base,      ONLY: na, nat, pmass, nax, nsp, rcmax
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE cell_base,      ONLY: ainv, a1, a2, a3
      USE cell_base,      ONLY: omega, alat
      USE cell_base,      ONLY: h, hold, deth, wmass, tpiba2
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx, nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, &
                          ONLY: nnrb => nnrbx, nr1b, nr2b, nr3b
      USE local_pseudo,   ONLY: vps, rhops
      USE io_global,      ONLY: io_global_start, stdout, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm, leg_ortho
      USE dener
      !USE derho
      USE cdvan
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE cg_module,      ONLY: ene_ok
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast, mp_root_sum

      USE cp_interfaces,  ONLY: rhoofr, dforce
      USE cg_module,      ONLY: itercg
      USE cp_main_variables, ONLY: distribute_lambda, descla, nlax, collect_lambda
      USE descriptors,       ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                   descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                   nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , la_me_ , la_nrl_
      USE dspev_module,   ONLY: pdspev_drv, dspev_drv
      USE twin_types !added:giovanni


      !
      IMPLICIT NONE

!input variables
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(kind=DP)            :: eigr( ngw, nat )
      COMPLEX(kind=DP)            :: c0( ngw, n )
      TYPE(twin_matrix)               :: bec!( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (kind=DP)           :: eigrb( ngb, nat )
      REAL(kind=DP)               :: rhor( nnr, nspin )
      REAL(kind=DP)               :: vpot( nnr, nspin )
      COMPLEX(kind=DP)            :: rhog( ngm, nspin )
      REAL(kind=DP)               :: rhos( nnrsx, nspin )
      REAL(kind=DP)               :: rhoc( nnr )
      COMPLEX(kind=DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(kind=DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(kind=DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(kind=DP)            :: sfac( ngs, nsp )  

!local variables
      REAL(kind=DP) :: atot0, atot1, atotl, atotmin
      REAL(kind=DP), ALLOCATABLE :: fion2(:,:)
      type(twin_matrix), dimension(:), allocatable :: c0hc0(:)!modified:giovanni
      REAL(kind=DP), ALLOCATABLE :: mtmp(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: mtmp_c(:,:)
      COMPLEX(kind=DP), ALLOCATABLE :: h0c0(:,:)
      INTEGER :: niter
      INTEGER :: i,k, is, nss, istart, ig, iss
      REAL(kind=DP) :: lambda, lambdap
      REAL(kind=DP), ALLOCATABLE :: epsi0(:,:)

      INTEGER :: np(2), coor_ip(2), ipr, ipc, nr, nc, ir, ic, ii, jj, root, j
      INTEGER :: desc_ip( descla_siz_ )
      INTEGER :: np_rot, me_rot, comm_rot, nrl
      !
      LOGICAL :: lgam !added:giovanni

      CALL start_clock( 'inner_loop')
       
      lgam=gamma_only.and..not.do_wf_cmplx !added:giovanni
      !
      allocate(fion2(3,nat))
!       allocate(c0hc0(nlax,nlax,nspin))
      allocate(h0c0(ngw,nx))
      !begin_added:giovanni
      allocate(c0hc0(nspin))
      do j=1,nspin
          call init_twin(c0hc0(j), lgam)
          call allocate_twin(c0hc0(j),nlax,nlax, lgam)
      enddo
      !end_added:giovanni
      lambdap=0.3d0!small step for free-energy calculation

      ! calculates the initial free energy if necessary
      IF( .not. ene_ok ) THEN

        ! calculates the overlaps bec between the wavefunctions c0
        ! and the beta functions
        CALL calbec( 1, nsp, eigr, c0, bec )
 
        ! rotates the wavefunctions c0 and the overlaps bec
        ! (the occupation matrix f_ij becomes diagonal f_i)      

        CALL rotate_twin( z0t, c0, bec, c0diag, becdiag, .false. )
  
        ! calculates the electronic charge density
        CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 )
        IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
        ! calculates the SCF potential, the total energy
        ! and the ionic forces
        vpot = rhor
        CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
       !entropy value already  been calculated

      END IF
     
      atot0=etot+entropy
      
    
!starts inner loop
      do niter=1,ninner
!calculates c0hc0, which defines the search line (1-\labda)* psihpsi+\labda*c0hc0

      ! calculateas the energy contribution associated with 
      ! the augmentation charges and the 
      ! corresponding contribution to the ionic force
       
         CALL newd( vpot, irb, eigrb, rhovan, fion2 )

         ! operates the Hamiltonian on the wavefunction c0
         h0c0( :, : )= 0.D0
         DO i= 1, n, 2                      
            CALL dforce( i, bec, betae, c0, h0c0(:,i), h0c0(:,i+1), rhos, nnrsx, ispin, f, n, nspin )
         END DO

    
          
         ! calculates the Hamiltonian matrix in the basis {c0}           
!          c0hc0(:,:,:)=0.d0
         do i=1,nspin
	    call set_twin(c0hc0(i), CMPLX(0.d0,0.d0))
         enddo
         !
         DO is= 1, nspin

            nss= nupdwn( is )
            istart= iupdwn( is )

            np(1) = descla( la_npr_ , is )
            np(2) = descla( la_npc_ , is )

            DO ipc = 1, np(2)
               DO ipr = 1, np(1)

                  coor_ip(1) = ipr - 1
                  coor_ip(2) = ipc - 1
                  CALL descla_init( desc_ip, descla( la_n_ , is ), descla( la_nx_ , is ), np, coor_ip, descla( la_comm_ , is ), 1 )

                  nr = desc_ip( nlar_ )
                  nc = desc_ip( nlac_ )
                  ir = desc_ip( ilar_ )
                  ic = desc_ip( ilac_ )

                  CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                    desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )
                  !
                  root = root * leg_ortho

                  IF(.not.c0hc0(iss)%iscmplx) THEN
		    ALLOCATE( mtmp( nr, nc ) )
		    mtmp = 0.0d0
		    CALL DGEMM( 'T', 'N', nr, nc, 2*ngw, - 2.0d0, c0( 1, istart + ir - 1 ), 2*ngw, &
                              h0c0( 1, istart + ic - 1 ), 2*ngw, 0.0d0, mtmp, nr )
		    IF (gstart == 2) THEN
		      DO jj = 1, nc
			  DO ii = 1, nr
			    i = ii + ir - 1
			    j = jj + ic - 1
			    mtmp(ii,jj) = mtmp(ii,jj) + DBLE( c0( 1, i + istart - 1 ) ) * DBLE( h0c0( 1, j + istart - 1 ) )
			  END DO
		      END DO
		    END IF
                  ELSE
		    ALLOCATE( mtmp_c( nr, nc ) )
		    mtmp_c = CMPLX(0.0d0,0.d0)
		    CALL ZGEMM( 'C', 'N', nr, nc, ngw, CMPLX(- 1.0d0,0.d0), c0( 1, istart + ir - 1 ),ngw, &
                              h0c0( 1, istart + ic - 1 ), ngw, CMPLX(0.0d0,0.d0), mtmp_c, nr )
                  ENDIF
                  
                  IF(.not.c0hc0(is)%iscmplx) THEN
		    CALL mp_root_sum( mtmp, c0hc0(is)%rvec(1:nr,1:nc), root, intra_image_comm )
		    DEALLOCATE( mtmp )
                  ELSE
		    CALL mp_root_sum( mtmp_c, c0hc0(is)%cvec(1:nr,1:nc), root, intra_image_comm )
		    DEALLOCATE( mtmp_c )
                  ENDIF
!                  IF( coor_ip(1) == descla( la_myr_ , is ) .AND. &
!                      coor_ip(2) == descla( la_myc_ , is ) .AND. descla( lambda_node_ , is ) > 0 ) THEN
!                     c0hc0(1:nr,1:nc,is) = mtmp
!                  END IF
               END DO
            END DO
         END DO


         if(mod(itercg,niter_cold_restart) == 0) then
!calculates free energy at lamda=1.
            CALL inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                 rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                 sfac, c0, bec, firstiter,psihpsi,c0hc0,1.d0,atot1, vpot)
!calculates free energy at lamda=lambdap
       
            CALL inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                 rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                 sfac, c0, bec, firstiter,psihpsi,c0hc0,lambdap,atotl, vpot)
!find minimum point lambda
        
            CALL three_point_min(atot0,atotl,atot1,lambdap,lambda,atotmin)
        
         else
            atotl=atot0
            atot1=atot0
            lambda=lambda_cold
         endif

!calculates free energy and rho at lambda
        

         ! calculates the new matrix psihpsi
         
         IF(lgam) THEN
	    DO is= 1, nspin
		psihpsi(is)%rvec(:,:) = (1.d0-lambda) * psihpsi(is)%rvec(:,:) + lambda * c0hc0(is)%rvec(:,:)
	    END DO
         ELSE
	    DO is= 1, nspin
		psihpsi(is)%cvec(:,:) = (1.d0-lambda) * psihpsi(is)%cvec(:,:) + lambda * c0hc0(is)%cvec(:,:)
	    END DO
         ENDIF

         ! diagonalize and calculates energies
         
         e0( : )= 0.D0 

         CALL  inner_loop_diag( c0, bec, psihpsi, z0t, e0 )

         !calculates fro e0 the new occupation and the entropy        

         CALL efermi( nelt, n, etemp, 1, f, ef, e0, entropy, ismear, nspin )


        !calculates new charge and new energy
        

        ! calculates the electronic charge density
         CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                     rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
        ! calculates the SCF potential, the total energy
        ! and the ionic forces
         vpot = rhor
         CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
                     tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                     tau0, fion2 )
       !entropy value already  been calculated
         if(ionode) then
           write(37,*) niter
           write(37,*) atot0,atotl,atot1
           write(37,*) lambda,atotmin,etot+entropy
         endif
         atotl=atot0
         atot0=etot+entropy

   
         if(lambda==0.d0) exit 


         
      enddo


      atot=atot0

!ATTENZIONE
!the following is of capital importance
      ene_ok= .TRUE.
!but it would be better to avoid it
      do is=1,nspin
       call deallocate_twin(c0hc0(is))
      enddo

      DEALLOCATE(fion2)
      DEALLOCATE(c0hc0)
      DEALLOCATE(h0c0)

      CALL stop_clock( 'inner_loop' )
      return
!====================================================================      
    END SUBROUTINE inner_loop_cold
!====================================================================


   SUBROUTINE inner_loop_lambda( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                          rhor, rhog, rhos, rhoc, ei1, ei2, ei3, &
                          sfac, c0, bec, firstiter,c0hc0,c1hc1,lambda,  &
                          free_energy, vpot )
    
!this subroutine for the energy matrix (1-lambda)c0hc0+labda*c1hc1
!calculates the corresponding free energy


      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, taurdr, &
                                tprnfor, ndr, ndw, nbeg, nomore, &
                                tsde, tortho, tnosee, tnosep, trane, &
                                tranp, tsdp, tcp, tcap, ampre, &
                                amprp, tnoseh
      USE core,           ONLY: nlcc_any
      USE energies,       ONLY: eht, epseu, exc, etot, eself, enl, &
                                ekin, atot, entropy, egrand
      USE electrons_base, ONLY: f, nspin, nel, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE ensemble_dft,   ONLY: tens,  ninner, ismear, etemp, &
                                 c0diag, becdiag
      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: gstart
      USE cvan,           ONLY: nvb, ish
      USE ions_base,      ONLY: na, nat, pmass, nax, nsp, rcmax
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE cell_base,      ONLY: ainv, a1, a2, a3
      USE cell_base,      ONLY: omega, alat
      USE cell_base,      ONLY: h, hold, deth, wmass, tpiba2
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx, nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, &
                          ONLY: nnrb => nnrbx, nr1b, nr2b, nr3b
      USE local_pseudo,   ONLY: vps, rhops
      USE io_global,      ONLY: io_global_start, stdout, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm
      USE dener
      !USE derho
      USE cdvan
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE cg_module,      ONLY: ene_ok
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast
      use cp_interfaces,  only: rhoofr, dforce
      USE cp_main_variables, ONLY: descla, nlax, nrlx

      !
      IMPLICIT NONE

!input variables
      INTEGER                :: nfi
      LOGICAL                :: tfirst 
      LOGICAL                :: tlast
      COMPLEX(kind=DP)            :: eigr( ngw, nat )
      COMPLEX(kind=DP)            :: c0( ngw, n )
      REAL(kind=DP)               :: bec( nhsa, n )
      LOGICAL                :: firstiter


      INTEGER                :: irb( 3, nat )
      COMPLEX (kind=DP)           :: eigrb( ngb, nat )
      REAL(kind=DP)               :: rhor( nnr, nspin )
      REAL(kind=DP)               :: vpot( nnr, nspin )
      COMPLEX(kind=DP)            :: rhog( ngm, nspin )
      REAL(kind=DP)               :: rhos( nnrsx, nspin )
      REAL(kind=DP)               :: rhoc( nnr )
      COMPLEX(kind=DP)            :: ei1( nr1:nr1, nat )
      COMPLEX(kind=DP)            :: ei2( nr2:nr2, nat )
      COMPLEX(kind=DP)            :: ei3( nr3:nr3, nat )
      COMPLEX(kind=DP)            :: sfac( ngs, nsp )
  
      REAL(kind=DP), INTENT(in)   :: c0hc0(nlax,nlax,nspin)
      REAL(kind=DP), INTENT(in)   :: c1hc1(nlax,nlax,nspin)
      REAL(kind=DP), INTENT(in)   :: lambda
      REAL(kind=DP), INTENT(out)  :: free_energy


!local variables
      REAL(kind=DP), ALLOCATABLE  :: clhcl(:,:,:), fion2(:,:)
      INTEGER :: i,k, is, nss, istart, ig
      REAL(kind=DP), ALLOCATABLE :: eaux(:), faux(:), zauxt(:,:,:)
      REAL(kind=DP) :: entropy_aux, ef_aux

      CALL start_clock( 'inner_lambda')
      
      allocate(clhcl(nlax,nlax,nspin))
      allocate(eaux(nx))
      allocate(faux(nx))
      allocate(zauxt(nrlx,nudx,nspin))
      allocate(fion2(3,nat))


!calculates the matrix clhcl
         
      DO is= 1, nspin
         clhcl(:,:,is)=(1.d0-lambda)*c0hc0(:,:,is)+lambda*c1hc1(:,:,is)
      END DO

      CALL inner_loop_diag( c0, bec, clhcl, zauxt, eaux )
      
      faux(:)=f(:)

      CALL efermi( nelt, n, etemp, 1, f, ef_aux, eaux, entropy_aux, ismear, nspin )


      ! calculates the electronic charge density
      CALL rhoofr( nfi, c0diag, irb, eigrb, becdiag, rhovan, &
                   rhor, rhog, rhos, enl, denl, ekin, dekin6 )
      IF(nlcc_any) CALL set_cc( irb, eigrb, rhoc )
  
      ! calculates the SCF potential, the total energy
      ! and the ionic forces
      vpot = rhor
      CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, &
                   tlast, ei1, ei2, ei3, irb, eigrb, sfac, &
                   tau0, fion2 )
      !entropy value already  been calculated

      free_energy=etot+entropy_aux
         
      f(:)=faux(:)

      deallocate(clhcl)
      deallocate(faux)
      deallocate(eaux)
      deallocate(zauxt)
      deallocate(fion2)

      CALL stop_clock( 'inner_lambda')

      return

    END SUBROUTINE inner_loop_lambda


    SUBROUTINE three_point_min(y0,yl,y1,l,lambda,amin)
!calculates the estimate for the minimum 

      USE kinds,  ONLY : DP


      implicit none

      REAL(kind=DP), INTENT(in) :: y0,yl,y1, l
      REAL(kind=DP), INTENT(out) :: lambda, amin


      REAL(kind=DP) :: a,b,c, x_min, y_min

! factors for f(x)=ax**2+b*x+c
      c=y0
      b=(yl-y0-y1*l**2.d0+y0*l**2.d0)/(l-l**2.d0)
      a=y1-y0-b
   
      
      x_min=-b/(2.d0*a)
      if( x_min <= 1.d0 .and. x_min >= 0.d0) then
         y_min=a*x_min**2.d0+b*x_min+c
         if(y_min <= y0 .and. y_min <= y1) then
            lambda=x_min
            amin=y_min
         else
            if(y0 < y1) then
                lambda=0.d0
                amin=y0
             else
                lambda=1.d0
                amin=y1
             endif
          endif
       else
          if(y0 < y1) then
             lambda=0.d0
             amin=y0
          else
             lambda=1.d0
             amin=y1
          endif
       endif


       return

     END SUBROUTINE three_point_min


!====================================================================
   SUBROUTINE inner_loop_diag( c0, bec, psihpsi, z0t, e0 )
!====================================================================
      !
      ! minimizes the total free energy
      ! using cold smearing,
      !
      ! declares modules
      USE kinds,          ONLY: dp
      USE control_flags,  ONLY: iprint, thdyn, tpre, iprsta, &
                                tfor, taurdr, &
                                tprnfor, ndr, ndw, nbeg, nomore, &
                                tsde, tortho, tnosee, tnosep, trane, &
                                tranp, tsdp, tcp, tcap, ampre, &
                                amprp, tnoseh
      USE energies,       ONLY: eht, epseu, exc, etot, eself, enl, &
                                ekin, atot, entropy, egrand
      USE electrons_base, ONLY: f, nspin, nel, iupdwn, nupdwn, nudx, &
                                nelt, nx => nbspx, n => nbsp, ispin 

      USE ensemble_dft,   ONLY: tens,  ninner, ismear, etemp, &
                                ef, c0diag, becdiag, &
                                compute_entropy2, &
                                compute_entropy_der, compute_entropy, &
                                niter_cold_restart, lambda_cold
      USE gvecp,          ONLY: ngm
      USE gvecs,          ONLY: ngs
      USE gvecb,          ONLY: ngb
      USE gvecw,          ONLY: ngw
      USE reciprocal_vectors, &
                          ONLY: gstart
      USE cvan,           ONLY: nvb, ish
      USE ions_base,      ONLY: na, nat, pmass, nax, nsp, rcmax
      USE grid_dimensions, &
                          ONLY: nnr => nnrx, nr1, nr2, nr3
      USE cell_base,      ONLY: ainv, a1, a2, a3
      USE cell_base,      ONLY: omega, alat
      USE cell_base,      ONLY: h, hold, deth, wmass, tpiba2
      USE smooth_grid_dimensions, &
                          ONLY: nnrsx, nr1s, nr2s, nr3s
      USE smallbox_grid_dimensions, &
                          ONLY: nnrb => nnrbx, nr1b, nr2b, nr3b
      USE local_pseudo,   ONLY: vps, rhops
      USE io_global,      ONLY: io_global_start, stdout, ionode, &
                                ionode_id
      USE mp_global,      ONLY: intra_image_comm
      USE dener
      !USE derho
      USE cdvan
      USE io_files,       ONLY: psfile, pseudo_dir, outdir
      USE uspp,           ONLY: nhsa=> nkb, betae => vkb, &
                                rhovan => becsum, deeq
      USE uspp_param,     ONLY: nh
      USE cg_module,      ONLY: ene_ok
      USE ions_positions, ONLY: tau0
      USE mp,             ONLY: mp_sum,mp_bcast, mp_root_sum

      USE cp_interfaces,  ONLY: rhoofr, dforce
      USE cg_module,      ONLY: itercg
      USE cp_main_variables, ONLY: distribute_lambda, descla, nlax, collect_lambda, nrlx
      USE descriptors,       ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                   descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                   nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , &
                                   la_me_ , la_nrl_ 
      USE dspev_module,   ONLY: pdspev_drv, dspev_drv


      !
      IMPLICIT NONE

      COMPLEX(kind=DP)            :: c0( ngw, n )
      REAL(kind=DP)               :: bec( nhsa, n )
      REAL(kind=DP)               :: psihpsi( nlax, nlax, nspin )
      REAL(kind=DP)               :: z0t( nrlx, nudx, nspin )
      REAL(kind=DP)               :: e0( nx )

      INTEGER :: i,k, is, nss, istart, ig
      REAL(kind=DP), ALLOCATABLE :: epsi0(:,:), dval(:), zaux(:,:), mtmp(:,:)

      INTEGER :: np(2), coor_ip(2), ipr, ipc, nr, nc, ir, ic, ii, jj, root, j
      INTEGER :: np_rot, me_rot, comm_rot, nrl, kk

      CALL start_clock( 'inner_diag')
      
      e0( : )= 0.D0 

      DO is = 1, nspin
 
            istart   = iupdwn( is )
            nss      = nupdwn( is )
            np_rot   = descla( la_npr_ , is )  * descla( la_npc_ , is )
            me_rot   = descla( la_me_ , is )
            nrl      = descla( la_nrl_ , is )
            comm_rot = descla( la_comm_ , is )

            allocate( dval( nx ) )

            dval = 0.0d0

            IF( descla( lambda_node_ , is ) > 0 ) THEN
               !
               ALLOCATE( epsi0( nrl, nss ), zaux( nrl, nss ) )

               CALL blk2cyc_redist( nss, epsi0, nrl, nss, psihpsi(1,1,is), SIZE(psihpsi,1), SIZE(psihpsi,2), descla(1,is) )

               CALL pdspev_drv( 'V', epsi0, nrl, dval, zaux, nrl, nrl, nss, np_rot, me_rot, comm_rot )
               !
               IF( me_rot /= 0 ) dval = 0.0d0
               !
            ELSE

               ALLOCATE( epsi0( 1, 1 ), zaux( 1, 1 ) )

            END IF
            
            CALL mp_sum( dval, intra_image_comm )

            DO i = 1, nss
               e0( i+istart-1 )= dval( i )
            END DO

            z0t(:,:,is) = 0.0d0

            IF( descla( lambda_node_ , is ) > 0 ) THEN
               !NB zaux is transposed
               !ALLOCATE( mtmp( nudx, nudx ) )
               z0t( 1:nrl , 1:nss, is ) = zaux( 1:nrl, 1:nss )
            END IF

            DEALLOCATE( epsi0 , zaux, dval )

   END DO

   ! rotates the wavefunctions c0 and the overlaps bec
   ! (the occupation matrix f_ij becomes diagonal f_i)      

   CALL rotate ( z0t, c0, bec, c0diag, becdiag, .false. )

   CALL stop_clock( 'inner_diag')
     
   RETURN
END SUBROUTINE
