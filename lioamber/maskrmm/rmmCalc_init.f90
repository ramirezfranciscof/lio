!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_init( qmpos )
!------------------------------------------------------------------------------!
!
!  This subroutine sets up all the things that depend on the atomic positions
!  like g2g-grid, the density matrix fitting coeficients (for the auxiliary
!  density, in int2) and the saved parts of coulomb integrals (int3mem).
!
!  The atomic positions may be provided as input (in which case are stored
!  inside of "r") or they must already be correctly stored there.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod , only: M, igrid2, nuc, natom, natomc, d, r, atmin &
                       &, rmax, jatc, nnps, nnpp, nnpd, nshell, MEMO &
                       &, kkind, kkinds, cool, cools, nsol, pc
   use faint_cpu77, only: int2, int3mem
   use tmpaux_SCF , only: neighbor_list_2e

   implicit none
   real*8, intent(in), optional :: qmpos(3,natom)

   real*8  :: alf, rexp
   integer :: igpu
   integer :: ii, jj, ti, tj, zij
!
!
!
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmCalc_init unknown')

   if ( present(qmpos) ) then
      do ii=1,natom
         r(ii,1) = qmpos(1,ii)
         r(ii,2) = qmpos(2,ii)
         r(ii,3) = qmpos(3,ii)
      end do
   end if

   if (.false.) then
   do ii=1,natom
      natomc(ii)=0
      do jj=1,natom
         d(ii,jj)=0.0d0
         d(ii,jj)=d(ii,jj)+(r(ii,1)-r(jj,1))**2
         d(ii,jj)=d(ii,jj)+(r(ii,2)-r(jj,2))**2
         d(ii,jj)=d(ii,jj)+(r(ii,3)-r(jj,3))**2
         zij=atmin(ii)+atmin(jj)
         ti=atmin(ii)/zij
         tj=atmin(jj)/zij
         alf=atmin(ii)*tj
         rexp=alf*d(ii,jj)
         if (rexp.lt.rmax) then
            natomc(ii)=natomc(ii)+1
            jatc(natomc(ii),ii)=jj
         endif
      enddo
   enddo

   do ii=nshell(0),1,-1
      nnps(nuc(ii))=ii
   enddo

   do ii=nshell(0)+nshell(1),nshell(0)+1,-1
      nnpp(nuc(ii))=ii
   enddo

   do ii=M,nshell(0)+nshell(1)+1,-1
      nnpd(nuc(ii))=ii
   enddo

   end if

! Nano: calculating neighbour list helps to make 2 electrons integral scale
! linearly with natoms/basis.
   call neighbor_list_2e()
   call g2g_timer_stop('rmmCalc_init unknown')
!
!
!
!------------------------------------------------------------------------------!
! -Create integration grid for XC here
! -Assign points to groups (spheres/cubes)
! -Assign significant functions to groups
! -Calculate point weights
!
! Precalculate two-index (density basis) "G" matrix used in density fitting
! here (S_ij in Dunlap, et al JCP 71(8) 1979) into RMM(M7)
! Also, pre-calculate G^-1 if G is not ill-conditioned into RMM(M9)
!
! Precalculate three-index (two in MO basis, one in density basis) matrix
! used in density fitting / Coulomb F element calculation here
! (t_i in Dunlap)
!
! Large elements of t_i put into double-precision cool here
! Size criteria based on size of pre-factor in Gaussian Product Theorem
! (applied to MO basis indices)
!
   call g2g_timer_start('rmmCalc_init g2g')
   call g2g_reload_atom_positions( igrid2 )
   call g2g_timer_stop('rmmCalc_init g2g')

   call g2g_timer_start('rmmCalc_init densfit')
   call int2()
   call g2g_timer_stop('rmmCalc_init densfit')

   call g2g_timer_start('rmmCalc_init aint')
   call aint_query_gpu_level( igpu )
   if (igpu.gt.1) call aint_new_step()
   if (igpu.gt.1) call aint_qmmm_init( nsol, r, pc)
   if (igpu.gt.2) call aint_coulomb_init()
   call g2g_timer_stop('rmmCalc_init aint')

   if (igpu.eq.5) MEMO = .false.
   if (MEMO) then
      call g2g_timer_start('rmmCalc_init intmem')
      if (allocated(kkind))  deallocate(kkind)
      if (allocated(kkinds)) deallocate(kkinds)
      if (allocated(cool))   deallocate(cool)
      if (allocated(cools))  deallocate(cools)
      call int3mem()
      call g2g_timer_stop('rmmCalc_init intmem')
   endif
!
!
end subroutine rmmCalc_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
