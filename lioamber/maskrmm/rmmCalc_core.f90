!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_core( Smat, Hmat, Enn, Ens, do_raw, do_ecp, do_sol )
!------------------------------------------------------------------------------!
!
!  This subroutine calculates all the things that depend exclusively on the
!  atomic positions (already setted up in rmmCalc_init). That is: the overlap
!  matrix, the core fock matrix, and the energy associated to the nuclear
!  interactions with themselves and the solven.
!
!  Notice that you can choose which parts of the core to calculate with each
!  call, so one may receive all contributions in a single matrix or one may
!  call this several times so as to get the different contributions in
!  different matrices (to, for example, obtain discriminated energies).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   use garcha_mod,  only: M, Md, igrid2, RMM

   use faint_cpu77, only: int1, intsol

   use ECP_mod    , only: ecpmode, term1e, VAAA, VAAB, VBAC, &
                       & FOCK_ECP_read, FOCK_ECP_write

   implicit none
   real*8 , intent(inout) :: Smat(M,M)
   real*8 , intent(inout) :: Hmat(M,M)
   logical, intent(in)    :: do_raw
   logical, intent(in)    :: do_ecp
   logical, intent(in)    :: do_sol
   real*8 , intent(out)   :: Enn ! Nuclear-nuclear energy
   real*8 , intent(out)   :: Ens ! Solvent-nuclear energy

   real*8  :: E1s ! Solvent-electron energy with bad dens
   integer :: igpu, idx0, MM, MMd, kk
!
!
!
!------------------------------------------------------------------------------!
   Enn = 0.0d0
   if (do_raw) then
      call g2g_timer_start('rmmCalc_core Hraw')
      call int1( Enn )
      call g2g_timer_stop('rmmCalc_core Hraw')
   else
      call rmmput_core( Hmat )
   end if
!
!
!
!------------------------------------------------------------------------------!
   if ( (ecpmode).and.(do_ecp) ) then
      call g2g_timer_start('rmmCalc_core Hecp')

      if (FOCK_ECP_read) then
!        Variable allocation and data read from ECP_restart 
         call intECP(0)

      else
!        Variable allocation and calculation of the N center terms by different
!        calls to subroutine intECP(N).
         call intECP(1)
         call intECP(2)
         call intECP(3)

      end if

      if (FOCK_ECP_write) call WRITE_ECP()
      call WRITE_POST(1)

      write(*,*) "Modifying Fock Matrix with ECP terms"
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      idx0=3*MM+2*MMd
      do kk = 1, MM
!        Backups 1e terms and modifies fock
         term1e( kk )   = RMM( idx0+kk )
         RMM( idx0+kk ) = RMM( idx0+kk ) + VAAA(kk) + VAAB(kk) + VBAC(kk)
      enddo

      call g2g_timer_stop('rmmCalc_core Hecp')
   end if
!
!
!
!------------------------------------------------------------------------------!
   Ens = 0.0d0
   call aint_query_gpu_level(igpu)
   if (do_sol) then
      if (igpu.le.1) then
         call g2g_timer_start('rmmCalc_core Hsol fort')
         call intsol( E1s, Ens, .true. )
         call g2g_timer_stop('rmmCalc_core Hsol fort')
      else
         call g2g_timer_start('rmmCalc_core Hsol aint')
         call aint_qmmm_fock( E1s, Ens )
         call g2g_timer_stop('rmmCalc_core Hsol aint')
      end if
   end if
!
!
!
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmCalc_core exit')
   call rmmget_fock( Smat )
   call rmmget_core( Hmat )
   call g2g_timer_stop('rmmCalc_core exit')
!
!
end subroutine rmmCalc_core
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
