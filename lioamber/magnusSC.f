!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE magnusSC
     > (M,NBCH,Rho0,Fock0,FockGuess0,Uinv,Linv,dt,Rho,Fock,Energy)
!------------------------------------------------------------------------------!
!
! OBS: ...Uinv,Linv,... = ...X,Xtrans,...
!
!------------------------------------------------------------------------------!
       IMPLICIT NONE
       INTEGER,INTENT(IN)     :: M,NBCH
       COMPLEX*16,INTENT(IN)  :: Rho0(M,M)
       REAL*8,INTENT(IN)      :: Fock0(M,M),FockGuess0(M,M),
     >                           Linv(M,M),Uinv(M,M)
       REAL*8,INTENT(IN)      :: dt
       COMPLEX*8,INTENT(OUT)  :: Rho(M,M)
       REAL*8,INTENT(OUT)     :: Fock(M,M),Energy
!
       LOGICAL :: docycle
       INTEGER :: step,step_max
       REAL*8  :: diff,diff_min,diff_tol
       REAL*8  :: Ecoul,Exc,Eelec,Efld
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: FockGuess,FockMid,RhoReal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
! SETUP
!----------------------------------------------------------!
       call g2g_timer_start('magnusSC')
       allocate(FockGuess(M,M),FockMid(M,M),RhoReal(M,M))
       diff_tol=1.0E-6
       diff_min=1.0E-10
       step_max=50
       step=0

!
!
!
! SELF CONSISTENT LOOP
!----------------------------------------------------------!
       docycle=.true.
       Fock=FockGuess0
       DO WHILE (docycle)
         FockGuess=Fock
         FockMid=(Fock0+FockGuess)/2
         call magnusBCH(M,FockMid,Rho0,Rho,NBCH,dt)

         RhoReal=real(Rho)
         call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoON  => RhoOA
         call rmm_put_rho(M,RhoReal)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,Fock)
         call MatrixTransform(M,Linv,Fock,Uinv)    ! FockOA => FockON

         nstep=nstep+1
         call matdist(M,FockGuess,Fock,diff)
         docycle=(docycle).and.(diff.gt.diff_max)
         docycle=(docycle).and.(step.lt.step_max)
       ENDDO
       Energy=Ecoul+Exc+Eelec+Efld
!
!
!
! CHECK FOR PROBLEMS 
!----------------------------------------------------------!
       if (diff.gt.diff_tol) then
         write(6,'(A)');write(6,'(A)');write(6,'(A)')
         write(6,*) 'CONVERGENCE NOT ACHIEVED:   DIFF = ',diff
         stop
       endif
!
!
! ENDING
!------------------------------------------------------------------------------!
       deallocate(FockGuess,FockMid,RhoReal)
       call g2g_timer_stop('magnusSC')
 200   format(A,F15.8)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
