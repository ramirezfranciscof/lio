!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE magnus_sc
     > (M,Rho0,Fock0,FockGuess0,Uinv,Linv,dt,Rho,Fock,Energy)
!------------------------------------------------------------------------------!
!
! OBS: ...Uinv,Linv,... = ...X,Xtrans,...
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY: NBCH
       IMPLICIT NONE
       INTEGER,INTENT(IN)    :: M
       COMPLEX*16,INTENT(IN)  :: Rho0(M,M)
       REAL*8,INTENT(IN)     :: Fock0(M,M),FockGuess0(M,M),
     >                          Linv(M,M),Uinv(M,M)
       REAL*8,INTENT(IN)     :: dt
       COMPLEX*16,INTENT(OUT) :: Rho(M,M)
       REAL*8,INTENT(OUT)    :: Fock(M,M),Energy
!
       INTEGER :: nstep,max_nstep
       REAL*8  :: Ecoul,Exc,Eelec,Efld,dif,tol_dif
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: FockGuess,FockMid,RhoReal
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
! SETUP VARIABLES
!----------------------------------------------------------!
       call g2g_timer_start('magnus_sc')
       Fock=FockGuess0
       allocate(FockGuess(M,M),FockMid(M,M),RhoReal(M,M))
       NBCH=10        !
       tol_dif=1.0E-8 !
       max_nstep=20   !
       nstep=0
       dif=tol_dif+1
!
!
!
! SELF CONSISTENT LOOP
!----------------------------------------------------------!
       DO WHILE ((dif.GE.tol_dif).AND.(nstep.LT.max_nstep))
         nstep=nstep+1
         FockGuess=Fock
         FockMid=(Fock0+FockGuess)/2
         call magnus_bch(M,FockMid,Rho0,Rho,NBCH,dt)

         RhoReal=real(Rho)
         call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoON  => RhoOA
         call rmm_put_rho(M,RhoReal)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,Fock)
         call MatrixTransform(M,Linv,Fock,Uinv)    ! FockOA => FockON
         call matdist(M,FockGuess,Fock,dif)
       ENDDO
       if (dif.gt.tol_dif) write(6,200) 'Convergence not achieved: ',dif

       Energy=Ecoul+Exc+Eelec+Efld
!       Fock=(Fock+FockGuess)/2
!
!
!
! ENDING
!----------------------------------------------------------!
       deallocate(FockGuess,FockMid,RhoReal)
       call g2g_timer_stop('magnus_sc')
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 200   FORMAT(A,F15.8)
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
