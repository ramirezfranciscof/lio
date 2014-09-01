!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine magnuspred
     > (M,Fock1a,Fock1b,Rho2,Fock5,NBCH,dt,Umat,Uinv,Lmat,Linv)
!------------------------------------------------------------------------------!
!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))      
!
! References for further information:
!   (*) Phys. Rev. B 2006, 74, 155112
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       integer,intent(in)    :: M,NBCH
       complex*8, intent(in) :: Rho2(M,M)
       real*8,intent(in)     :: Fock1a(M,M),Fock1b(M,M),dt
       real*8,intent(in)     :: Umat(M,M),Uinv(M,M),Lmat(M,M),Linv(M,M)
       real*8,intent(out)    :: Fock5(M,M)

       complex*8,allocatable :: Rho4(:,:)
       real*8,allocatable    :: Fock3(:,:),RhoReal(:,:)
       integer :: i,j,k,kk
       real*8 :: Ecoul,Exc,Eelec,Efld,dtmid
!------------------------------------------------------------------------------!
       call g2g_timer_start('magnuspred')
       ALLOCATE(Rho4(M,M),Fock3(M,M),RhoReal(M,M))
       dtmid=dt*0.50D0

! Step 1: Matrix F1a and F1b are used to extrapolate F3
       Fock3=(7.0d0/4.0d0)*Fock1b-(3.0d0/4.0d0)*Fock1a

! Step2: F3 is used to propagate rho2 to rho4
       call magnus(M,Fock3,Rho2,Rho4,NBCH,dtmid)

! Step3: rho4 is transformed and copied to RMM
       RhoReal=REAL(Rho4)
       call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoOM  => RhoOA
       ! DIVERGENCE
       call rmm_put_rho(M,RhoReal)

! Step4: Density matrix 4 is used to calculate F5
       call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
       call rmm_get_sfm(M,Fock5)
       call MatrixTransform(M,Linv,Fock5,Uinv) ! FockOA => FockOM
       DEALLOCATE(Rho4,Fock3,RhoReal)
       call g2g_timer_stop('magnuspred')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
