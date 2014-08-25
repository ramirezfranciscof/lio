!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE magnus_bch(M,Fock,RhoOld,RhoNew,N,dt)
!------------------------------------------------------------------------------!
!
! Magnus propagator (N order) using BCH formula
! Entrada: Fock(t+(deltat/2)), rho(t)
! Salida:  rho6=rho(t+deltat)
!
!
! Reference:  J. Chem. Theory Comput. 2011, 7, 1344–1355
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                    :: M,N
       REAL*8,INTENT(IN)                     :: Fock(M,M),dt
       COMPLEX*16,INTENT(IN)                 :: RhoOld(M,M)
       COMPLEX*16,INTENT(OUT)                :: RhoNew(M,M)

       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: Conmut,TermPos,TermNeg
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: Omega1
       INTEGER                               :: kk
       REAL*8                                :: invfact

       COMPLEX*16,PARAMETER :: icmplx=CMPLX(0.0D0,1.0D0)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       call g2g_timer_start('magnus_bch')
       allocate(Conmut(M,M),TermPos(M,M),TermNeg(M,M))
       allocate(Omega1(M,M))
!
       Omega1=(-1)*(icmplx)*(Fock)*(dt)
       RhoNew=RhoOld
       Conmut=RhoOld
       invfact=1
!
       do kk=1,N
         TermPos=MATMUL(Omega1,Conmut)
         TermNeg=MATMUL(Conmut,Omega1)
         Conmut=TermPos-TermNeg
         invfact=invfact/kk
         RhoNew=RhoNew+invfact*Conmut
       enddo
!
       deallocate(Omega1)
       deallocate(Conmut,TermPos,TermNeg)
       call g2g_timer_stop('magnus_bch')
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
