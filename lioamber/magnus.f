!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE magnus(M,FockMid,RhoOld,RhoNew,N,dt)
!------------------------------------------------------------------------------!
! This subroutine applies the magnus propagator using the
! BCH formula of order N to propagate the density matrix
! RhoOld into RhoNew, using FockMid. It assumes that both
! input matrixes are transformed (aka: "in the orthonormal
! basis") and returns RhoNew in the same way.
!
! RhoOld  => t0         [in]
! FockMid => t0+(dt/2)  [in]
! RhoNew  => t0+dt      [out]
!
! References for further information:
!   (*) J. Chem. Theory Comput. 2011, 7, 1344–1355
!
! 07/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                    :: M,N
       REAL*8,INTENT(IN)                     :: FockMid(M,M),dt
       COMPLEX*8,INTENT(IN)                  :: RhoOld(M,M)
       COMPLEX*8,INTENT(OUT)                 :: RhoNew(M,M)
!
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:)  :: Conmut,TermPos,TermNeg
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:)  :: Omega1
       INTEGER                               :: kk
       REAL*8                                :: invfact
!
       COMPLEX*8,PARAMETER :: icmplx=CMPLX(0.0D0,1.0D0)
!
!------------------------------------------------------------------------------!
!
       call g2g_timer_start('magnus')
       allocate(Conmut(M,M),TermPos(M,M),TermNeg(M,M))
       allocate(Omega1(M,M))
!
       Omega1=(-1)*(icmplx)*(FockMid)*(dt)
       RhoNew=RhoOld
       Conmut=RhoOld
!
       invfact=1
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
       call g2g_timer_stop('magnus')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
