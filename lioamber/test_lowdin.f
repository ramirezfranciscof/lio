!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE test_lowdin(nat,NM,Fock,UInv,Rho,Eigvec,Eigval)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:DSX,DSY,DSZ,RMM
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: nat,NM
       REAL*8,INTENT(IN)  :: Fock(NM,NM),UInv(NM,NM),Rho(NM,NM),
     >                       Eigvec(NM,NM),Eigval(NM)

       REAL*8,ALLOCATABLE :: Force(:,:),Vector(:)
       REAL*8,ALLOCATABLE :: DSM(:,:),DSV(:,:,:)
       REAL*8             :: scratch
       LOGICAL            :: dototal,dotest(1)
       INTEGER            :: NV,ndir,iat,ii,jj,kk
!------------------------------------------------------------------------------!
       NV=NM*(NM+1)/2
       ALLOCATE(Force(nat,3),Vector(NV))
       ALLOCATE(DSM(NM,NM),DSV(NV,nat,3))
       dototal=.FALSE.
       dotest(:)=.TRUE.
!
!----------------------------------------------------------!
! FFR-Q: LODWIN PASO A PASO
!
       IF (dotest(1)) THEN
         CALL g2g_timer_start('lowdin-test1')
         OPEN(UNIT=501,FILE='fflowdin.dat')
         CALL calcDSM(Force)
         Force=0.0d0

         DO iat=1,nat
         DO ndir=1,3

           DO ii=1,NM
           DO jj=1,NM
             IF (ndir.EQ.1) DSM(ii,jj)=DSX(ii,jj,iat)
             IF (ndir.EQ.2) DSM(ii,jj)=DSY(ii,jj,iat)
             IF (ndir.EQ.3) DSM(ii,jj)=DSZ(ii,jj,iat)
           ENDDO
           ENDDO

           CALL lowdin_force(scratch,NM,Fock,DSM,UInv,Rho,Eigvec,Eigval)
           Force(iat,ndir)=scratch

         ENDDO
         ENDDO

         IF (dototal) CALL int1G(Force)
         IF (dototal) CALL int3G(Force,.TRUE.)
         DO iat=1,nat
           WRITE(501,*) Force(iat,1),Force(iat,2),Force(iat,3)
         ENDDO

         CLOSE(UNIT=501)
         CALL g2g_timer_stop('lowdin-test1')
       ENDIF
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
