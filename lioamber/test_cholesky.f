!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE test_cholesky(nat,NM,rmmptr,Fock,LInv,UInv,Rho)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:DSX,DSY,DSZ,RMM
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: nat,NM,rmmptr
       REAL*8,INTENT(IN)  :: Fock(NM,NM),LInv(NM,NM),UInv(NM,NM),
     >                       Rho(NM,NM)

       REAL*8,ALLOCATABLE :: Force(:,:),Vector(:)
       REAL*8,ALLOCATABLE :: DSM(:,:),DSV(:,:,:)
       REAL*8             :: scratch
       LOGICAL            :: dototal,dolast,dotest(1)
       INTEGER            :: NV,ndir,iat,ii,jj,kk
!------------------------------------------------------------------------------!
       NV=NM*(NM+1)/2
       ALLOCATE(Force(nat,3),Vector(NV))
       ALLOCATE(DSM(NM,NM),DSV(NV,nat,3))
       dototal=.FALSE.
       dolast=.TRUE.
       dotest(:)=.TRUE.
!
!----------------------------------------------------------!
! FFR-Q: PASO A PASO
!
!
       IF (dotest(1)) THEN
         CALL g2g_timer_start('cholesky-test1')
         OPEN(UNIT=501,FILE='ffcholmano.dat')

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

           CALL cholesky_force(scratch,NM,Fock,LInv,DSM,UInv,Rho)
           Force(iat,ndir)=scratch

         ENDDO
         ENDDO

         IF (dototal) CALL int1G(Force)
         IF (dototal) CALL int3G(Force,.TRUE.)
         DO iat=1,nat
           WRITE(501,*) Force(iat,1),Force(iat,2),Force(iat,3)
         ENDDO
         CLOSE(UNIT=501)
         CALL g2g_timer_stop('cholesky-test1')
       ENDIF
!
!----------------------------------------------------------!
! FFR-Q: LA OTRA PROPAGACION VIA CODIGO VIEJO
!
       IF (dolast) THEN
         CALL g2g_timer_start('cholesky-test2')
         OPEN(UNIT=501,FILE='ffcholcode.dat')

         DO kk=1,NV
           Vector(kk)=RMM(rmmptr+kk-1)
         ENDDO

         CALL prepQMM(NM,RMM(rmmptr),UInv,Rho,Fock,LInv)
         Force=0.0d0
         CALL intSG(Force)
         IF (dototal) CALL int1G(Force)
         IF (dototal) CALL int3G(Force,.TRUE.)
         DO iat=1,nat
           WRITE(501,*) Force(iat,1),Force(iat,2),Force(iat,3)
         ENDDO

         DO kk=1,NV
           RMM(rmmptr+kk-1)=Vector(kk)
         ENDDO

         CLOSE(UNIT=501)
         CALL g2g_timer_stop('cholesky-test2')
       ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
