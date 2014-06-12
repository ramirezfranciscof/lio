!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE test_intSG(nat,NM,QM)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:DSX,DSY,DSZ
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: nat,NM
       REAL*8,INTENT(IN)  :: QM(NM,NM)

       REAL*8,ALLOCATABLE :: Force(:,:),Vector(:)
       REAL*8,ALLOCATABLE :: DSM(:,:),DSV(:,:,:)
       LOGICAL            :: dototal,dotest(2)
       INTEGER            :: NV,ndir,iat,ii,jj,kk
!------------------------------------------------------------------------------!
       NV=NM*(NM+1)/2
       ALLOCATE(Force(nat,3),Vector(NV))
       ALLOCATE(DSM(NM,NM),DSV(NV,nat,3))
       dotest(:)=.TRUE.
       dototal=.FALSE.
!
!----------------------------------------------------------!
! (1) ¿ESTOY SACANDO DS?
!----------------------------------------------------------!
!
       IF (dotest(1)) THEN
         OPEN(UNIT=501,FILE='DSM.dat')

         CALL calcDSV(DSV)
         Force=0.0d0
         CALL calcDSM(Force)

         DO iat=1,nat
         DO ndir=1,3

           DO kk=1,NV
             Vector(kk)=DSV(kk,iat,ndir)
           ENDDO
           CALL spunpack('L',NM,Vector,DSM)

           WRITE(501,*) 'Atom: ',iat,'   Dir: ',ndir
           WRITE(501,*)
           DO ii=1,NM
           DO jj=1,NM
             IF (ndir.EQ.1) WRITE(501,*)
     >       DSM(ii,jj),DSX(ii,jj,iat),DSM(ii,jj)-DSX(ii,jj,iat)
             IF (ndir.EQ.2) WRITE(501,*) 
     >       DSM(ii,jj),DSY(ii,jj,iat),DSM(ii,jj)-DSY(ii,jj,iat)
             IF (ndir.EQ.3) WRITE(501,*) 
     >       DSM(ii,jj),DSZ(ii,jj,iat),DSM(ii,jj)-DSZ(ii,jj,iat)
           ENDDO
           ENDDO
           WRITE(501,*) '----------------------------------------'

         ENDDO
         ENDDO

         CLOSE(UNIT=501)
       ENDIF
!
!----------------------------------------------------------!
! (2) ¿ESTOY INTERPRETANDO BIEN EL PRODUCTO?
!----------------------------------------------------------!
!
       IF (dotest(2)) THEN
         OPEN(UNIT=501,FILE='ff0code.dat')
         OPEN(UNIT=502,FILE='ff0mano.dat')

         Force=0.0d0
         CALL intSG(Force)
         IF (dototal) CALL int1G(Force)
         IF (dototal) CALL int3G(Force,.TRUE.)
         DO iat=1,nat
           WRITE(501,*) Force(iat,1),Force(iat,2),Force(iat,3)
         ENDDO

         Force=0.0d0
         CALL calcDSM(Force)
         Force=0.0d0
         DO iat=1,nat
           DO ii=1,NM
           DO jj=1,NM
             Force(iat,1)=Force(iat,1)+QM(ii,jj)*DSX(jj,ii,iat)
             Force(iat,2)=Force(iat,2)+QM(ii,jj)*DSY(jj,ii,iat)
             Force(iat,3)=Force(iat,3)+QM(ii,jj)*DSZ(jj,ii,iat)
           ENDDO
           ENDDO
         ENDDO

         IF (dototal) CALL int1G(Force)
         IF (dototal) CALL int3G(Force,.TRUE.)
         DO iat=1,nat
           WRITE(502,*) Force(iat,1),Force(iat,2),Force(iat,3)
         ENDDO

         CLOSE(UNIT=501)
         CLOSE(UNIT=502)
       ENDIF
!
!----------------------------------------------------------!
       DEALLOCATE(Force,Vector)
       DEALLOCATE(DSM,DSV)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
