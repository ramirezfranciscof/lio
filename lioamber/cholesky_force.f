!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cholesky_force(ff,M,Fock,LInv,DSM,UInv,Rho)
!------------------------------------------------------------------------------!
!
! Subrutina de apoyo para el calculo de las contribuciones a
! la fuerza sobre los atomos que aparecen cuando el estado
! cuantico no es el fundamental.
!
! Las matrices de entrada deben estar en la base de orbitales
! atomicos (no ortonormal).
!
!
! Fuentes:
! (X) JChemPhys 130 224106 (2009)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       REAL*8,INTENT(OUT)               :: ff
       INTEGER,INTENT(IN)               :: M
       REAL*8,INTENT(IN),DIMENSION(M,M) :: Fock,LInv,DSM,UInv,Rho


       REAL*8,ALLOCATABLE :: USU(:,:),DUU(:,:),TRA(:,:),PFM(:,:)
       REAL*8             :: suma
       INTEGER            :: aa,bb,ii,jj,kk
       INTEGER            :: TestID
!
!------------------------------------------------------------------------------!
!
       TestID=1
       ff=0.0d0
       ALLOCATE(USU(M,M),DUU(M,M),TRA(M,M),PFM(M,M))
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.1) THEN
         PFM=MATMUL(DSM,UInv)
         USU=MATMUL(LInv,PFM)

         DO ii=1,M
         DO jj=1,M
           IF (ii.LT.jj) DUU(ii,jj)=USU(ii,jj)
           IF (ii.EQ.jj) DUU(ii,jj)=USU(ii,jj)/2
           IF (ii.GT.jj) DUU(ii,jj)=0.0d0
         ENDDO
         ENDDO

         PFM=MATMUL(DUU,Rho)
         TRA=MATMUL(Fock,PFM)

         DO kk=1,M
           ff=ff+2*TRA(kk,kk)
         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.2) THEN
         PFM=MATMUL(DSM,UInv)
         USU=MATMUL(LInv,PFM)

         DO ii=1,M
         DO jj=1,M
           IF (ii.LT.jj) DUU(ii,jj)=USU(ii,jj)
           IF (ii.EQ.jj) DUU(ii,jj)=USU(ii,jj)/2
           IF (ii.GT.jj) DUU(ii,jj)=0.0d0
         ENDDO
         ENDDO

         PFM=MATMUL(Rho,Fock)

         DO ii=1,M
         DO jj=1,ii
           ff=ff+2*PFM(ii,jj)*DUU(jj,ii)
         ENDDO
         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.3) THEN
         PFM=MATMUL(DSM,UInv)
         USU=MATMUL(LInv,PFM)
         PFM=MATMUL(Rho,Fock)

         DO ii=1,M
           DO jj=1,ii-1
             ff=ff+2*PFM(ii,jj)*USU(jj,ii)
           ENDDO
           ff=ff+PFM(ii,ii)*USU(ii,ii)
         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.4) THEN
         PFM=MATMUL(Rho,Fock)

         DO ii=1,M

           suma=0.0d0
           DO aa=1,M
           DO bb=1,M
             suma=suma+Linv(ii,aa)*DSM(aa,bb)*UInv(bb,ii)
           ENDDO
           ENDDO
           ff=ff+suma*PFM(ii,ii)

           DO jj=1,ii-1
             suma=0.0d0
             DO aa=1,M
             DO bb=1,M
               suma=suma+Linv(jj,aa)*DSM(aa,bb)*UInv(bb,ii)
             ENDDO
             ENDDO
             ff=ff+2*suma*PFM(ii,jj)
           ENDDO

         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.5) THEN
         PFM=MATMUL(Rho,Fock)

         DO aa=1,M
         DO bb=1,M

           suma=0.0d0

           DO ii=1,M
             suma=suma+PFM(ii,ii)*LInv(ii,aa)*UInv(bb,ii)
             DO jj=1,ii-1
               suma=suma+2*PFM(ii,jj)*LInv(jj,aa)*UInv(bb,ii)
             ENDDO
           ENDDO

           ff=ff+suma*DSM(aa,bb)

         ENDDO
         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       IF (TestID.EQ.6) THEN
         PFM=MATMUL(Rho,Fock)

         DO aa=1,M
         DO bb=1,M
           suma=0.0d0

           DO ii=1,M
             suma=suma+PFM(ii,ii)*LInv(ii,aa)*UInv(bb,ii)
             DO jj=1,ii-1
               suma=suma+2*PFM(ii,jj)*LInv(jj,aa)*UInv(bb,ii)
             ENDDO
           ENDDO

           TRA(bb,aa)=suma
         ENDDO
         ENDDO


         DO ii=1,M
         DO jj=1,M
           ff=ff+TRA(ii,jj)*DSM(jj,ii)
         ENDDO
         ENDDO
       ENDIF
!
!----------------------------------------------------------!
!
       DEALLOCATE(USU,DUU,TRA,PFM)
20     FORMAT(4(E15.8,2x))
       RETURN;END SUBROUTINE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
