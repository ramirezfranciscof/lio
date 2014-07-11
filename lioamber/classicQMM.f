!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE classicQMM(M,QMM,NCO,EigVal,CoefMat)
!------------------------------------------------------------------------------!
!
! Subrutina de apoyo para el calculo de las contribuciones a
! la fuerza sobre los atomos que aparecen cuando el estado
! cuantico es el fundamental.
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
       INTEGER,INTENT(IN) :: M,NCO
       REAL*8,INTENT(OUT) :: QMM(M*(M+1)/2)
       REAL*8,INTENT(IN)  :: EigVal(NCO)
       REAL*8,INTENT(IN)  :: CoefMat(M,NCO)

       REAL*8             :: New,factor
       INTEGER            :: ii,jj,kk,idx
!
!------------------------------------------------------------------------------!
!
       DO jj=1,M
       DO ii=jj,M

         factor=4.0d0
         if (ii.eq.jj) factor=2.0d0
         idx=ii+(2*M-jj)*(jj-1)/2

         QMM(idx)=0.0d0
         DO kk=1,NCO
           New=factor
           New=New*EigVal(kk)
           New=New*CoefMat(ii,kk)*CoefMat(jj,kk)
           QMM(idx)=QMM(idx)-New
         ENDDO

       ENDDO
       ENDDO

       RETURN;END SUBROUTINE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
