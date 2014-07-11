!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE prepQM(M,QM,UInv,Rho,Fock,LInv)
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
       INTEGER,INTENT(IN)       :: M
       REAL*8,INTENT(IN)        :: UInv(M,M),LInv(M,M)
       REAL*8,INTENT(IN)        :: Rho(M,M),Fock(M,M)
       REAL*8,INTENT(INOUT)     :: QM(M,M)

       REAL*8,ALLOCATABLE       :: PFM(:,:)
       REAL*8                   :: suma
       INTEGER                  :: ii,jj,aa,bb
!
!------------------------------------------------------------------------------!
!
       ALLOCATE(PFM(M,M))
       PFM=MATMUL(Rho,Fock)


       DO aa=1,M
       DO bb=1,M

         suma=0.0d0
         DO ii=1,M
           suma=suma+UInv(aa,ii)*PFM(ii,ii)*LInv(ii,bb)
           IF (ii.GT.1) THEN
           DO jj=1,ii-1
             suma=suma+2*UInv(aa,ii)*PFM(ii,jj)*LInv(jj,bb)
           ENDDO
           ENDIF
         ENDDO
         QM(aa,bb)=-suma

       ENDDO
       ENDDO


       DEALLOCATE(PFM)
       RETURN;END SUBROUTINE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
