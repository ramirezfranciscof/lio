!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE prepQMM(M,QMM,UInv,Rho,Fock,LInv)
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
       REAL*8,INTENT(IN)        :: Uinv(M,M),Linv(M,M)
       REAL*8,INTENT(IN)        :: Rho(M,M),Fock(M,M)
       REAL*8,INTENT(INOUT)     :: QMM(M*(M+1)/2)

       REAL*8,ALLOCATABLE       :: Qaux(:,:)
       REAL*8                   :: scratch
       INTEGER                  :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       ALLOCATE(Qaux(M,M))
       CALL prepQM(M,Qaux,UInv,Rho,Fock,LInv)

       DO ii=1,M
       DO jj=1,ii
           idx=ii+(2*M-jj)*(jj-1)/2
           QMM(idx)=Qaux(ii,jj)+Qaux(jj,ii)
           IF (ii.EQ.jj) QMM(idx)=QMM(idx)/2
       ENDDO
       ENDDO

       DEALLOCATE(Qaux)
       RETURN;END SUBROUTINE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
