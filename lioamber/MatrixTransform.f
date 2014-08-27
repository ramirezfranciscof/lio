!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE MatrixTransform(NM,Cdir,Matrix,Cinv)
!------------------------------------------------------------------------------!
!
! This subroutine transforms the input Matrix using Cdir
! and Cinv. The usual transformations this will be applied
! to are the following:
!
! Fock(OM) = Linv * Fock(OA) * Uinv = Xt * Fock(OA) * X
! Fock(OA) = Lmat * Fock(OM) * Umat = Y  * Fock(OM) * Yt
! Rho(OM)  = Umat * Rho(OA)  * Lmat = Yt * Rho(OA)  * Y
! Rho(OA)  = Uinv * Rho(OM)  * Linv = X  * Rho(OM)  * Xt
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       INTEGER,INTENT(IN)   :: NM
       REAL*8,INTENT(IN)    :: Cdir(NM,NM),Cinv(NM,NM)
       REAL*8,INTENT(INOUT) :: Matrix(NM,NM)
       REAL*8,ALLOCATABLE   :: Output(:,:)
       REAL*8               :: Termino
       INTEGER              :: aa,bb,ii,jj
!
!------------------------------------------------------------------------------!
!
       ALLOCATE(Output(NM,NM))

       DO aa=1,NM
       DO bb=1,NM

         Output(aa,bb)=0.0d0
         DO ii=1,NM
         DO jj=1,NM
           Termino=Cdir(aa,ii)*Matrix(ii,jj)*Cinv(jj,bb)
           Output(aa,bb)=Output(aa,bb)+Termino
         ENDDO
         ENDDO

       ENDDO
       ENDDO

       Matrix=Output
       DEALLOCATE(Output)
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
