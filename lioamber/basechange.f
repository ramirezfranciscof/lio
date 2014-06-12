!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE basechange(NM,Cdir,Matrix,Cinv)
!------------------------------------------------------------------------------!
!
! This subroutine changes the expression of the input
! Matrix from base B1 to base B2. It uses the base change
! matrixes Cdir and Cinv, where Cdir is the matrix that
! goes from base B1 to base B2 (and has in its columns
! the vectors of B1 written in the base B2) and Cinv is
! its inverse (which coincides with the matrix that goes
! from base B2 to base B1).
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
