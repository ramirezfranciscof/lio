!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE messrho(NM,Matrix)
!------------------------------------------------------------------------------!
!
! When taken out of RMM packed storage (SP), rho matrix has doubled
! elements in non diagonal positions. So before packing rho, one has
! to mess it up.
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Matrix(NM,NM)
       INTEGER            :: ii,jj
!
!------------------------------------------------------------------------------!
!
       DO ii=1,NM
       DO jj=1,NM
         IF (ii.ne.jj) THEN
           Matrix(ii,jj)=2*Matrix(ii,jj)
         ENDIF
       ENDDO
       ENDDO
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
