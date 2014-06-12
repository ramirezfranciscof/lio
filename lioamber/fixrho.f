!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE fixrho(NM,Matrix)
!------------------------------------------------------------------------------!
!
! When taken out of RMM packed storage (SP), rho matrix has doubled
! elements in non diagonal positions. This subroutine fixes that.
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
           Matrix(ii,jj)=Matrix(ii,jj)/2
         ENDIF
       ENDDO
       ENDDO
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
