!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! PACKED STORAGE
!------------------------------------------------------------------------------!
!
! This two subroutines are intended to facilitate the passage
! from a matrix and a vector with the same information (the
! vector has it in what is known as "packed storage" or SP).
! The matrix must be symetric for packed storage to be used,
! and this routines only support real values.
!
! NOTE: If Vector is originally of size > NM*(NM+1)/2, then
! calling the subroutine with "Vector(i0)" will make it
! work with the positions of the original vector that go
! from i0 to i0+NM*(NM+1)/2-1.
! (something similar may work with the variable 'Matrix'
! but it will probably be more intrincated)
!
!-----------------------!
! STORAGE SCHEME UPLO=U
!-----------------------!
!  | 1  7  8 |
!  |(7) 2  9 |   <=>   ( 1 , 7 , 2 , 8 , 9 , 3 )
!  |(8)(9) 3 |
!
!-----------------------!
! STORAGE SCHEME UPLO=L
!-----------------------!
!  | 1 (7)(8)|
!  | 7  2 (9)|   <=>   ( 1 , 7 , 8 , 2 , 9 , 3 )
!  | 8  9  3 |
!
!----------------------------------------------------------!
! For more information, visit:
!   (*) http://www.netlib.org/lapack/lug/node123.html
!   (*) http://www.netlib.org/lapack/lug/node24.html
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE spunpack(UPLO,NM,Vector,Matrix)
       IMPLICIT NONE
       CHARACTER(LEN=1)   :: UPLO
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(IN)  :: Vector(NM*(NM+1)/2)
       REAL*8,INTENT(OUT) :: Matrix(NM,NM)
       INTEGER            :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       IF (UPLO.EQ.'U') THEN
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           Matrix(ii,jj)=Vector(idx)
           Matrix(jj,ii)=Vector(idx)
         ENDDO
         ENDDO

       ELSE IF (UPLO.EQ.'L') THEN
         DO ii=1,NM
         DO jj=1,ii
           idx=ii+(2*NM-jj)*(jj-1)/2
           Matrix(ii,jj)=Vector(idx)
           Matrix(jj,ii)=Vector(idx)
         ENDDO
         ENDDO

       ELSE
         PRINT*,'NOT GOOD INPUT FOR UPLO'

       ENDIF
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE sprepack(UPLO,NM,Vector,Matrix)
       IMPLICIT NONE
       CHARACTER(LEN=1)   :: UPLO
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Vector(NM*(NM+1)/2)
       REAL*8,INTENT(IN)  :: Matrix(NM,NM)
       INTEGER            :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       IF (UPLO.EQ.'U') THEN
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           Vector(idx)=Matrix(ii,jj)
         ENDDO
         ENDDO

       ELSE IF (UPLO.EQ.'L') THEN
         DO ii=1,NM
         DO jj=1,ii
           idx=ii+(2*NM-jj)*(jj-1)/2
           Vector(idx)=Matrix(ii,jj)
         ENDDO
         ENDDO

       ELSE
         PRINT*,'NOT GOOD INPUT FOR UPLO'

       ENDIF
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE spunpack_rtc(UPLO,NM,Vector,Matrix)
       IMPLICIT NONE
       CHARACTER(LEN=1)   :: UPLO
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(IN)  :: Vector(NM*(NM+1)/2)
#ifdef TD_SIMPLE
       COMPLEX*8,INTENT(OUT) :: Matrix(NM,NM)
#else
       COMPLEX*16,INTENT(OUT) :: Matrix(NM,NM)
#endif
       INTEGER            :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       IF (UPLO.EQ.'U') THEN
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           IF(jj.eq.ii) then
              Matrix(jj,jj)=cmplx(Vector(idx),0.0D0)
           ELSE
              Matrix(ii,jj)=cmplx(Vector(idx),0.0D0)
              Matrix(ii,jj)=Matrix(ii,jj)*0.50D0      ! RHO CASE
              Matrix(jj,ii)=cmplx(Vector(idx),0.0D0)
              Matrix(jj,ii)=Matrix(jj,ii)*0.50D0      ! RHO CASE
           ENDIF
         ENDDO
         ENDDO

       ELSE IF (UPLO.EQ.'L') THEN
         DO ii=1,NM
         DO jj=1,ii
           idx=ii+(2*NM-jj)*(jj-1)/2
           IF(jj.eq.ii) THEN
             Matrix(ii,ii)=cmplx(Vector(idx),0.0D0)
           ELSE
             Matrix(ii,jj)=cmplx(Vector(idx),0.0D0)
             Matrix(ii,jj)=Matrix(ii,jj)*0.50D0
             Matrix(jj,ii)=cmplx(Vector(idx),0.0D0)
             Matrix(jj,ii)=Matrix(jj,ii)*0.50D0
           ENDIF
         ENDDO
         ENDDO

       ELSE
         PRINT*,'NOT GOOD INPUT FOR UPLO'

       ENDIF
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE sprepack_ctr(UPLO,NM,Vector,Matrix)
       IMPLICIT NONE
       CHARACTER(LEN=1)   :: UPLO
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Vector(NM*(NM+1)/2)
#ifdef TD_SIMPLE
       COMPLEX*8,INTENT(IN) :: Matrix(NM,NM)
#else
       COMPLEX*16,INTENT(IN) :: Matrix(NM,NM)
#endif
       INTEGER            :: ii,jj,idx
!
!------------------------------------------------------------------------------!
!
       IF (UPLO.EQ.'U') THEN
         DO jj=1,NM
         DO ii=1,jj
           idx=ii+(jj*(jj-1)/2)
           if(ii.eq.jj) then
              Vector(idx)=REAL(Matrix(ii,jj))
           else
              Vector(idx)=REAL(Matrix(ii,jj))
              Vector(idx)=Vector(idx)*2.0D0
           endif
         ENDDO           
         ENDDO

       ELSE IF (UPLO.EQ.'L') THEN
         DO ii=1,NM
         DO jj=1,ii
           idx=ii+(2*NM-jj)*(jj-1)/2
           IF(ii.eq.jj) THEN
             Vector(idx)=REAL(Matrix(ii,jj))
           ELSE
             Vector(idx)=REAL(Matrix(ii,jj))
             Vector(idx)=Vector(idx)*2.0D0
           ENDIF
         ENDDO
         ENDDO

       ELSE
         PRINT*,'NOT GOOD INPUT FOR UPLO'

       ENDIF
!
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

