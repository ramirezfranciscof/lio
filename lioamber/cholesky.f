!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cholesky(M,Smat,Lmat,Linv,Umat,Uinv)
!------------------------------------------------------------------------------!
! This subroutine receives a simetric matrix and returns
! the matrixes of the Cholesky LU decomposition, such that:
!  Id   = Xt * Smat * X
!  Smat = Y  *  Id  * Yt = Lmat * Umat
!
! References for further information:
!   (*) JChemPhys 130 224106 (2009)
!   (*) Modern Quantum Chemistry, Attila Szabo
!
! 07/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                :: M
       REAL*8,INTENT(IN),DIMENSION(M,M)  :: Smat
       REAL*8,INTENT(OUT),DIMENSION(M,M) :: Lmat,Linv,Umat,Uinv
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: Scpy
       CHARACTER(LEN=1)                  :: UPLO
       INTEGER                           :: ErrID,ii,jj
!
!------------------------------------------------------------------------------!
!
       CALL g2g_timer_start('cholesky')
       ALLOCATE(Scpy(M,M))
       Scpy = Smat
       Lmat = 0.0d0
       Linv = 0.0d0
       Umat = 0.0d0
       Uinv = 0.0d0
       UPLO='L'

       CALL g2g_timer_start('DPOTRF')
       CALL DPOTRF(UPLO,M,Scpy,M,ErrID)
       CALL g2g_timer_stop('DPOTRF')

       CALL g2g_timer_start('depack')
       do ii=1,M
       do jj=1,ii
         Lmat(ii,jj)=Scpy(ii,jj)
         Umat(jj,ii)=Scpy(ii,jj)
       enddo
       enddo
       CALL g2g_timer_stop('depack')


       CALL g2g_timer_start('DTPTRI')
       CALL DTRTRI(UPLO,'N',M,Scpy,M,ErrID)
       CALL g2g_timer_stop('DTPTRI')

       CALL g2g_timer_start('depack')
       do ii=1,M
       do jj=1,ii
         Linv(ii,jj)=Scpy(ii,jj)
         Uinv(jj,ii)=Scpy(ii,jj)
       enddo
       enddo
       CALL g2g_timer_stop('depack')

       DEALLOCATE(Scpy)
       CALL g2g_timer_stop('cholesky')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! AUXILIARIES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cholesky_vec(M,Svec,Lmat,Linv,Umat,Uinv)
       IMPLICIT NONE
       INTEGER,INTENT(IN)                :: M
       REAL*8,INTENT(IN),DIMENSION(M,M)  :: Svec(M*(M+1)/2)
       REAL*8,INTENT(OUT),DIMENSION(M,M) :: Lmat,Linv,Umat,Uinv
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: Smat
!------------------------------------------------------------------------------!
       ALLOCATE(Smat(M,M))
       CALL spunpack('L',M,Svec,Smat)
       CALL cholesky(M,Smat,Lmat,Linv,Umat,Uinv)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
