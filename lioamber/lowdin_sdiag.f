!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE lowdin_sdiag(NM,SMvec,X,Xt,Y,Yt,eigvec_out,eigval_out)
!------------------------------------------------------------------------------!
!
!
! For more information, read:
!   (*) Modern Quantum Chemistry, Attila Szabo
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                  :: NM
       REAL*8,INTENT(IN)                   :: SMvec(NM*(NM+1)/2)
       REAL*8,INTENT(OUT),DIMENSION(NM,NM) :: X,Xt,Y,Yt
       REAL*8,INTENT(OUT),OPTIONAL         :: eigvec_out(NM,NM),
     >                                        eigval_out(NM)

       REAL*8,ALLOCATABLE :: eigvec(:,:),eigval(:)
       REAL*8,ALLOCATABLE :: SMcpy(:),work(:),MAT(:,:)
       REAL*8             :: tolerance,num
       INTEGER            :: NV,ErrID,ii,jj
!
!------------------------------------------------------------------------------!
!
       NV=NM*(NM+1)/2
       ALLOCATE(eigvec(NM,NM),eigval(NM))
       ALLOCATE(SMcpy(NV),work(NV),MAT(NM,NM))

       SMcpy=SMvec
       CALL DSPEV('V','L',NM,SMcpy,eigval,eigvec,NM,work,ErrID)
       IF (PRESENT(eigval_out)) eigval_out=eigval
       IF (PRESENT(eigvec_out)) eigvec_out=eigvec

!-----(LINEAR DEPENDENCY ELIMINATION)
!       tolerance=1.0D-06
!       DO jj=1,NM
!         IF (eigval(jj).LT.tolerance) THEN
!           WRITE(*,*) 'LINEAR DEPENDENCY DETECTED'
!           DO ii=1,NM
!             X(ii,jj)=0.0D0
!             Y(ii,jj)=0.0D0
!           ENDDO
!         ENDIF
!       ENDDO
!------

       DO jj=1,NM
         num=SQRT(eigval(jj))
         DO ii=1,NM
           ! THE FOLLOWING ARE BEING USED AS INTERMEDIATES
           Xt(ii,jj)  = eigvec(ii,jj)/num
           Yt(ii,jj)  = eigvec(ii,jj)*num
           MAT(ii,jj) = eigvec(jj,ii)
         ENDDO
       ENDDO
       X  = MATMUL(Xt,MAT)
       Y  = MATMUL(Yt,MAT)
       Xt = X
       Yt = Y

       DEALLOCATE(SMcpy,work,MAT)
       DEALLOCATE(eigvec,eigval)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
