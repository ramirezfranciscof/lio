            subroutine conmutc(A,B,C,M)
!========================================================================!
!!!!!!!!  Hace C=A*B-B*A
!------------------------------------------------------------------------!
       REAL*8 , intent(in)  :: A(M,M)
       INTEGER,intent(in) :: M
#ifdef TD_SIMPLE
       COMPLEX*8 :: alpha, beta
       COMPLEX*8 , intent(in) :: B(M,M)
       COMPLEX*8 , intent(out) :: C(M,M)
       COMPLEX*8, ALLOCATABLE :: scratch(:,:)
#else
       COMPLEX*16 ::  alpha, beta
       COMPLEX*16 , intent(in) :: B(M,M)
       COMPLEX*16 , intent(out) :: C(M,M)
       COMPLEX*16, ALLOCATABLE :: scratch(:,:)
#endif
!------------------------------------------------------------------------!
       ALLOCATE(scratch(M,M))
       DO i=1,M
          DO j=1,M
             scratch(i,j)=CMPLX(A(i,j),0.0D0)
          ENDDO
       ENDDO
       alpha=CMPLX(1.0D0,0.0D0)
       beta=CMPLX(0.0D0,0.0D0)
#ifdef TD_SIMPLE
       call CGEMM('N','N',M,M,M,alpha,B,M,scratch,M,beta,C,M)
       beta=CMPLX(-1.0D0,0.0D0)
       call CGEMM('N','N',M,M,M,alpha,scratch,M,B,M,beta,C,M)
#else       
       call ZGEMM('N','N',M,M,M,alpha,B,M,scratch,M,beta,C,M)
       beta=CMPLX(-1.0D0,0.0D0)
       call ZGEMM('N','N',M,M,M,alpha,scratch,M,B,M,beta,C,M)
#endif
!----------Testing-------------------------!
!        scratch=matmul(B,A)
!        scratch1=matmul(A,B)
!        C=scratch1-scratch
!----------End-Testing---------------------!
       DEALLOCATE(scratch)
       return
       end




















