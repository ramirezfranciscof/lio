            subroutine conmutcc(A,B,C,M)
!======================================================================!
!!!!!!!!  Hace C=A*B-B*A
!----------------------------------------------------------------------!
       INTEGER,intent(in) :: M
#ifdef TD_SIMPLE
       COMPLEX*8 :: alpha, beta
       COMPLEX*8 , intent(inout) :: B(M,M), C(M,M), A(M,M)
#else
       COMPLEX*16 ::  alpha, beta
       COMPLEX*16 , intent(inout) :: B(M,M), C(M,M), A(M,M)
#endif
!---------------------------------------------------------------------!
        alpha=CMPLX(1.0D0,0.0D0)
        beta=CMPLX(0.0D0,0.0D0)
#ifdef TD_SIMPLE
        call CGEMM('N','N',M,M,M,alpha,B,M,A,M,beta,C,M)
        beta=CMPLX(-1.0D0,0.0D0)
        call CGEMM('N','N',M,M,M,alpha,A,M,B,M,beta,C,M)
#else       
        call ZGEMM('N','N',M,M,M,alpha,B,M,A,M,beta,C,M)
        beta=CMPLX(-1.0D0,0.0D0)
        call ZGEMM('N','N',M,M,M,alpha,A,M,B,M,beta,C,M)
#endif
       return
       end












