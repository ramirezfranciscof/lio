            subroutine conmut(A,B,C,M)
!!!!!!!!  Hace C=A*B-B*A
       REAL*8 , intent(in)  :: A(M,M), B(M,M)
       REAL*8 , intent(out)  :: C(M,M)
       REAL*8 :: alpha, beta
!---------------------------------------------------------------------!
       alpha=1.0D0
       beta=0.0D0
       call DGEMM('N','N',M,M,M,alpha,B,M,A,M,beta,C,M)
       beta=-1.0D0
       call DGEMM('N','N',M,M,M,alpha,A,M,B,M,beta,C,M)
       return
       end












