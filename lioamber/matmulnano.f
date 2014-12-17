            subroutine matmulnano(A,B,C,M)
!!!!!!!!  Hace C=Bt*(A*B) para matrices cuadradas
       REAL*8 , intent(in)  :: A(M,M), B(M,M)
       REAL*8,intent(out) :: C(M,M)
       INTEGER*8,intent(in) :: M
       logical ta, tb

       real*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
       
        call DGEMM('T','N',M,M,M,1.0D0,B,M,A,M,0.0D0,scratch,M)
        call DGEMM('N','N',M,M,M,1.0D0,scratch,M,B,M,0.0D0,C,M)
       deallocate (scratch)
       return


       end












