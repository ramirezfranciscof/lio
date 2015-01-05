             subroutine fock_ao_to_on(A,X,C,M)
!!!!!!!!  Hace C=Bt*(A*B)
       integer, intent(in) :: M
#ifdef CUBLAS
       integer , intent(in)  :: X
       REAL*8, intent(in) :: A(M,M)
       REAL*8, intent(out) :: C(M,M)
       REAL*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
       call cumxtf(A,X,scratch,M)
       call cumfx(scratch,X,C,M)
#else
       REAL*8 , intent(in)  :: A(M,M), X(M,M)
       REAL*8,intent(out) :: C(M,M)
       real*8, dimension (:,:), ALLOCATABLE :: scratch
       allocate (scratch(M,M))
        call DGEMM('T','N',M,M,M,1.0D0,X,M,A,M,0.0D0,scratch,M)
        call DGEMM('N','N',M,M,M,1.0D0,scratch,M,X,M,0.0D0,C,M)
#endif
       deallocate (scratch)
       return
       end














