            subroutine cumatmulnanoc(A,devPtrX,C,M)
!!!!!!!!  Hace C=Bt*(A*B)
       integer , intent(in)  :: devPtrX
       integer, intent(in) :: M
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: A(M,M)
       COMPLEX*8, intent(out) :: C(M,M)
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
#else
       COMPLEX*16, intent(in) :: A(M,M)
       COMPLEX*16, intent(out) :: C(M,M)
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
#endif
       allocate (scratch(M,M))
       call cumxp(A,devPtrX,scratch,M)
       call cumpxt(scratch,devPtrX,C,M)
       deallocate (scratch)
       return
       end












