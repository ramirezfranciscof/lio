            subroutine rho_transform(A,X,C,M)
!!!!!!!!  Hace C=Xt*(A*X)
#ifdef CUBLAS
       integer*8 , intent(in)  :: X
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
       call cumxp(A,X,scratch,M)
       call cumpxt(scratch,X,C,M)
       deallocate (scratch)
#else
       REAL*8 , intent(in)  :: X(M,M)
       integer, intent(in) :: M
       logical ta, tb
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: A(M,M)
       COMPLEX*8, intent(out) :: C(M,M)
       COMPLEX*8 :: alpha,beta
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#else
       COMPLEX*16, intent(in) :: A(M,M)
       COMPLEX*16, intent(out) :: C(M,M)
       COMPLEX*16 :: alpha,beta
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#endif
       allocate (scratch(M,M),scratch1(M,M))
       alpha=cmplx(1.0D0,0.0D0)
       beta=cmplx(0.0D0,0.0D0)
       DO i=1,M
          DO j=1,M
             scratch1(i,j)=cmplx(X(i,j),0.0D0)
          ENDDO
       ENDDO
#ifdef TD_SIMPLE
        call CGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call CGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#else       
        call ZGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call ZGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#endif
      deallocate(scratch,scratch1)
#endif
       return
       end
       subroutine complex_rho_on_to_ao(A,X,C,M)
!!!!!!!!  Hace C=Xt*(A*X)
#ifdef CUBLAS
       integer*8 , intent(in)  :: X
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
       call cumxp(A,X,scratch,M)
       call cumpxt(scratch,X,C,M)
       deallocate (scratch)
#else
       REAL*8 , intent(in)  :: X(M,M)
       integer, intent(in) :: M
       logical ta, tb
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: A(M,M)
       COMPLEX*8, intent(out) :: C(M,M)
       COMPLEX*8 :: alpha,beta
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#else
       COMPLEX*16, intent(in) :: A(M,M)
       COMPLEX*16, intent(out) :: C(M,M)
       COMPLEX*16 :: alpha,beta
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#endif
       allocate (scratch(M,M),scratch1(M,M))
       alpha=cmplx(1.0D0,0.0D0)
       beta=cmplx(0.0D0,0.0D0)
       DO i=1,M
          DO j=1,M
             scratch1(i,j)=cmplx(X(i,j),0.0D0)
          ENDDO
       ENDDO
#ifdef TD_SIMPLE
        call CGEMM('N','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call CGEMM('N','T',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#else       
        call ZGEMM('N','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call ZGEMM('N','T',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#endif
      deallocate(scratch,scratch1)
#endif
       return
       end

       subroutine complex_rho_ao_to_on(A,X,C,M)
!!!!!!!!  Hace C=Xt*(A*X)
#ifdef CUBLAS
       integer*8 , intent(in)  :: X
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
       call cumxtp(A,X,scratch,M)
       call cumpx(scratch,X,C,M)
       deallocate (scratch)
#else
       REAL*8 , intent(in)  :: X(M,M)
       integer, intent(in) :: M
       logical ta, tb
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: A(M,M)
       COMPLEX*8, intent(out) :: C(M,M)
       COMPLEX*8 :: alpha,beta
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#else
       COMPLEX*16, intent(in) :: A(M,M)
       COMPLEX*16, intent(out) :: C(M,M)
       COMPLEX*16 :: alpha,beta
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch,scratch1
#endif
       allocate (scratch(M,M),scratch1(M,M))
       alpha=cmplx(1.0D0,0.0D0)
       beta=cmplx(0.0D0,0.0D0)
       DO i=1,M
          DO j=1,M
             scratch1(i,j)=cmplx(X(i,j),0.0D0)
          ENDDO
       ENDDO
#ifdef TD_SIMPLE
        call CGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call CGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#else       
        call ZGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call ZGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#endif
      deallocate(scratch,scratch1)
#endif
       return
       end




















