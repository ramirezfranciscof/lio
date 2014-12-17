            subroutine matmulnanoc(A,B,C,M)
!!!!!!!!  Hace C=Bt*(A*B)
       REAL*8 , intent(in)  :: B(M,M)
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
             scratch1(i,j)=cmplx(B(i,j),0.0D0)
          ENDDO
       ENDDO
#ifdef TD_SIMPLE
        call CGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call CGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#else       
        call ZGEMM('T','N',M,M,M,alpha,scratch1,M,A,M,beta,scratch,M)
        call ZGEMM('N','N',M,M,M,alpha,scratch,M,scratch1,M,beta,C,M)
#endif
!
!        do i=1,M
!        do j=1,M
!         scratch(i,j)=A(j,i)
!        enddo
!        enddo
!         scratch1=0
!        do i=1,M
!        do j=1,M
!         do k= 1,M
!         scratch1(i,j)= scratch1(i,j) + scratch(k,i)*B(k,j)
!         enddo
!        enddo
!        enddo
!         C=0
!        do i=1,M
!        do j=1,M
!         do k= 1,M
!        C(i,j)= C(i,j) + B(k,i)*scratch1(k,j)
!         enddo
!       enddo
!       enddo
!       DO i=1,M
!          DO j=1,M
!             write(10001,*) C(i,j)
!          ENDDO
!       ENDDO
       return
       end












