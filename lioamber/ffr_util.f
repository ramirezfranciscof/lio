!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  FFR_UTIL
!------------------------------------------------------------------------------!
!
! General use subroutines
!
! ffr_print_dmatrix(nunit,M,dmatrix)
! ffr_print_conmut(nunit,M,amat,bmat)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE ffr_print_dmatrix(nunit,M,dmatrix)
       implicit none
       integer,intent(in) :: nunit,M
       real*8,intent(in)  :: dmatrix(M,M)
       integer            :: ii,jj
!------------------------------------------------------------------------------!
       write(nunit,'(A)') ' Printing Matrix:'
       do ii=1,M;do jj=1,M
         write(nunit,200) ii,jj,dmatrix(ii,jj)
       enddo;enddo
       write(nunit,'(A)') '----------'
  200  FORMAT('(',I3,',',I3,')',3X,E18.10)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE ffr_print_conmut(nunit,M,amat,bmat)
       implicit none
       integer,intent(in) :: nunit,M
       real*8,intent(in)  :: amat(M,M),bmat(M,M)
       real*8,allocatable :: matrix1(:,:),matrix2(:,:),matrix(:,:)
       integer            :: ii,jj
!------------------------------------------------------------------------------!
       allocate(matrix1(M,M),matrix2(M,M),matrix(M,M))
       matrix1=matmul(amat,bmat)
       matrix2=matmul(bmat,amat)
       matrix=matrix1-matrix2
       call ffr_print_dmatrix(nunit,M,matrix)
       deallocate(matrix1,matrix2,matrix)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
