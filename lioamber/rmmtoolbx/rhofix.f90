!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  RhoFix // RhoMess
!------------------------------------------------------------------------------!
!
! In the RMM vector, Rho is stored in packed storage but
! each non diagonal position has dobule the real value.
! So before putting a density matrix in RMM, non diagonal
! positions need to be x2 (messrho) and when taking it
! out of RMM, the non diagonal positions need to be
! divided by two (fixrho).
!
! Note that fixrho modifies the matrix, whereas mess rho
! modifies the vector.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  SUBROUTINE rhofix_mat(Matrix)
!------------------------------------------------------------------------------!
  implicit none
  real*8,intent(inout) :: Matrix(M,M)
  integer              :: ii,jj
!
  do ii=1,size(Matrix,1)
  do jj=1,size(Matrix,2)
    if (ii.ne.jj) Matrix(ii,jj)=Matrix(ii,jj)/(2.0d0)
  enddo
  enddo
!
  return;end subroutine
!
!
!------------------------------------------------------------------------------!
  SUBROUTINE rhofix_vec(M,Vector)
!------------------------------------------------------------------------------!
  implicit none
  integer,intent(in)   :: M
  real*8,intent(inout) :: Vector((M*(M+1))/2)
  integer              :: nn,ii,jj,idx
!
  nn=sqrt(8*size(vector)+1)
  nn=(nn-1)/2
  do ii=1,M
  do jj=1,ii
    idx=ii+(2*M-jj)*(jj-1)/2
    if (ii.ne.jj) Vector(idx)=Vector(idx)/(2.0d0)
  enddo
  enddo
!
  return;end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
