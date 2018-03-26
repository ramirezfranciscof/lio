!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_ener( Msize, Hmat, Pmat, Energy )
!------------------------------------------------------------------------------!
!
!  This subroutine simply calculates the energy associated to a given core
!  fock. This core fock matrix can be calculated independently from the 
!  density matrix and only depends on the position of the nuclei, but the
!  energy does depend on the electrons so it must be re-calculated more
!  frequently than the core fock matrix.
!
!  Depending if this Hmat core matrix contains just one part (1e, ecp, solv)
!  or many of them, so will the energy be associated with only one of those
!  parts or many.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(in)  :: Msize
   real*8 , intent(in)  :: Hmat( Msize, Msize )
   real*8 , intent(in)  :: Pmat( Msize, Msize )
   real*8 , intent(out) :: Energy
   integer              :: ii, jj

   Energy = 0.0d0
   do jj = 1, Msize
   do ii = 1, Msize
      Energy = Energy + Hmat(ii,jj) * Pmat(ii,jj)
   end do
   end do

end subroutine rmmCalc_ener
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
