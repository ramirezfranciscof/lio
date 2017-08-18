!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_updatevels( dtin, natoms, nucmass, nucfors, nucvels )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8,  intent(in)    :: dtin
   integer, intent(in)    :: natoms
   real*8,  intent(in)    :: nucmass(natoms)
   real*8,  intent(in)    :: nucfors(3, natoms)
   real*8,  intent(inout) :: nucvels(3, natoms)

   real*8  :: nucvels_update
   integer :: nn, kk

   do nn = 1, natoms
   do kk = 1, 3
      nucvels_update = dtin * nucfors(kk,nn) / nucmass(nn)
      nucvels(kk,nn) = nucvels(kk,nn) + nucvels_update
   enddo
   enddo

end subroutine ehrenaux_updatevels
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
