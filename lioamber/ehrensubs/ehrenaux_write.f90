!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! EHRENFEST WRITE SUBROUTINES
!
! * ehren_write_dip (handles dipole moment printing)
! * ehren_write_for (handles grandient printing to output)
! * ehren_write_pop (handles population/charge printing to output)
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehren_write_dip( time, dipole, dipmod, unitid, header )
   implicit none
   real*8 , intent(in) :: time
   real*8 , intent(in) :: dipole(3)
   real*8 , intent(in) :: dipmod
   integer, intent(in) :: unitid
   logical, intent(in) :: header
    

    open( unit=unitid, file="dipole_moment" )
    if (header) then
       write(unit=unitid, fmt=100) ""
       write(unit=unitid, fmt=100) &
       &'#DIPOLE MOMENT, X Y Z COMPONENTS AND NORM (DEBYES)'
       write(unit=unitid, fmt=100) ""
    else
       write(unit=unitid, fmt=200) &
       & time, dipole(1), dipole(2), dipole(3), dipmod
    endif

100 format(A)
200 format(4x, 5(F12.8, 2x))
end subroutine ehren_write_dip


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehren_write_for(natoms, offset, time, forces, unitid)
   implicit none
   integer, intent(in) :: natom
   integer, intent(in) :: offset
   real*8 , intent(in) :: time
   real*8 , intent(in) :: forces(3, natom+offset)
   integer, intent(in) :: unitid

   integer             :: kk

   write(unit=unitid, fmt=100) "Time = ",time
   do kk=offset+1, offset+natom
      write(unit=unitid, fmt=200) kk, dxyz(1, kk), dxyz(2, kk), dxyz(3, kk)
   enddo

100 format (A,    1x, F12.8)
200 format (I6, 3(2x, f10.6))
end subroutine ehren_write_for


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
