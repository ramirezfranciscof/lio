!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine add_commut( this, mata, matb, rcoef, icoef, stat_e )

   implicit none
   class(liomat)      , intent(inout) :: this
   class(liomat)      , intent(in)    :: mata
   class(liomat)      , intent(in)    :: matb
   real(wpk)          , intent(in)    :: rcoef
   real(wpk), optional, intent(in)    :: icoef
   integer  , optional, intent(inout) :: stat_e
   character(len=*)   , intent(in)    :: myname="add_commut"

   if ( this%my_size == 0 ) then
      call check_stat( myname, 14, 1, stat_e )
   end if

   if ( this%my_size /= mata%my_size ) then
      call check_stat( myname, 18, 2, stat_e )
   end if

   if ( this%my_size /= matb%my_size ) then
      call check_stat( myname, 22, 2, stat_e )
   end if

   call this%add_matmul( mata, matb,  rcoef,  icoef, stat_e )
   call this%add_matmul( matb, mata, -rcoef, -icoef, stat_e )

end subroutine add_commut
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
