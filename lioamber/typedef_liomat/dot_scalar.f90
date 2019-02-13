!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dot_scalar( this, rval, ival, stat_e )

   implicit none
   class(liomat)      , intent(inout) :: this
   real(wpk)          , intent(in)    :: rval
   real(wpk)          , intent(in)    :: ival
   integer  , optional, intent(inout) :: stat_e
   character(len=*)   , parameter     :: myname="dot_scalar"

   if ( this%my_size == 0 ) then
      call check_stat( myname, 12, 1, stat_e )
   end if

   do jj = 1, this%my_size
   do ii = 1, this%my_size
      this%rvals(ii,jj) = rval * this%rvals(ii,jj) - ival * this%ivals(ii,jj)
      this%ivals(ii,jj) = ival * this%rvals(ii,jj) + rval * this%ivals(ii,jj)
   end do
   end do

end subroutine dot_scalar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
