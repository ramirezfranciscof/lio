!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine add_direct( this, matin, rcoef, icoef, stat_e )

   implicit none
   class(liomat)      , intent(inout) :: this
   class(liomat)      , intent(in)    :: matin
   real(wpk)          , intent(in)    :: rcoef
   real(wpk), optional, intent(in)    :: icoef
   integer  , optional, intent(inout) :: stat_e
   character(len=*)   , parameter     :: myname="add_direct"


   if ( this%my_size == 0 ) then
      call check_stat( myname, 14, 1, stat_e )
   end if

   if ( this%my_size /= matin%my_size ) then
      call check_stat( myname, 18, 2, stat_e )
   end if


   this%rvals(:,:) = this%rvals(:,:) + rcoef * matin%rvals(:,:)
   this%ivals(:,:) = this%ivals(:,:) + rcoef * matin%ivals(:,:)

   if ( present(icoef) ) then
      this%rvals(:,:) = this%rvals(:,:) - icoef * matin%ivals(:,:)
      this%ivals(:,:) = this%ivals(:,:) + icoef * matin%rvals(:,:)
   end if

end subroutine add_direct
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
