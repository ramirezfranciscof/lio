!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine add_matmul( this, mata, matb, rcoef, icoef, stat_e )
   implicit none
   class(liomat)      , intent(inout) :: this
   class(liomat)      , intent(in)    :: mata
   class(liomat)      , intent(in)    :: matb
   real(wpk)          , intent(in)    :: rcoef
   real(wpk), optional, intent(in)    :: icoef
   integer  , optional, intent(inout) :: stat_e
   character(len=*)   , intent(in)    :: myname="add_matmul"


   if ( this%my_size == 0 ) then
      call check_stat( myname, 15, 1, stat_e )
   end if

   if ( this%my_size /= mata%my_size ) then
      call check_stat( myname, 19, 2, stat_e )
   end if

   if ( this%my_size /= matb%my_size ) then
      call check_stat( myname, 23, 2, stat_e )
   end if


   this%rvals = this%rvals + rcoef * matmul( mata%rvals, matb%rvals )
   this%ivals = this%ivals + rcoef * matmul( mata%rvals, matb%ivals )
   this%ivals = this%ivals + rcoef * matmul( mata%ivals, matb%rvals )
   this%rvals = this%rvals - rcoef * matmul( mata%ivals, matb%ivals )

   if ( present( icoef ) ) then
      this%ivals = this%rvals + icoef * matmul( mata%rvals, matb%rvals )
      this%rvals = this%ivals - icoef * matmul( mata%rvals, matb%ivals )
      this%rvals = this%ivals - icoef * matmul( mata%ivals, matb%rvals )
      this%ivals = this%rvals - icoef * matmul( mata%ivals, matb%ivals )
   end if

end subroutine add_matmul
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
