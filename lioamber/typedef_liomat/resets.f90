!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine unset_0( this )

   implicit none
   character(len=40), parameter     :: myname="unset_0"
   class(liomat)    , intent(inout) :: this

   this%my_size = 0
   if ( allocated(this%rvals) ) deallocate(this%rvals)
   if ( allocated(this%ivals) ) deallocate(this%ivals)

end subroutine unset_0
!
!
!------------------------------------------------------------------------------!
subroutine reset_0( this, set_size )

   implicit none
   character(len=40), parameter     :: myname="reset_0"
   class(liomat)    , intent(inout) :: this
   integer          , intent(in)    :: set_size

   call this%unset_0()

   this%my_size = set_size
   allocate( this%rvals(set_size, set_size) )
   allocate( this%ivals(set_size, set_size) )

   this%rvals(:,:) = 0.0_wpk
   this%ivals(:,:) = 0.0_wpk
 
end subroutine reset_0
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
