!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine setup_r( this, inpmat )

   implicit none
   character(len=40), parameter     :: myname="setup_r"
   class(liomat)    , intent(inout) :: this
   real(wpk)        , intent(in)    :: inpmat(:,:)
   integer                          :: mystat

   if ( size(inpmat,1) /= this%my_size ) then 
      call check_stat( myname, 53, 1)
   end if

   if ( size(inpmat,2) /= this%my_size ) then
      call check_stat( myname, 53, 1)
   end if

   this%rvals(:,:) = inpmat(:,:)

end subroutine setup_r
!
!
!------------------------------------------------------------------------------!
subroutine setup_i( this, inpmat )

   implicit none
   character(len=40), parameter     :: myname="setup_i"
   class(liomat)    , intent(inout) :: this
   real(wpk)        , intent(in)    :: inpmat(:,:)
   integer                          :: mystat

   if ( size(inpmat,1) /= this%my_size ) then
      call check_stat( myname, 53, 1)
   end if

   if ( size(inpmat,2) /= this%my_size ) then
      call check_stat( myname, 53, 1)
   end if

   this%ivals(:,:) = inpmat(:,:)

end subroutine setup_i
!
!
!------------------------------------------------------------------------------!
subroutine setup_c( this, inpmat_r, inpmat_i )

   implicit none
   character(len=40), parameter     :: myname="setup_c"
   class(liomat)    , intent(inout) :: this
   real(wpk)        , intent(in)    :: inpmat_r(:,:)
   real(wpk)        , intent(in)    :: inpmat_i(:,:)
   integer                          :: mystat

   if ( size(inpmat,1) /= this%my_size ) then
      call check_stat( myname, 53, 1)
   end if

   if ( size(inpmat,2) /= this%my_size ) then
      call check_stat( myname, 53, 1)
   end if

   this%rvals(:,:) = inpmat_r(:,:)
   this%ivals(:,:) = inpmat_i(:,:)

end subroutine setup_c
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
