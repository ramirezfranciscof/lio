!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine find_free_unit(file_unit)
!--------------------------------------------------------------------!
  implicit none
  integer,intent(out) :: file_unit
  integer             :: errorid
  logical             :: keep_looking
  integer,parameter   :: first_unit=100
  integer,parameter   :: limit_unit=10000


  file_unit=first_unit
  keep_looking=.true.
  do while (keep_looking)
     file_unit=file_unit+1
     inquire(unit=file_unit,opened=keep_looking,iostat=errorid)
     if (errorid.ne.0)             keep_looking=.false.
     if (file_unit.ge.limit_unit)  keep_looking=.false.
  enddo


  if (file_unit.ge.limit_unit) then
     print*,'ALL UNITS ARE OCCUPIED'
     stop
  endif


  if (errorid.ne.0) then
     print*,'AN ERROR OCCURRED ON UNIT ',fileunit
     print*,'(iostat= ',errorid,' )'
     stop
  endif


  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
