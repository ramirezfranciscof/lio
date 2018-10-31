!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module error_log
!------------------------------------------------------------------------------!
!
!  This module was designed as a simple way to keep track of errors and to
!  help identify their cause and exact location in the code. To do so, it
!  offers two subroutines:
!
!
!  (1) check_stat( subname, subline, statidn, passidn )
!
!      Call this procedure to check if <statidn> is equal to zero. If it
!      isn't, then it will record the <subname> and <subline> into the
!      internal data arrays of the log. Then it will check the presence
!      of passidn to decide whether to print the log of errors and stop
!      the program or return control to the caller procedure. If the
!      control returns to caller, it will do so by copying the statidn
!      in the passidn.
!
!
!  (2) print_log()
!
!      This procedure prints to screen the list of all errors recorded in
!      the log, from the first one found to the last one.
!
!
!------------------------------------------------------------------------------!
!
!  The module keeps an internal record of all errors registered by the use
!  of the procedure <check_stat> inside the following variables:
!
!     logdata_subnames: All the name of the caller procedures that detected 
!                       the error.
!
!     logdata_sublines: All the lines of the corresponding procedures where
!                       the error was detected.
!
!     logdata_statidns: All statidns registered.
!
!
!------------------------------------------------------------------------------!
!
!  It is worth noting that the way in which check_stat was implemented, allows
!  it to receive an external stat that is itself optional for the caller sub,
!  further pushing the responsability of handling the error further down the
!  calling chain. Check the following example, which will return control to
!  the <example_sub> if the <extern_stat> is present (which in return will do
!  the same to its caller) and will print and stop if it is not. Also, if the
!  control is returned, the <extern_stat> will automatically have the same
!  value as the <intern_stat>.
!
!
!     subroutine example_sub( ... argsA1 ..., extern_stat )
!        (...)
!        use error_log, only: check_stat
!        (...)
!        integer, intent(out), optional :: extern_stat
!        integer                        :: intern_stat
!        (...)
!        call internal_sub( ...args... , intern_stat )
!   =>   call check_stat( "example_sub", line_number, intern_stat, extern_stat)
!   =>   if ( intern_stat /= 0 ) return
!        (...)
!     end subroutine
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   private
   integer                   , parameter   :: SUBNAME_LEN = 80

   integer                                 :: logdata_size = 0
   character(len=SUBNAME_LEN), allocatable :: logdata_subnames(:)
   integer                   , allocatable :: logdata_sublines(:)
   integer                   , allocatable :: logdata_statidns(:)

   public :: check_stat
   public :: print_log

contains
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_stat( subname, subline, statidn, passidn )

   implicit none
   character(len=*), intent(in)            :: subname
   integer,          intent(in)            :: subline
   integer,          intent(in)            :: statidn
   integer,          intent(out), optional :: passidn

   integer                                 :: nn, pi, pf
   character(len=SUBNAME_LEN), allocatable :: logtemp_subnames(:)
   integer                   , allocatable :: logtemp_sublines(:)
   integer                   , allocatable :: logtemp_statidns(:)


   if ( present(passidn) ) passidn = statidn
   if ( statidn == 0 ) return


   if ( logdata_size > 0 ) then
      allocate( logtemp_subnames(logdata_size) )
      allocate( logtemp_sublines(logdata_size) )
      allocate( logtemp_statidns(logdata_size) )

      logtemp_subnames = logdata_subnames
      logtemp_sublines = logdata_sublines
      logtemp_statidns = logdata_statidns

      deallocate( logdata_subnames )
      deallocate( logdata_sublines )
      deallocate( logdata_statidns )
   end if

   logdata_size = logdata_size + 1
   allocate( logdata_subnames(logdata_size) )
   allocate( logdata_sublines(logdata_size) )
   allocate( logdata_statidns(logdata_size) )

   do nn = 1, logdata_size-1
      logdata_subnames(nn) = logtemp_subnames(nn)
      logdata_sublines(nn) = logtemp_sublines(nn)
      logdata_statidns(nn) = logtemp_statidns(nn)
   end do

   pi = 1
   pf = min( SUBNAME_LEN, len(subname) )
   logdata_subnames(logdata_size) = subname( pi : pf )
   logdata_sublines(logdata_size) = subline
   logdata_statidns(logdata_size) = statidn


   if ( .not.present(passidn) ) then
      print "(A)", ""
      print "(A)", "CRITICAL ERROR DETECTED!"
      call print_log()
      print "(A)", "Aborting run..."
      print "(A)", ""
      stop

   end if

end subroutine check_stat
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine print_log()

   implicit none
   integer :: nn

   print "(A)", "Printing log of errors from first found to last..."
   print "(A)", ""
   do nn = 1, logdata_size
      print "(A,I6,A,I6,2A)", " => stat=", logdata_statidns(nn) &
                           &, " found in line ", logdata_sublines(nn) &
                           &, " of procedure ", adjustl(logdata_subnames(nn))
   end do
   print "(A)", ""
   print "(A)", "If you don't know and can't find out what caused any of this"
   print "(A)", "these errors, please report it back to us through our github"
   print "(A)", "website: https://github.com/MALBECC/lio"
   print "(A)", ""

end subroutine print_log
!
!
!
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
