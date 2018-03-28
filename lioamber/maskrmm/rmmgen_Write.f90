!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmgen_Write( funit_out, rmm_pointer_nid, matinp, vecinp, ndfact )
   use garcha_mod, only: M, Md, NCO, Nunp, natom, Nang, RMM, rhoalpha, rhobeta
   implicit none
   integer, intent(in)           :: funit_out
   integer, intent(in)           :: rmm_pointer_nid
   real*8 , intent(in), optional :: matinp(:,:), vecinp(:), ndfact

   integer :: M1, M3, M5, M7, M9, M11, M13, M15, M17, M18
   integer :: M18b, M19, M20, M21, M22, M23
   integer :: nn_i, nn_f, nn, ii, jj
   integer :: size_compare, rmm_size
   integer :: MM, MMd, NCOa, NCOb
   real*8  :: fact
   real*8, allocatable :: rmmcopy(:)

   MM   = M  * (M+1)  / 2
   MMd  = Md * (Md+1) / 2
   NCOa = NCO
   NCOb = NCO + Nunp

   M1   = 1                    ! Density Matrix
   M3   = M1   + MM            ! Density Matrix New
   M5   = M3   + MM            ! Overlap or Fock[a] matrices
   M7   = M5   + MM            ! Density Fitting G
   M9   = M7   + MMd           ! Density Fitting Gm
   M11  = M9   + MMd           ! Fock Core Matrix H
   M13  = M11  + MM            ! Eigenvalues[a] +  space used in least-squares
   M15  = M13  + M             ! Aux vector for ESSl
   M17  = M15  + MM            ! Least squares
   M18  = M17  + MMd           ! vectors of MO[a]
   M18b = M18  + M*NCOa        ! vectors of MO[b]
   M19  = M18b + M*NCOb        ! Weights
   M20  = M19  + natom*50*Nang ! New Fock matrix[a]
   M21  = M20  + MM            ! New Fock matrix[b]
   M22  = M21  + MM            ! Eigenvalues[b]
   M23  = M22  + M             ! RAM storage of two-electron integrals (MEMO=T)

   if ( rmm_pointer_nid < 0 ) then
   end if


   select case(rmm_pointer_nid)
      case(-1)
         write(unit=funit_out, fmt=*) "Writing Alpha Density Matrix"
         nn_i = 1
         nn_f = MM
      case(-2)
         write(unit=funit_out, fmt=*) "Writing Beta Matrix"
         nn_i = 1
         nn_f = MM
      case(1)
         write(unit=funit_out, fmt=*) "Writing RMM(M1): Density Matrix"
         nn_i = M1
         nn_f = M3-1
      case(3)
         write(unit=funit_out, fmt=*) "Writing RMM(M3): Density Matrix New"
         nn_i = M3
         nn_f = M5-1
      case(5)
         write(unit=funit_out, fmt=*) "Writing RMM(M5): Overlap or Fock[a] matrices"
         nn_i = M5
         nn_f = M7-1
      case(7)
         write(unit=funit_out, fmt=*) "Writing RMM(M7): Density Fitting G"
         nn_i = M7
         nn_f = M9-1
      case(9)
         write(unit=funit_out, fmt=*) "Writing RMM(M9): Density Fitting Gm"
         nn_i = M9
         nn_f = M11-1
      case(11)
         write(unit=funit_out, fmt=*) "Writing RMM(M11): Fock Core Matrix H"
         nn_i = M11
         nn_f = M13-1
      case(13)
         write(unit=funit_out, fmt=*) "Writing RMM(M13): Eigenvalues[a] +  space used in least-squares"
         nn_i = M13
         nn_f = M15-1
      case(15)
         write(unit=funit_out, fmt=*) "Writing RMM(M15): Aux vector for ESSl"
         nn_i = M15
         nn_f = M17-1
      case(17)
         write(unit=funit_out, fmt=*) "Writing RMM(M17): Least squares"
         nn_i = M17
         nn_f = M18-1
      case(18)
         write(unit=funit_out, fmt=*) "Writing RMM(M18): vectors of MO[a]"
         nn_i = M18
         nn_f = M18b-1
      case(181)
         write(unit=funit_out, fmt=*) "Writing RMM(M18b): vectors of MO[b]"
         nn_i = M18b
         nn_f = M19-1
      case(19)
         write(unit=funit_out, fmt=*) "Writing RMM(M19): Weights"
         nn_i = M19
         nn_f = M20-1
      case(20)
         write(unit=funit_out, fmt=*) "Writing RMM(M20): New Fock matrix[a]"
         nn_i = M20
         nn_f = M21-1
      case(21)
         write(unit=funit_out, fmt=*) "Writing RMM(M21): New Fock matrix[b]"
         nn_i = M21
         nn_f = M22-1
      case(22)
         write(unit=funit_out, fmt=*) "Writing RMM(M22): Eigenvalues[b]"
         nn_i = M22
         nn_f = M23-1
      case(23)
         write(unit=funit_out, fmt=*) "Writing RMM(M23): RAM storage of two-electron integrals"
         nn_i = M23
         nn_f = size(RMM)
      case default
         write(unit=funit_out, fmt=*) "Writing all RMM"
         nn_i = 1
         nn_f = size(RMM)
   end select

   rmm_size = nn_f - nn_i + 1
   allocate( rmmcopy(rmm_size) )
   if ( rmm_pointer_nid >= 0 ) then
      do nn = 1, rmm_size
         rmmcopy(nn) = RMM(nn_i+nn-1)
      end do
   else if ( rmm_pointer_nid == -1 ) then
      do nn = 1, rmm_size
         rmmcopy(nn) = rhoalpha(nn)
      end do
   else if ( rmm_pointer_nid == -2 ) then
      do nn = 1, rmm_size
         rmmcopy(nn) = rhobeta(nn)
      end do
   else
      write(unit=funit_out, fmt=*) "Bad rmm_pointer_nid...", rmm_pointer_nid
      return
   end if

   if ( rmm_pointer_nid >= 0 ) then
   if ( nn_f > size(RMM) ) then
      write(unit=funit_out, fmt=*) "RMM is not that big..."
      write(unit=funit_out, fmt=*) "Maybe you requested an OS pointer in a CS calc?"
      return
   end if
   end if

   if ( present(matinp) ) then
      size_compare = size(matinp,1) * (size(matinp,2)+1) / 2
      if ( size_compare /= rmm_size ) then
         write(unit=funit_out, fmt=*) "I cant compare this matrix..."
         write(unit=funit_out, fmt=*) "matinp 1:       ", size(matinp,1)
         write(unit=funit_out, fmt=*) "matinp 2:       ", size(matinp,1)
         write(unit=funit_out, fmt=*) "matinp MM:      ", size_compare
         write(unit=funit_out, fmt=*) "RMM requested: ", rmm_size
         return
      end if
   end if

   if ( present(vecinp) ) then
      if ( size(vecinp) /= rmm_size ) then
         write(unit=funit_out, fmt=*) "I cant compare this vector..."
         write(unit=funit_out, fmt=*) "vecinp:         ", size(vecinp)
         write(unit=funit_out, fmt=*) "RMM requested: ", rmm_size
         return
      end if
   end if

   ii = 0
   jj = 1
   do nn = nn_i, nn_f

      if ( present(matinp) ) then
         ii = ii + 1
         if ( ii > size(matinp,1) ) then
            ii = 1
            jj = jj+1
         end if
         fact = 1.0d0
         if ( present(ndfact) .and. (ii /= jj) ) fact=ndfact
      end if

      if ( present(matinp) .and. present(vecinp) ) then
            write(unit=funit_out, fmt=*) &
            & nn, rmmcopy(nn), vecinp(nn), &
            & ii, jj, matinp(ii,jj)*fact, matinp(jj,ii)*fact
      else
         if ( present(matinp) ) then
            write(unit=funit_out, fmt=*) &
            & nn, rmmcopy(nn), &
            & ii, jj, matinp(ii,jj)*fact, matinp(jj,ii)*fact
         else if ( present(vecinp) ) then
            write(unit=funit_out, fmt=*) &
            & nn, rmmcopy(nn), vecinp(nn)
         else
            write(unit=funit_out, fmt=*) &
            & nn, rmmcopy(nn)
         end if
      end if

   end do

end subroutine rmmgen_Write
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
