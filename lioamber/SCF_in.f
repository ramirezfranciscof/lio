!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      Subroutine SCF_in(E,qmcoords,qmvels,clcoords,nsolin,dipxyz)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      use garcha_mod
      use ehrenfest
      REAL*8 , intent(in)    :: qmcoords(3,natom)
      REAL*8 , intent(in)    :: qmvels(3,natom)
      REAL*8 , intent(in)    :: clcoords(4,nsolin)
      REAL*8 , intent(inout) :: E,dipxyz(3)
      integer :: nn,kk,ii,jj
      nsol=nsolin
      ntatom=nsol+natom

      call g2g_timer_sum_start("Total")
!      do kk=1,natom
!        write(666,123) qmvels(1:3,kk)
!      enddo
! 123  FORMAT(3(3X,f12.8))

      deallocate (r,v,Em,Rm,pc,nnat)

      allocate (r(ntatom,3),v(ntatom,3),Em(ntatom)
     >,Rm(ntatom),pc(ntatom))
      allocate (nnat(100))
       ngDyn=natom*ng0
       
      if(writexyz) write(18,*) natom
      if(writexyz) write(18,*)

      do ii=1,nsol
        nn=natom+ii
        pc(nn)=clcoords(4,ii)

        do kk=1,3
          r(nn,kk)=clcoords(kk,ii)/0.529177D0
        enddo
c        write(18,345) 8,r(nn,1),r(nn,2),r(nn,3)
      enddo


      if (allocated(nucpos)) deallocate(nucpos)
      if (allocated(nucvel)) deallocate(nucvel)
      allocate(nucpos(3,natom))
      allocate(nucvel(3,natom))
      do ii=1,natom
       do kk=1,3
          r(ii,kk)   = qmcoords(kk,ii)/0.529177D0
          rqm(ii,kk) = qmcoords(kk,ii)/0.529177D0
! velocity units in angstrom per 1/20.455 pico-second must go to atomic units
          nucpos(kk,ii) = r(ii,kk)
          nucvel(kk,ii) = qmvels(kk,ii)*(20.455)*(2.418884326505E-5)
          nucvel(kk,ii) = nucvel(kk,ii)/(0.529177d0)
       enddo
       write(18,345) Iz(ii),qmcoords(:,ii)
      enddo

!%%%%-FFR-START-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! I am not sure this should be here, but it is the only
! place to put it right now
       if (allocated(atom_mass)) deallocate(atom_mass)
       allocate(atom_mass(natom))
       call ehren_masses(natom,Iz,atom_mass)

       call liomain()
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
!--------------------------------------------------------------------!
      if (do_ehrenfest) then
        if (first_step) then
          call SCF(E,dipxyz)
        endif
        call ehrendyn(E,dipxyz)
      else
        if (OPEN) then 
          call SCFOP(E,dipxyz)
        else
          call SCF(E,dipxyz)
        endif
      endif
! OLD CODE:
!       if (OPEN) then
!          call SCFOP(E,dipxyz)
!       else
!          call SCF(E,dipxyz)
!       endif
!%%%%-FFR-START-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!
!--------------------------------------------------------------------!
 345  format(2x,I2,2x,3(f10.6,2x))
      return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
