!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      Subroutine SCF_in(E,qmcoords,qmvels,clcoords,nsolin,dipxyz)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      use garcha_mod
      REAL*8 , intent(in)  :: qmcoords(3,natom)
      REAL*8 , intent(in)  :: qmvels(3,natom)
      REAL*8 , intent(in)  :: clcoords(4,nsolin)
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
          nucpos(kk,ii) = r(ii,kk)
          nucvel(kk,ii) = qmvels(kk,ii)
c          write(89,*) ii, kk, qmcoords(ii,kk)
c          write(87,*) ii, kk, r(ii,kk)
       enddo
       write(18,345) Iz(ii),qmcoords(:,ii)
      enddo

!--------------------------------------------------------------------!
! I am not sure this should be here, but it is the only
! place to put it right now (FFR)
       call liomain()
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
!--------------------------------------------------------------------!
      if (do_ehrenfest) then

        if (first_step) then
          call SCF(E,dipxyz)
        else
          call SCF(E,dipxyz)
          do ii=1,M
          do jj=1,M
            write(666,*) RhoSave(ii,jj),RhoCero(ii,jj)
!     >                     ABS(RhoSave(ii,jj)+RhoCero(ii,jj)),
!     >                     ABS(RhoSave(ii,jj)+RhoCero(ii,jj))/2
          enddo
          enddo
          write(666,*) ''
          write(666,*) ''
          call ehrendyn(E,dipxyz)
        endif
      else
        if (OPEN) then 
          call SCFOP(E,dipxyz)
        else
          call SCF(E,dipxyz)
        endif
      endif

!
!--------------------------------------------------------------------!
 345  format(2x,I2,2x,3(f10.6,2x))
      return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
