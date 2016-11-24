! MAIN SUBROUTINE ----------------------------------------------------
! DFT calculation with gaussian basis sets
!---------------------------------------------------------------------
      program liosolo

      use garcha_mod, only : natom, nsol, OPEN, NMAX, rqm, pc, Smat, RealRho,  &
                             Iz, r, a0, allnml, basis, basis_set, cube_dens,   &
                             cube_dens_file, cube_elec, cube_elec_file,        &
                             cube_orb, cube_orb_file, cube_res, cube_sel,      &     
                             cubegen_only, dens, dgtrig, DIIS, epsilon, etold, &
                             exter, field, fitting_set, frestartin, frestart,  &
                             fx, fy, fz, gold, good_cut, hybrid_converg,       &
                             idip, iexch, igrid, igrid2, int_basis, integ,     &
                             intsoldouble, m, nbch, ntatom,ntdstep, nunp,      &
                             omit_bas, predcoef, propagator, rmax, rmaxs,      &
                             style, tdrestart, tdstep, timedep, told, vcinp,   &
                             verbose, writedens, writeforces, writexyz, ndiis, &
                             driving_rate,transport_calc,                      &
                             generate_rho0,save_charge_freq ! agregadas las cosas de transporte. Juanderboy.

      use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,     &
                             local_nonlocal, ecp_debug, ecp_full_range_int,    &
                             verbose_ECP, FOCK_ECP_read, FOCK_ECP_write,       &
                             Fulltimer_ECP, cut2_0, cut3_0
      use mathsubs
#ifdef CUBLAS
      use cublasmath
#endif
      implicit none 
      character(len=20)   :: argument, inpfile, inpbasis, inpcoords
      integer             :: charge, ifind, ierr, i, k, narg, ios
      real*8              :: dipxyz(3), escf
      logical             :: filexist
      real*8, allocatable :: dxyzqm(:,:), dxyzcl(:,:)
      namelist /lio/ natom, nsol, charge, OPEN, NMAX, Nunp, VCINP, frestartin, &
                     frestart, GOLD, told, Etold, rmax, rmaxs, predcoef, idip, &
                     writexyz, intsoldouble, DIIS, ndiis, dgtrig, Iexch, integ,&
                     dens, igrid, igrid2, timedep, tdstep, ntdstep, propagator,&
                     NBCH, field, a0, epsilon, exter, Fx, Fy, Fz, tdrestart,   &
                     writedens, basis_set, fitting_set, int_basis, writeforces,&
                     ecpmode, ecptypes, tipeECP, ZlistECP, local_nonlocal,     &
                     cutECP, ecp_debug, ecp_full_range_int, verbose_ECP,       &
                     verbose, FOCK_ECP_read, FOCK_ECP_write, Fulltimer_ECP,    &
                     cut2_0, cut3_0, hybrid_converg, good_cut, style, allnml,  &
                     cube_elec, cube_dens, cube_dens_file, cube_orb_file,      &
                     cube_orb, cube_sel, cubegen_only, cube_res,               &
                     cube_elec_file,driving_rate,transport_calc,               &
                     generate_rho0,save_charge_freq ! agregadas las cosas de transporte. Juanderboy.
!
      ! Calls default values for variables.
      call lio_defaults()
     !defaults
!      basis='input'  ! name of the base file
!      basis_set='DZVP'
!      fitting_set='DZVP Coulomb Fitting'
!      int_basis=.false.
!      cubegen_only=.false.
!      cube_res=40
!      cube_dens=.false.
!      cube_dens_file='dens.cube'
!      cube_orb=.false.
!      cube_sel=0
!      cube_orb_file="orb.cube"
!      cube_elec=.false.
!      cube_elec_file="field.cube"
!      restart_freq=1
!      energy_freq=1 
!      output='output'
!      fcoord='qm.xyz'
!      fmulliken='mulliken'
!      frestart='restart.out'
!      frestartin='restart.in'
!      verbose=.false.
!      OPEN=.false. 
!      NMAX= 100
!      NUNP= 0
!      VCINP= .false.
!      GOLD= 10.
!      told=1.0D-6
!      rmax=16
!      rmaxs=5
!      predcoef=.false.
!      idip=1
!      writexyz=.true.
!      intsoldouble=.true.
!      DIIS=.true.
!      ndiis=30
!      dgtrig=100.
!      Iexch=9
!      integ=.true.
!      DENS = .true.
!      IGRID = 2
!      IGRID2 = 2
!      timedep = 0
!      tdstep = 2.D-3
!      field=.false.
!      a0=1000.0
!      epsilon=1.D0
!      Fx=0.05
!      Fy=0.05
!      Fz=0.05
!      NBCH=10
!      propagator=1
!      tdrestart=.false.
!      writedens=.true.
!      writeforces=.false.
!
! TRANSPORT STUFF
      driving_rate=0.0D0
      transport_calc=.false.
      generate_rho0=.false.
      save_charge_freq=1

 
      ! Reads command line arguments for LIO.
      narg=command_argument_count()
      do i=1, narg
          call get_command_argument(i,argument)
          select case(adjustl(argument))
              case("-i")
                  call get_command_argument(i+1,inpfile)
              case("-b")
                  call get_command_argument(i+1,basis)
              case("-bs")
                 omit_bas=.true.
                 call get_command_argument(i+1,basis_set)
              case("-fs")
                 call get_command_argument(i+1,fitting_set)
              case("-ib")
                 int_basis=.true.
              case("-c")
                 call get_command_argument(i+1,inpcoords)
              case("-v")
                 verbose=.true.
              case default
          endselect
      enddo

      ! Checks input files' existence.
      call g2g_timer_sum_start("Total")

      inquire(file=inpfile,exist=filexist)
      if(filexist) then
          open(unit=100,file=inpfile,iostat=ios)
      else
          write(*,*) 'input file ',adjustl(inpfile),' not found'
          stop
      endif
      read(100,nml=lio,iostat=ierr)
      if(ierr.gt.0) stop 'input error in lio namelist'

      inquire(file=inpcoords,exist=filexist)
      if(filexist) then
          open(unit=101,file=inpcoords,iostat=ios)
      else
          write(*,*) 'input file ',adjustl(inpcoords),' not found'
          stop
      endif
      ! Reads coordinates file.
      ntatom = natom + nsol
      allocate (iz(natom), r(ntatom,3), rqm(natom,3), pc(ntatom))
      do i=1,natom
          read(101,*) iz(i), r(i,1:3)
          rqm(i,1:3) = r(i,1:3)
      enddo
      do i=natom+1,ntatom
          read(101,*) pc(i), r(i,1:3)
      enddo
      r  = r   / 0.529177D0
      rqm= rqm / 0.529177D0
!--------------------------------------------------------
!       call lio_init()   !initialize lio
       call liomain()
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
!--------------------------------------------------------
       if(OPEN) then
         call SCFOP(escf,dipxyz)
       else
         call SCF(escf,dipxyz)
       endif
!--------------------------------------------------------

       write(*,*) 'SCF ENRGY=',escf
!---------------------------------------------------------!
      ! The last argument indicates LIO is being used alone.
      call init_lio_common(natom, Iz, nsol, charge, 0)


 
      call liomain()       !no hace nada!!!!!!
      if (.not.allocated(Smat))    allocate(Smat(M,M))
      if (.not.allocated(RealRho)) allocate(RealRho(M,M))

      if(OPEN) then
          if (ecpmode) stop "ECP is unavailable for Open Shell calculations."         
          call SCFOP(escf,dipxyz)
      else
          call SCF(escf,dipxyz)
      endif


      if(writeforces) then
          if (ecpmode) stop "ECP does not currently feature forces calculation."
 
          open(unit=123,file='forces')

          allocate ( dxyzqm(3, natom) )
          dxyzqm = 0.0

          if(nsol.gt.0) then
              allocate ( dxyzcl(3, natom+nsol) )
              dxyzcl = 0.0
          endif
! juanderboy-------------------------------------------------------------!
       call dft_get_qm_forces(dxyzqm)
       if (nsol.gt.0) then
         call dft_get_mm_forces(dxyzcl,dxyzqm)
       endif
! juanderboy -------------------------------------------------------------!

          call dft_get_qm_forces(dxyzqm)
          if (nsol.gt.0) then
              call dft_get_mm_forces(dxyzcl, dxyzqm)
          endif

          do k=1,natom
              write(123,100) k, dxyzqm(k,1), dxyzqm(k,2), dxyzqm(k,3)
          enddo

          if(nsol.gt.0) then
              do k=natom,natom+nsol
                  write(123,100) k, dxyzcl(k,1), dxyzcl(k,2), dxyzcl(k,3)
              enddo
          endif
       
          deallocate (dxyzqm)
          if(nsol.gt.0) deallocate (dxyzcl)
       endif

       call lio_finalize()

100    format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
       end program

