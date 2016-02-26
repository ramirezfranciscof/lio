!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine dft_get_qm_forces(dxyzqm)
!--------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod,only:natom,nsol,cubegen_only,qm_forces_ds,
     >                     do_ehrenfest,first_step
       implicit none
       real*8,intent(out) :: dxyzqm(3,natom)
       real*8,allocatable :: ff1G(:,:),ffSG(:,:),ff3G(:,:)
       real*8             :: factor
       integer            :: fileunit,kk,ii,igpu
       logical            :: print_forces

!--------------------------------------------------------------------!
       if(cubegen_only) return
       call g2g_timer_sum_start('Forces')
       allocate(ff1G(natom,3),ffSG(natom,3),ff3G(natom,3))

       call g2g_timer_start('int1G')
       ff1G=0.0d0
       call aint_query_gpu_level(igpu)
       if (igpu.lt.4) then
         call g2g_timer_sum_start('Nuclear attraction gradients')
         call int1G(ff1G)
         call g2g_timer_sum_stop('Nuclear attraction gradients')
       elseif (nsol.le.0) then
         call g2g_timer_sum_start('Nuclear attraction gradients')
         call int1G(ff1G)
         call aint_qmmm_forces(ff1G,0)
         call g2g_timer_sum_stop('Nuclear attraction gradients')
       endif
       call g2g_timer_stop('int1G')
       do ii=1,natom
         kk=610+ii
         write(kk,101) ff1G(ii,:)
       enddo
!       write(671,*) ''


       call g2g_timer_start('intSG')
       call g2g_timer_sum_start('Overlap gradients')
       ffSG=0.0d0
       call intSG(ffSG)
!       write(666,101) ffSG
!       write(666,*) ''
!       if (do_ehrenfest) then
!         ffSG=-transpose(qm_forces_ds)
!       endif
       do ii=1,natom
         kk=620+ii
         write(kk,101) ffSG(ii,:)
       enddo
!       write(672,*) ''
101    format(3(3X,E20.13))
       call g2g_timer_stop('intSG')
       call g2g_timer_sum_stop('Overlap gradients')

       call g2g_timer_start('int3G')
       call g2g_timer_sum_start('Coulomb+Exchange-correlation')
       ff3G=0.0d0
       call int3G(ff3G,.true.)
       call g2g_timer_stop('int3G')
       call g2g_timer_sum_stop('Coulomb+Exchange-correlation')
       do ii=1,natom
         kk=630+ii
         write(kk,101) ff3G(ii,:)
       enddo
!       write(673,*) ''


       factor=1.D0
c       factor=627.509391D0/0.5291772108D0
       do kk=1,natom
       do ii=1,3
         dxyzqm(ii,kk)=ff1G(kk,ii)+ffSG(kk,ii)+ff3G(kk,ii)
         dxyzqm(ii,kk)=dxyzqm(ii,kk)*factor
!         dxyzqm(ii,kk)=0.0d0
       enddo
       enddo

       do ii=1,natom
         kk=660+ii
         write(kk,101) dxyzqm(:,ii)
       enddo


!--------------------------------------------------------------------!
       print_forces=.true.
       if (print_forces) then
         fileunit=3242
         open(unit=fileunit,file='Forces.log',access='APPEND')

         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         do kk=1,natom
           write(fileunit,200) 'TOTS',kk,
     >       ff1G(kk,1)+ffSG(kk,1)+ff3G(kk,1),
     >       ff1G(kk,2)+ffSG(kk,2)+ff3G(kk,2),
     >       ff1G(kk,3)+ffSG(kk,3)+ff3G(kk,3)
         enddo
         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         write(fileunit,*)

         do kk=1,natom
           write(fileunit,200) 
     >       'FF1G',kk,ff1G(kk,1),ff1G(kk,2),ff1G(kk,3)
           write(fileunit,200)
     >       'FFSG',kk,ffSG(kk,1),ffSG(kk,2),ffSG(kk,3)
           write(fileunit,200)
     >       'FF3G',kk,ff3G(kk,1),ff3G(kk,2),ff3G(kk,3)
           write(fileunit,*)
         enddo
         write(fileunit,*)

         close(fileunit)
       endif

       if (nsol.le.0) then
         call g2g_timer_sum_stop('Forces')
         call g2g_timer_sum_stop("Total")
         call g2g_timer_summary()
         call g2g_timer_clear()
       endif

       if (first_step) first_step=.false.

!--------------------------------------------------------------------!
       deallocate(ff1G,ffSG,ff3G)
 200   format(1X,A4,1X,I4,3(2X,E14.7))
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
