      subroutine dft_get_mm_forces(dxyzcl,dxyzqm)
      use garcha_mod
c      use qmmm_module, only : qmmm_struct
      implicit real*8 (a-h,o-z)
       REAL*8 , intent(inout) :: dxyzqm(3,natom)
       REAL*8 , intent(inout) :: dxyzcl(3,nsol)
       real*8, dimension (:,:), ALLOCATABLE :: ff,ffcl
c       real*8, dimension (:,:), ALLOCATABLE :: ffs,ffcls
       allocate(ff(natom,3), ffcl(ntatom,3))
c       allocate(ffs(natom,3), ffcls(ntatom,3))
c       real*8 ftot(3)

         ffcl=0
         ff=0
        if(noconverge.eq.0) then
        call g2g_timer_start('intsolG')
        call intsolG(ff,ffcl)
        call g2g_timer_stop('intsolG')
        endif
       do i=1,natom 
        do j=1,3
       dxyzqm(j,i)=ff(i,j)+dxyzqm(j,i)
        enddo
         enddo
       do jj=1,nsol
        do j=1,3
        dxyzcl(j,jj)=ffcl(natom+jj,j) !+dxyzcl(j,jj)
       enddo
       enddo         
       
       deallocate (ff,ffcl) 
897    format (F17.11)

      end
c---------------------------------------------------------------------
