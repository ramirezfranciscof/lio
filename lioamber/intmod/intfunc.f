!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine intfunc(ni,nj,ri,rj,ai,aj,samenuc,outint)
!--------------------------------------------------------------------!
!
!
!
!
!--------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: ni,nj       ! Type of orbitals
       real*8,intent(in)  :: ri(3),rj(3) ! Position of nuc
       real*8,intent(in)  :: ai,aj       ! Exponents of Gauss
       logical,intent(in) :: samenuc     ! The nuc is the same
       real*8,intent(out) :: outint(3)   ! Output integral

       real*8     :: IMAT(3,4,4)           ! Integral Matrix
       real*8     :: cij,cta,ct2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       if (samenuc) then
         return
       endif
       call setim(ni,nj+1,ai,aj,ri,rj,IMAT)



! Case <  s | s' >
!--------------------------------------------------------------------!
       if ((ni.eq.0).and.(nj.eq.0)) then
         cta=2*aj

         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,1) ! <s|xs>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,1) ! <s|ys>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,2) ! <s|zs>

       endif


! Case <  p | s' >
!--------------------------------------------------------------------!
       if ((ni.eq.1).and.(nj.eq.0)) then
         cta=2*aj

         orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,1) ! <px|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,1) ! <py|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,1) ! <pz|xs>

         orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,1) ! <px|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,1) ! <py|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,1) ! <pz|ys>

         orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,2) ! <px|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,2) ! <py|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,2) ! <pz|zs>

       endif


! Case <  d | s' >
!--------------------------------------------------------------------!
       if ((ni.eq.2).and.(nj.eq.0)) then
         cta=2*aj

         orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,1) ! <dxx|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,1) ! <dyy|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,1) ! <dzz|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,1) ! <dxy|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,1) ! <dyz|xs>
         orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,1) ! <dzx|xs>

         orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,1) ! <dxx|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,1) ! <dyy|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,1) ! <dzz|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,1) ! <dxy|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,1) ! <dyz|ys>
         orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,1) ! <dzx|ys>

         orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,2) ! <dxx|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,2) ! <dyy|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,2) ! <dzz|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,2) ! <dxy|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,2) ! <dyz|zs>
         orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,2) ! <dzx|zs>

       endif


! Case <  s | p' >
!--------------------------------------------------------------------!
       if ((ni.eq.0).and.(nj.eq.1)) then
         cta=2*aj

         orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,1,1) ! <s|xpx>
         orbint(1)=orbint(1)- 1 *IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,1) !
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <s|xpy> 
         orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <s|xpz>

         orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <s|ypx>
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,1,1) ! <s|ypy>
         orbint(2)=orbint(2)- 1 *IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,1) !
         orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <s|ypz>

         orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <s|zpx>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <s|zpy>
         orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,3) ! <s|zpz>
         orbint(3)=orbint(3)- 1 *IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,1) !

       endif


! Case <  p | p' >
!--------------------------------------------------------------------!
       if ((ni.eq.1).and.(nj.eq.1)) then
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMAT)


            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,1,1) ! <px|xpx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,1,1) ! <py|xpx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,2,1) ! <pz|xpx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <px|xpy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,1,1) ! <py|xpy> 
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,2,1) ! <pz|xpy>

            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <px|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,2) ! <py|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,2) ! <pz|xpz>


            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <px|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,1,1) ! <py|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,2,1) ! <pz|ypx>

            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,1,1) ! <px|ypy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,1,1) ! <py|ypy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,2,1) ! <pz|ypy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <px|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,2) ! <py|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,2) ! <pz|ypz>

            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <px|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,2) ! <py|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,2) ! <pz|zpx>

            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <px|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,2) ! <py|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,2) ! <pz|zpy>

            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,3) ! <px|zpz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,3) ! <py|zpz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,3) ! <pz|zpz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,1) !

       endif


! Case <  d | p' >
!--------------------------------------------------------------------!
       if ((ni.eq.2).and.(nj.eq.1)) then
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMAT)

            orbint(1)=orbint(1)+cta*IMAT(1,3,3)*IMAT(2,1,1)*IMAT(3,1,1) ! <dxx|xpx>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,3,1)*IMAT(3,1,1) ! <dyy|xpx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,3,1) ! <dzz|xpx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,2,1)*IMAT(3,1,1) ! <dxy|xpx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,2,1) ! <dyz|xpx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,2,1) ! <dzx|xpx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,3,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <dxx|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,3,2)*IMAT(3,1,1) ! <dyy|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,3,1) ! <dzz|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,2,2)*IMAT(3,1,1) ! <dxy|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,2,1) ! <dyz|ypx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,2,1) ! <dzx|ypx>

            orbint(3)=orbint(3)+cta*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <dxx|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,2) ! <dyy|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,2) ! <dzz|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,2) ! <dxy|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,2) ! <dyz|zpx>
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,2) ! <dzx|zpx>

            orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,2)*IMAT(3,1,1) ! <dxx|xpy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,2)*IMAT(3,1,1) ! <dyy|xpy> 
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,3,1) ! <dzz|xpy> 
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,2)*IMAT(3,1,1) ! <dxy|xpy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,2,1) ! <dyz|xpy> 
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,2,1) ! <dzx|xpy> 

            orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,3)*IMAT(3,1,1) ! <dxx|ypy>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,3)*IMAT(3,1,1) ! <dyy|ypy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,3,1) ! <dzz|ypy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,3)*IMAT(3,1,1) ! <dxy|ypy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,2,1) ! <dyz|ypy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,2,1) ! <dzx|ypy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <dxx|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,2) ! <dyy|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,2) ! <dzz|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,2) ! <dxy|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,2) ! <dyz|zpy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,2) ! <dzx|zpy>

            orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,2) ! <dxx|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,2) ! <dyy|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,2) ! <dzz|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,2) ! <dxy|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,2) ! <dyz|xpz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,2) ! <dzx|xpz>

            orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,2) ! <dxx|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,2) ! <dyy|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,2) ! <dzz|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,2) ! <dxy|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,2) ! <dyz|ypz>
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,2) ! <dzx|ypz>

            orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,3) ! <dxx|zpz>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,3) ! <dyy|zpz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,3) ! <dzz|zpz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,3) ! <dxy|zpz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,3) ! <dyz|zpz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,3) ! <dzx|zpz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,1) !


       endif


! Case <  s | d' >
!--------------------------------------------------------------------!
       if ((ni.eq.0).and.(nj.eq.2)) then
            ct2=2*cij
            cta=ct2*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMAT)

            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,1,1)*IMAT(3,1,1) ! <s|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <s|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <s|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,2)*IMAT(3,1,1) ! <s|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <s|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <s|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,2) !

            orbint(2)=orbint(2)+cta*IMAT(1,1,3)*IMAT(2,1,2)*IMAT(3,1,1) ! <s|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,4)*IMAT(3,1,1) ! <s|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <s|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <s|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <s|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <s|ydzx>

            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <s|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <s|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,4) ! <s|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <s|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <s|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <s|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,1,1) !


       endif


! Case <  p | d' >
!--------------------------------------------------------------------!
       if ((ni.eq.1).and.(nj.eq.2)) then
            cij=coefs(ni,ii)*coefs(nj,jj)
            ct2=2*cij
            cta=ct2*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMAT)

            orbint(1)=orbint(1)+cta*IMAT(1,2,4)*IMAT(2,1,1)*IMAT(3,1,1) ! <px|xdxx>
     >                         -ct2*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,2,1)*IMAT(3,1,1) ! <py|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,1,1)*IMAT(3,2,1) ! <pz|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <px|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,3)*IMAT(3,1,1) ! <py|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,2,1) ! <pz|xdyy>

            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <px|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,3) ! <py|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,3) ! <pz|xdzz>

            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,2)*IMAT(3,1,1) ! <px|xdxy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,2)*IMAT(3,1,1) ! <py|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,2)*IMAT(3,2,1) ! <pz|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <px|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <py|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <pz|xdyz>

            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <px|xdzx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,1,2) ! <py|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,2,2) ! <pz|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,2) !


            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,4)*IMAT(3,1,1) ! <px|ydyy>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,4)*IMAT(3,1,1) ! <py|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,4)*IMAT(3,2,1) ! <pz|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <px|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,3) ! <py|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,3) ! <pz|ydzz>

            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <px|ydxy>
     >                         -cij*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,3)*IMAT(3,1,1) ! <py|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,2,1) ! <pz|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <px|ydyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,1,2) ! <py|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,2,2) ! <pz|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,2) !

            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <px|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <py|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <pz|ydzx>

            orbint(3)=orbint(3)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <px|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,1,2) ! <py|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,2,2) ! <pz|zdxx>

            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <px|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,1,2) ! <py|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,2,2) ! <pz|zdyy>

            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,4) ! <px|zdzz>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,4) ! <py|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,4) ! <pz|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,2,2) !

            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <px|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <py|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <pz|zdxy>

            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <px|zdyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,3) ! <py|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,3) ! <pz|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <px|zdzx>
     >                         -cij*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,3) ! <py|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,3) ! <pz|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,2,1) !

       endif


! Case <  d | d' >
!--------------------------------------------------------------------!
       if ((ni.eq.2).and.(nj.eq.2)) then
            ct2=2*cij
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMAT)

            orbint(1)=orbint(1)+cta*IMAT(1,3,4)*IMAT(2,1,1)*IMAT(3,1,1) ! <dxx|xdxx>
     >                         -ct2*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,3,1)*IMAT(3,1,1) ! <dyy|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,1,1)*IMAT(3,3,1) ! <dzz|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,4)*IMAT(2,2,1)*IMAT(3,1,1) ! <dxy|xdxx>
     >                         -ct2*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,4)*IMAT(2,2,1)*IMAT(3,2,1) ! <dyz|xdxx>
     >                         -ct2*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,4)*IMAT(2,1,1)*IMAT(3,2,1) ! <dzx|xdxx>
     >                         -ct2*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <dxx|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,3)*IMAT(3,1,1) ! <dyy|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,3,1) ! <dzz|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,3)*IMAT(3,1,1) ! <dxy|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,3)*IMAT(3,2,1) ! <dyz|xdyy>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,3)*IMAT(3,2,1) ! <dzx|xdyy>

            orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <dxx|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,3) ! <dyy|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,3) ! <dzz|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,3) ! <dxy|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,3) ! <dyz|xdzz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,3) ! <dzx|xdzz>

            orbint(1)=orbint(1)+cta*IMAT(1,3,3)*IMAT(2,1,2)*IMAT(3,1,1) ! <dxx|xdxy>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,3,2)*IMAT(3,1,1) ! <dyy|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,2)*IMAT(3,3,1) ! <dzz|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,2,2)*IMAT(3,1,1) ! <dxy|xdxy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,2)*IMAT(3,2,1) ! <dyz|xdxy>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,1) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,2)*IMAT(3,2,1) ! <dzx|xdxy>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(1)=orbint(1)+cta*IMAT(1,3,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <dxx|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,3,2)*IMAT(3,1,2) ! <dyy|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,3,2) ! <dzz|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <dxy|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,2,2) ! <dyz|xdyz>
            orbint(1)=orbint(1)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <dzx|xdyz>

            orbint(1)=orbint(1)+cta*IMAT(1,3,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <dxx|xdzx>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,3,1)*IMAT(3,1,2) ! <dyy|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,3,2) ! <dzz|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,2,1)*IMAT(3,1,2) ! <dxy|xdzx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,2,2) ! <dyz|xdzx>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,2) !
            orbint(1)=orbint(1)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,2,2) ! <dzx|xdzx>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,2) !


            orbint(2)=orbint(2)+cta*IMAT(1,3,3)*IMAT(2,1,2)*IMAT(3,1,1) ! <dxx|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,3)*IMAT(2,3,2)*IMAT(3,1,1) ! <dyy|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,3)*IMAT(2,1,2)*IMAT(3,3,1) ! <dzz|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,3)*IMAT(2,2,2)*IMAT(3,1,1) ! <dxy|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,3)*IMAT(2,2,2)*IMAT(3,2,1) ! <dyz|ydxx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,3)*IMAT(2,1,2)*IMAT(3,2,1) ! <dzx|ydxx>

            orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,4)*IMAT(3,1,1) ! <dxx|ydyy>
     >                         -ct2*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,4)*IMAT(3,1,1) ! <dyy|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,4)*IMAT(3,3,1) ! <dzz|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,4)*IMAT(3,1,1) ! <dxy|ydyy>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,4)*IMAT(3,2,1) ! <dyz|ydyy>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,4)*IMAT(3,2,1) ! <dzx|ydyy>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <dxx|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,3) ! <dyy|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,3) ! <dzz|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,3) ! <dxy|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,3) ! <dyz|ydzz>
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,3) ! <dzx|ydzz>

            orbint(2)=orbint(2)+cta*IMAT(1,3,2)*IMAT(2,1,3)*IMAT(3,1,1) ! <dxx|ydxy>
     >                         -cij*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,3,3)*IMAT(3,1,1) ! <dyy|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,3)*IMAT(3,3,1) ! <dzz|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,2,3)*IMAT(3,1,1) ! <dxy|ydxy>
     >                         -cij*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,3)*IMAT(3,2,1) ! <dyz|ydxy>
     >                         -cij*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,3)*IMAT(3,2,1) ! <dzx|ydxy>
     >                         -cij*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,1) !

            orbint(2)=orbint(2)+cta*IMAT(1,3,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <dxx|ydyz>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,3,3)*IMAT(3,1,2) ! <dyy|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,3,2) ! <dzz|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,2,3)*IMAT(3,1,2) ! <dxy|ydyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,2,2) ! <dyz|ydyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,2) !
            orbint(2)=orbint(2)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,2,2) ! <dzx|ydyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,2) !

            orbint(2)=orbint(2)+cta*IMAT(1,3,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <dxx|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,3,2)*IMAT(3,1,2) ! <dyy|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,3,2) ! <dzz|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <dxy|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,2,2) ! <dyz|ydzx>
            orbint(2)=orbint(2)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <dzx|ydzx>


            orbint(3)=orbint(3)+cta*IMAT(1,3,3)*IMAT(2,1,1)*IMAT(3,1,2) ! <dxx|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,3,1)*IMAT(3,1,2) ! <dyy|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,1,1)*IMAT(3,3,2) ! <dzz|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,2,3)*IMAT(2,2,1)*IMAT(3,1,2) ! <dxy|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,1,3)*IMAT(2,2,1)*IMAT(3,2,2) ! <dyz|zdxx>
            orbint(3)=orbint(3)+cta*IMAT(1,2,3)*IMAT(2,1,1)*IMAT(3,2,2) ! <dzx|zdxx>

            orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,3)*IMAT(3,1,2) ! <dxx|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,3)*IMAT(3,1,2) ! <dyy|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,3)*IMAT(3,3,2) ! <dzz|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,3)*IMAT(3,1,2) ! <dxy|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,3)*IMAT(3,2,2) ! <dyz|zdyy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,3)*IMAT(3,2,2) ! <dzx|zdyy>

            orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,4) ! <dxx|zdzz>
     >                         -ct2*IMAT(1,3,1)*IMAT(2,1,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,4) ! <dyy|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,3,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,4) ! <dzz|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,1,1)*IMAT(3,3,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,4) ! <dxy|zdzz>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,2,1)*IMAT(3,1,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,4) ! <dyz|zdzz>
     >                         -ct2*IMAT(1,1,1)*IMAT(2,2,1)*IMAT(3,2,2) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,4) ! <dzx|zdzz>
     >                         -ct2*IMAT(1,2,1)*IMAT(2,1,1)*IMAT(3,2,2) !

            orbint(3)=orbint(3)+cta*IMAT(1,3,2)*IMAT(2,1,2)*IMAT(3,1,2) ! <dxx|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,3,2)*IMAT(3,1,2) ! <dyy|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,2)*IMAT(3,3,2) ! <dzz|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,2,2)*IMAT(3,1,2) ! <dxy|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,2)*IMAT(3,2,2) ! <dyz|zdxy>
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,2)*IMAT(3,2,2) ! <dzx|zdxy>

            orbint(3)=orbint(3)+cta*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,3) ! <dxx|zdyz>
     >                         -cij*IMAT(1,3,1)*IMAT(2,1,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,3) ! <dyy|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,3,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,3) ! <dzz|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,1,2)*IMAT(3,3,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,3) ! <dxy|zdyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,2,2)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,3) ! <dyz|zdyz>
     >                         -cij*IMAT(1,1,1)*IMAT(2,2,2)*IMAT(3,2,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,3) ! <dzx|zdyz>
     >                         -cij*IMAT(1,2,1)*IMAT(2,1,2)*IMAT(3,2,1) !

            orbint(3)=orbint(3)+cta*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,3) ! <dxx|zdzx>
     >                         -cij*IMAT(1,3,2)*IMAT(2,1,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,3) ! <dyy|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,3,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,3) ! <dzz|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,1,1)*IMAT(3,3,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,3) ! <dxy|zdzx>
     >                         -cij*IMAT(1,2,2)*IMAT(2,2,1)*IMAT(3,1,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,3) ! <dyz|zdzx>
     >                         -cij*IMAT(1,1,2)*IMAT(2,2,1)*IMAT(3,2,1) !
            orbint(3)=orbint(3)+cta*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,3) ! <dzx|zdzx>
     >                         -cij*IMAT(1,2,2)*IMAT(2,1,1)*IMAT(3,2,1) !


       endif


       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
