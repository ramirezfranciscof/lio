!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_finalize()
!--------------------------------------------------------------------!
! DEALLOCATION OF GLOBAL VARIABLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       implicit none

       integer igpu
!--------------------------------------------------------------------!
       if (allocated(Smat))    deallocate(Smat)
       if (allocated(RealRho)) deallocate(RealRho)
!--------------------------------------------------------------------!
       deallocate(r,v,rqm, Em, Rm, pc,Iz, nnat,
     > af,c,a,cx,ax,cd,ad,B,Nuc,ncont,Nucx,ncontx,Nucd
     > ,ncontd, indexii, indexiid, RMM, X, XX)
c       deallocate(old1,old2,old3)
       deallocate(natomc,nnps,nnpp,nnpd,nns)
       deallocate(nnd,nnp,atmin,jatc,d)
       call g2g_deinit()
       call aint_query_gpu_level(igpu)
       if (igpu.gt.1) then
         call aint_deinit()
       endif
!--------------------------------------------------------------------!
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
