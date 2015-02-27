!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine testforce(Sinv,Fmtx,Pmtx)
!--------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod, only:M,a,c,r,nshell,natom,ncont,nl,nuc
       use intmod
       implicit none
       real*8,intent(in)      :: Sinv(M,M)
       real*8,intent(in)      :: Fmtx(M,M)
       complex*16,intent(in)  :: Pmtx(M,M)

       real*8,allocatable     :: nr(:,:),nv(:,:),as(:,:),cs(:,:)
       real*8,allocatable     :: ffold(:,:)
       complex*16,allocatable :: ffnew(:,:)
       complex*16,allocatable :: MatT(:,:),MatD(:,:),Bmat(:,:)
       integer,allocatable    :: nucof(:)
       integer                :: kk,ii,jj

!--------------------------------------------------------------------!
       allocate(nr(3,natom),nv(3,natom),as(nl,M),cs(nl,M),nucof(M))
       allocate(MatT(M,M),MatD(M,M),Bmat(M,M))
       allocate(ffold(natom,3),ffnew(3,natom))

       call g2g_timer_start('intold_')
       ffold=0.0d0
       call intSG(ffold)
       call g2g_timer_stop('intold_')

!      Crear MatT y MatD
       MatD=matmul(Fmtx,Pmtx)
       MatD=matmul(Sinv,MatD)
       MatT=transpose(MatD)

!      Transpone a y c
       do ii=1,M
         nucof(ii)=nuc(ii)
         do jj=1,nl
           as(jj,ii)=a(ii,jj)
           cs(jj,ii)=c(ii,jj)
         enddo
       enddo

!      Transpone r y setea nv=0
       do ii=1,natom
       do jj=1,3
         nr(jj,ii)=r(ii,jj)
         nv(jj,ii)=0.0d0
       enddo
       enddo
       call g2g_timer_start('intnew_')
       ffnew=CMPLX(0.0d0,0.0d0)
       print*,'About to get in'
       call fzaDS2(natom,M,nshell(0),nshell(1),ncont,nl,
     >             MatT,MatD,nr,nv,as,cs,nucof,Bmat,ffnew)
       call g2g_timer_stop('intnew_')

       do ii=1,natom
       do kk=1,3
         print*,ffold(ii,kk),ffnew(kk,ii)
       enddo
       enddo

!--------------------------------------------------------------------!
       deallocate(nr,nv,as,cs,nucof)
       deallocate(MatT,MatD,Bmat)
       deallocate(ffold,ffnew)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
