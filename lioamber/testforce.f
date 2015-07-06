!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine testforce(Sinv,Fmtx,Pmtx)
!--------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod, only:M,a,c,r,nshell,natom,ncont,nl,nuc,RMM,Md
       use intmod
       use testmod
       implicit none
       real*8,intent(in)      :: Sinv(M,M)
       real*8,intent(in)      :: Fmtx(M,M)
       complex*16,intent(in)  :: Pmtx(M,M)

       real*8,allocatable     :: nr(:,:),nv(:,:),as(:,:),cs(:,:)
       real*8,allocatable     :: ffold(:,:)
       complex*16,allocatable :: ffnew(:,:)
       complex*16,allocatable :: AuxMat(:,:), Mdir(:,:),Mtrp(:,:)
       complex*16,allocatable :: Bmat(:,:),Mato(:,:)
       integer,allocatable    :: nucof(:)
       integer                :: kk,ii,jj,unitid,nk,nn
       integer                :: MM,MMd,indx

!--------------------------------------------------------------------!
       allocate(nr(3,natom),nv(3,natom),as(nl,M),cs(nl,M),nucof(M))
       allocate(AuxMat(M,M),Mdir(M,M),Mtrp(M,M),Bmat(M,M))
       allocate(ffold(natom,3),ffnew(3,natom))
       allocate(Mato(M,M))


       call g2g_timer_start('intold_')
       ffold=0.0d0
       call testinit()
       DSX(:,:,:)=0.0d0
       DSY(:,:,:)=0.0d0
       DSZ(:,:,:)=0.0d0
!       call calcDSM(ffold)
       call intSG(ffold)
       call g2g_timer_stop('intold_')
       do kk=1,natom
        do ii=1,M
        do jj=1,M
          DSXc(ii,jj,kk)=DSX(ii,jj,kk)
          DSYc(ii,jj,kk)=DSY(ii,jj,kk)
          DSZc(ii,jj,kk)=DSZ(ii,jj,kk)
        enddo
        enddo
       enddo


!      Crear Mdir y Mtrp
       Mdir=matmul(Fmtx,Pmtx)
       Mdir=matmul(Sinv,Mdir)
       Mato=Mdir
       Mtrp=matmul(Pmtx,Fmtx)
       Mtrp=matmul(Mtrp,Sinv)
       Mtrp=transpose(Mtrp)
       AuxMat=Mdir+Mtrp

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

!      Calcula las fuerzas
       call g2g_timer_start('intnew_')
       DSX(:,:,:)=0.0d0
       DSY(:,:,:)=0.0d0
       DSZ(:,:,:)=0.0d0
       ffnew=CMPLX(0.0d0,0.0d0)
       call fzaDS2(natom,M,nshell(0),nshell(1),ncont,nl,
     >             AuxMat,nr,nv,as,cs,nucof,Bmat,ffnew)
       call g2g_timer_stop('intnew_')
       do kk=1,natom
        do ii=1,M
        do jj=1,M
          DSXt(ii,jj,kk)=DSX(jj,ii,kk)
          DSYt(ii,jj,kk)=DSY(jj,ii,kk)
          DSZt(ii,jj,kk)=DSZ(jj,ii,kk)
        enddo
        enddo
       enddo

       DSX=DSX+DSXt
       DSY=DSY+DSYt
       DSZ=DSZ+DSZt
       unitid=300
       do kk=1,natom
        do ii=1,M
        do jj=1,M
          write(unitid+10+kk,*) ii,jj,DSXc(ii,jj,kk),DSX(ii,jj,kk)
          write(unitid+20+kk,*) ii,jj,DSYc(ii,jj,kk),DSY(ii,jj,kk)
          write(unitid+30+kk,*) ii,jj,DSZc(ii,jj,kk),DSZ(ii,jj,kk)
        enddo
        enddo
       enddo

       do ii=1,M
       do jj=1,M
         write(300,*) ii,jj
         write(300,*) r(nuc(ii),1),r(nuc(ii),2),r(nuc(ii),3)
         do kk=1,ncont(ii)
           write(300,*) c(ii,kk),a(ii,kk)
         enddo
         write(300,*)
         write(300,*) r(nuc(jj),1),r(nuc(jj),2),r(nuc(jj),3)
         do kk=1,ncont(jj)
           write(300,*) c(jj,kk),a(jj,kk)
         enddo
         write(300,*)
         write(300,*) DSX(ii,jj,1),DSX(ii,jj,2),DSX(ii,jj,3)
         write(300,*) DSXc(ii,jj,1),DSXc(ii,jj,2),DSXc(ii,jj,3)
         write(300,*)
         write(300,*) DSY(ii,jj,1),DSY(ii,jj,2),DSY(ii,jj,3)
         write(300,*) DSYc(ii,jj,1),DSYc(ii,jj,2),DSYc(ii,jj,3)
         write(300,*) 
         write(300,*) DSZ(ii,jj,1),DSZ(ii,jj,2),DSZ(ii,jj,3)
         write(300,*) DSZc(ii,jj,1),DSZc(ii,jj,2),DSZc(ii,jj,3)
         write(300,*) 
         write(300,*)
         write(300,*)
       enddo
       enddo



       print*,'     natom         dir'
       do ii=1,natom
       do kk=1,3
         print*,ii,kk,ffold(ii,kk),ffnew(kk,ii)
       enddo
       enddo
       MM=M*(M+1)/2
       MMd=Md*(Md+1)/2
       indx=1+M+2*MMd+4*MM
!       call sprepack('L',M,RMM(indx),Mato)
!       call messrho(M,RMM(indx))
       ffold(:,:)=0.0
       call intSG(ffold)
!       call calcDSM(ffold)
       print*,'     natom         dir'
       do ii=1,natom
       do kk=1,3
         print*,ii,kk,ffold(ii,kk),ffnew(kk,ii)
       enddo
       enddo


!--------------------------------------------------------------------!
       deallocate(Mato)
       deallocate(nr,nv,as,cs,nucof)
       deallocate(AuxMat,Mdir,Mtrp,Bmat)
       deallocate(ffold,ffnew)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
