!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine fzaDS2(Natoms,Nbasis,nbs,nbp,Nconts,Ncontm,AuxMat, &
               nucpos,nucvel,alpha,coefs,nucof,Bmat,force)
!--------------------------------------------------------------------!
!
! Ya fueron probados los < s / ds > y < p / ds >
!
!
!--------------------------------------------------------------------!
  use testmod
  implicit none
  integer,intent(in)     :: Natoms          ! Number of atoms
  integer,intent(in)     :: Nbasis          ! Number of basis
  integer,intent(in)     :: nbs,nbp         !
  integer,intent(in)     :: Nconts(Nbasis)  ! Number of contractionbs
  integer,intent(in)     :: Ncontm          ! Max contractionbs

  complex*16,intent(in)  :: AuxMat(Nbasis,Nbasis)
  real*8,intent(in)      :: nucpos(3,Natoms),nucvel(3,Natoms)
  real*8,intent(in)      :: alpha(Ncontm,Nbasis)
  real*8,intent(in)      :: coefs(Ncontm,Nbasis)
  integer,intent(in)     :: nucof(Nbasis)

  complex*16,intent(out) :: Bmat(Nbasis,Nbasis)
  complex*16,intent(out) :: force(3,Natoms)

  real*8     :: orbint(3)
  real*8     :: IMTX(3,4,4)
  real*8     :: posj(3),posi(3),anj,ani
  real*8     :: cij,cta,ct2
  integer    :: ii,ni,nki,jj,nj,nkj,kk
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  force=CMPLX(0.0d0,0.0d0)
  Bmat=CMPLX(0.0d0,0.0d0)
  DSX(:,:,:)=0.0d0
  DSY(:,:,:)=0.0d0
  DSZ(:,:,:)=0.0d0


! Integrals of <  s | s' >
!--------------------------------------------------------------------!
  do jj=1,nbs
  do ii=1,nbs
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj

      call setim(0,1,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xs>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ys>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zs>

    enddo
    enddo

    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo


! Integrals of <  p | s' >
!--------------------------------------------------------------------!
  do jj=1,nbs
  do ii=nbs+1,nbs+nbp,3
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)

      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(1,1,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xs>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ys>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zs>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo


! Integrals of <  d | s' >
!--------------------------------------------------------------------!
  do jj=1,nbs
  do ii=nbs+nbp+1,Nbasis,6
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)

    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(2,1,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xs>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xs>

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ys>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ys>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zs>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zs>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo




! Integrals of <  s | p' >
!--------------------------------------------------------------------!
  do jj=nbs+1,nbs+nbp,3
  do ii=1,nbs
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(0,2,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xpy> 
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xpz>

      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ypz>

      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|zpz>

!     write(500,*) ii,ni,jj,nj
!     write(500,*) coefs(ni,ii),coefs(nj,jj)
!     write(500,*) ani,anj
!     write(500,*) posi
!     write(500,*) posj
!     write(500,*) orbint
!     write(500,*) 

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo


! Integrals of <  p | p' >
!--------------------------------------------------------------------!
  do jj=nbs+1,nbs+nbp,3
  do ii=nbs+1,nbs+nbp,3
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(1,2,ani,anj,posi,posj,IMTX)


      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xpx>

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xpy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xpy> 
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xpy>

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xpz>


      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ypx>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|ypy>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ypz>


      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zpx>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zpy>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|zpz>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo



! Integrals of <  d | p' >
!--------------------------------------------------------------------!
  do jj=nbs+1,nbs+nbp,3
  do ii=nbs+nbp+1,Nbasis,6
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(2,2,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xpx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xpx>

      orbint(2)=orbint(2)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ypx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ypx>

      orbint(3)=orbint(3)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zpx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zpx>


      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|xpy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|xpy> 
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|xpy> 
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|xpy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|xpy> 
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|xpy> 

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|ypy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|ypy>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|zpy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|zpy>

      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|xpz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|xpz>

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|ypz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|ypz>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|zpz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|zpz>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo



! Integrals of <  s | d' >
!--------------------------------------------------------------------!
  do jj=nbs+nbp+1,Nbasis,6
  do ii=1,nbs
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(0,3,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xdzx>

      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <s|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ydzx>

      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <s|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|zdzx>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo



! Integrals of <  p | d' >
!--------------------------------------------------------------------!
  do jj=nbs+nbp+1,Nbasis,6
  do ii=nbs+1,nbs+nbp,3
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(1,3,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,2,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xdxx>

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|xdyy>

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|xdzz>

      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xdxy>

      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|xdyz>

      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xdzx>


      orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ydxx>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,2,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ydyy>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <px|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <py|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <pz|ydzz>

      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|ydxy>

      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|ydyz>

      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ydzx>


      orbint(3)=orbint(3)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zdxx>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <px|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <py|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <pz|zdyy>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zdzz>

      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zdxy>

      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|zdyz>

      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|zdzx>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo



! Integrals of <  d | d' >
!--------------------------------------------------------------------!
  do jj=nbs+nbp+1,Nbasis,6
  do ii=nbs+nbp+1,Nbasis,6
    nkj=nucof(jj)
    nki=nucof(ii)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)


    orbint(:)=0.0d0
    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      cij=coefs(ni,ii)*coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj
      call setim(2,3,ani,anj,posi,posj,IMTX)

      orbint(1)=orbint(1)+cta*IMTX(1,3,4)*IMTX(2,1,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,3,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,3,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,2,1)*IMTX(3,1,1) &
                         -ct2*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,2,1) &
                         -ct2*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xdxx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,2,1) &
                         -ct2*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xdxx>

      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <dxx|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,3)*IMTX(3,1,1) ! <dyy|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,3,1) ! <dzz|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <dxy|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,2,1) ! <dyz|xdyy>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <dzx|xdyy>

      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|xdzz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|xdzz>

      orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|xdxy>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|xdxy>

      orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|xdyz>
      orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|xdyz>

      orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|xdzx>
      orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|xdzx>


      orbint(2)=orbint(2)+cta*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ydxx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ydxx>

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,3,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,4)*IMTX(3,1,1) &
                         -ct2*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,2,1) &
                         -ct2*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ydyy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,2,1) &
                         -ct2*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ydyy>

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <dxx|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) ! <dyy|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) ! <dzz|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <dxy|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) ! <dyz|ydzz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <dzx|ydzz>

      orbint(2)=orbint(2)+cta*IMTX(1,3,2)*IMTX(2,1,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,3,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,3,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,2,3)*IMTX(3,1,1) &
                         -cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|ydxy>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,2,1) &
                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|ydxy>

      orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,2) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|ydyz>
      orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,2) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|ydyz>

      orbint(2)=orbint(2)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|ydzx>
      orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|ydzx>


      orbint(3)=orbint(3)+cta*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zdxx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zdxx>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <dxx|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,2) ! <dyy|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,2) ! <dzz|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <dxy|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,2) ! <dyz|zdyy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <dzx|zdyy>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,4) &
                         -ct2*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,4) &
                         -ct2*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zdzz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,4) &
                         -ct2*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zdzz>

      orbint(3)=orbint(3)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|zdxy>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|zdxy>

      orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|zdyz>
      orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) &
                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|zdyz>

      orbint(3)=orbint(3)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,3) &
                         -cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|zdzx>
      orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,3) &
                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|zdzx>

    enddo
    enddo
    do kk=1,3
      force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
      Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
    enddo

    DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+orbint(1)
    DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+orbint(2)
    DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+orbint(3)
  enddo
  enddo


  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
