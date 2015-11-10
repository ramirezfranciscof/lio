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
  real*8,intent(inout)   :: alpha(Ncontm,Nbasis)
  real*8,intent(inout)   :: coefs(Ncontm,Nbasis)
  integer,intent(in)     :: nucof(Nbasis)

  complex*16,intent(out) :: Bmat(Nbasis,Nbasis)
  complex*16,intent(out) :: force(3,Natoms)

  real*8,allocatable     :: orbpot(:,:)

  real*8     :: IMTX(3,4,4)
  real*8     :: posj(3),posi(3),anj,ani
  real*8     :: cij,cta,ct2
  real*8     :: intx,inty,intz,term1,term2
  integer    :: opi(3),opj(3),pi(3),pj(3)
  integer    :: ii,ni,nki,jj,nj,nkj,kk
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  force=DCMPLX(0.0d0,0.0d0)
  Bmat=DCMPLX(0.0d0,0.0d0)
  DSX(:,:,:)=0.0d0
  DSY(:,:,:)=0.0d0
  DSZ(:,:,:)=0.0d0

  allocate(orbpot(3,Nbasis))

! Prepare variables
!--------------------------------------------------------------------!
  orbpot(:,:)=0

  do kk=nbs+1,nbs+nbp,3
     orbpot(1,kk+0)=1  ! px
     orbpot(2,kk+1)=1  ! py
     orbpot(3,kk+2)=1  ! pz

     alpha(:,kk+1)=alpha(:,kk)
     alpha(:,kk+2)=alpha(:,kk)
     coefs(:,kk+1)=coefs(:,kk)
     coefs(:,kk+2)=coefs(:,kk)
  enddo

  do kk=nbs+nbp+1,Nbasis,6
     orbpot(1,kk+0)=2  ! dxx (x)
     orbpot(1,kk+1)=1  ! dxy (x)
     orbpot(2,kk+1)=1  ! dxy (y)
     orbpot(2,kk+2)=2  ! dyy (y)
     orbpot(1,kk+3)=1  ! dxz (x)
     orbpot(3,kk+3)=1  ! dxz (z)
     orbpot(2,kk+4)=1  ! dyz (y)
     orbpot(3,kk+4)=1  ! dyz (z)
     orbpot(3,kk+5)=2  ! dzz (z)

     alpha(:,kk+1)=alpha(:,kk)
     alpha(:,kk+2)=alpha(:,kk)
     alpha(:,kk+3)=alpha(:,kk)
     alpha(:,kk+4)=alpha(:,kk)
     alpha(:,kk+5)=alpha(:,kk)

     coefs(:,kk+1)=coefs(:,kk)
     coefs(:,kk+2)=coefs(:,kk)/SQRT(3.0)
     coefs(:,kk+3)=coefs(:,kk)
     coefs(:,kk+4)=coefs(:,kk)
     coefs(:,kk+5)=coefs(:,kk)/SQRT(3.0)
     coefs(:,kk+0)=coefs(:,kk)/SQRT(3.0)
  enddo


! Calculate Integrals
!--------------------------------------------------------------------!
  do jj=1,Nbasis
  do ii=1,Nbasis
    nkj=nucof(jj)
    nki=nucof(ii)
    opi(:)=orbpot(:,ii)
    opj(:)=orbpot(:,jj)
    posj(:)=nucpos(:,nkj)
    posi(:)=nucpos(:,nki)

    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      call setim(0,1,ani,anj,posi,posj,IMTX)

      cij=coefs(ni,ii)*coefs(nj,jj)
      write(555,*) ii,coefs(ni,ii),' .... ',jj,coefs(nj,jj)
      ct2=cij*2
      cta=ct2*anj

      do kk=1,3
         pi(:)=opi(:)
         pj(:)=opj(:)
         pj(kk)=pj(kk)+1
         intx=IMTX(1,1+pi(1),1+pj(1))
         inty=IMTX(2,1+pi(2),1+pj(2))
         intz=IMTX(3,1+pi(3),1+pj(3))
         term1=cta*intx*inty*intz

         term2=0.0d0
         if (pj(kk).ge.2) then
           pj(kk)=pj(kk)-2
           intx=IMTX(1,1+pi(1),1+pj(1))
           inty=IMTX(2,1+pi(2),1+pj(2))
           intz=IMTX(3,1+pi(3),1+pj(3))
           if (pj(kk).eq.0) term2=-cij*intx*inty*intz
           if (pj(kk).eq.1) term2=-ct2*intx*inty*intz
         endif

         Bmat(ii,jj)=Bmat(ii,jj)+NucVel(kk,nkj)*(term1+term2)
         Force(kk,nkj)=Force(kk,nkj)+AuxMat(ii,jj)*(term1+term2)

         if (kk.eq.1) DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+(term1+term2)
         if (kk.eq.2) DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+(term1+term2)
         if (kk.eq.3) DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+(term1+term2)

!      enddo
      if (nki.eq.1) then
      if (((ii.eq.1).and.(jj.eq.13)).or.((ii.eq.13).and.(jj.eq.1))) then
        write(666,*) ii,jj
        write(666,*) cta*intx*inty*intz,-cij*intx*inty*intz,-ct2*intx*inty*intz
        write(666,*) term1,term2
        write(666,*) DSX(ii,jj,nkj),DSY(ii,jj,nkj),DSZ(ii,jj,nkj)
        write(666,*)
        write(666,*)
      endif
      endif

      enddo

    enddo
    enddo

  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
