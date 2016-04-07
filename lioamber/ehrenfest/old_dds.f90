!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine fzaDDS(Natoms,Nbasis,nbs,nbp,Nconts,Ncontm,AuxMat, &
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
  integer    :: ii,ni,nki,jj,nj,nkj,ki,kj,kk
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
     orbpot(1,kk+0)=2   ! dxx (x)
     orbpot(1,kk+1)=1   ! dxy (x)
     orbpot(2,kk+1)=1   ! dxy (y)
     orbpot(2,kk+2)=2   ! dyy (y)
     orbpot(1,kk+3)=1   ! dxz (x)
     orbpot(3,kk+3)=1   ! dxz (z)
     orbpot(2,kk+4)=1   ! dyz (y)
     orbpot(3,kk+4)=1   ! dyz (z)
     orbpot(3,kk+5)=2   ! dzz (z)

     alpha(:,kk+1)=alpha(:,kk)
     alpha(:,kk+2)=alpha(:,kk)
     alpha(:,kk+3)=alpha(:,kk)
     alpha(:,kk+4)=alpha(:,kk)
     alpha(:,kk+5)=alpha(:,kk)

     coefs(:,kk+1)=coefs(:,kk)             ! dxy
     coefs(:,kk+2)=coefs(:,kk)/SQRT(3.0)   ! dyy
     coefs(:,kk+3)=coefs(:,kk)             ! dxz
     coefs(:,kk+4)=coefs(:,kk)             ! dyz
     coefs(:,kk+5)=coefs(:,kk)/SQRT(3.0)   ! dzz
     coefs(:,kk+0)=coefs(:,kk)/SQRT(3.0)   ! dxx
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

      do kj=1,3
      do ki=1,3
         pi(:)=opi(:)
         pj(:)=opj(:)

         term1=0.0d0
         pi(ki)=opi(ki)+1
         pj(kj)=opj(kj)+1
         intx=IMTX(1,1+pi(1),1+pj(1))
         inty=IMTX(2,1+pi(2),1+pj(2))
         intz=IMTX(3,1+pi(3),1+pj(3))
         term1=cij*ani*anj*intx*inty*intz

         term2=0.0d0
         if (opj(kj).gt.0) then
           pi(ki)=opi(ki)+1
           pj(kj)=opj(kj)-1
           intx=IMTX(1,1+pi(1),1+pj(1))
           inty=IMTX(2,1+pi(2),1+pj(2))
           intz=IMTX(3,1+pi(3),1+pj(3))
           cj2=-1.0d0
           if (pj(kj).eq.2) cj2=-2.0d0
           term2=cij*ani*cj2*intx*inty*intz
         endif

         term3=0.0d0
         if (opi(ki).gt.0) then
           pi(ki)=opi(ki)-1
           pj(kj)=opj(kj)+1
           intx=IMTX(1,1+pi(1),1+pj(1))
           inty=IMTX(2,1+pi(2),1+pj(2))
           intz=IMTX(3,1+pi(3),1+pj(3))
           ci2=-1.0d0
           if (pi(ki).eq.2) ci2=-2.0d0
           term3=cij*ci2*anj*intx*inty*intz
         endif

         term4=0.0d0
         if ((opi(ki).gt.0).and.(opj(kj).gt.0)) then
           pi(ki)=opi(ki)-1
           pj(kj)=opj(kj)-1
           intx=IMTX(1,1+pi(1),1+pj(1))
           inty=IMTX(2,1+pi(2),1+pj(2))
           intz=IMTX(3,1+pi(3),1+pj(3))
           term4=cij*ci2*cj2*intx*inty*intz
         endif

         term=term1+term2+term3+term4
         Force(kj,nkj)=Force(kj,nkj)+AuxMat(ii,jj)*NucVel(ki,nki)*term
      enddo
      enddo

    enddo
    enddo

  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
