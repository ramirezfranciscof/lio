!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine forceDS_dss(natoms,nbasis,pos,vel,Matin,MatB,fterm)
!--------------------------------------------------------------------!
!
!
!--------------------------------------------------------------------!
  use testmod
  use basis_copy
  implicit none
  integer,intent(in)    :: natoms          ! Number of atoms
  integer,intent(in)    :: nbasis          ! Number of basis
  real*8,intent(in)     :: pos(3,natoms)
  real*8,intent(in)     :: vel(3,natoms)
  complex*16,intent(in) :: Matin(nbasis,nbasis)
  real*8,intent(out)    :: MatB(nbasis,nbasis)
  complex*8,intent(out) :: ftermSinv(3,natoms)

  real*8     :: IMTX(3,4,4)
  real*8     :: posj(3),posi(3),anj,ani
  real*8     :: cij,cta,ct2
  real*8     :: intx,inty,intz,term1,term2
  integer    :: opi(3),opj(3),pi(3),pj(3)
  integer    :: ii,ni,nki,jj,nj,nkj,kk
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  fterm=dcmplx(0.0d0,0.0d0)

  do jj=1,Nbasis
  do ii=1,Nbasis
    nkj=patent_atom(jj)
    nki=parent_atom(ii)
    opi(:)=orbpot(:,ii)
    opj(:)=orbpot(:,jj)
    posj(:)=pos(:,nkj)
    posi(:)=pos(:,nki)

    do nj=1,Nconts(jj)
    do ni=1,Nconts(ii)
      anj=alpha(nj,jj)
      ani=alpha(ni,ii)
      call setim(0,1,ani,anj,posi,posj,IMTX)

      cij=coefs(ni,ii)*coefs(nj,jj)
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
         forceDS(kk,nkj)=forceDS(kk,nkj)+AuxMat(ii,jj)*(term1+term2)

         if (kk.eq.1) DSX(ii,jj,nkj)=DSX(ii,jj,nkj)+(term1+term2)
         if (kk.eq.2) DSY(ii,jj,nkj)=DSY(ii,jj,nkj)+(term1+term2)
         if (kk.eq.3) DSZ(ii,jj,nkj)=DSZ(ii,jj,nkj)+(term1+term2)
      enddo

    enddo
    enddo

  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
