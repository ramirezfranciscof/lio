!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine fzaDS2(Natoms,Nbasis,nbs,nbp,Nconts,Ncontm,
     >                   AuxMat,nucpos,nucvel,alpha,coefs,nucof,
     >                   Bmat,force)
!--------------------------------------------------------------------!
!
!
!
!
!--------------------------------------------------------------------!
       implicit none
       integer,intent(in)     :: Natoms  ! Number of atoms
       integer,intent(in)     :: Nbasis  ! Number of basis
       integer,intent(in)     :: nbs,nbp ! 
       integer,intent(in)     :: Nconts(Nbasis)
                                         ! Number of contractionbs
       integer,intent(in)     :: Ncontm  ! Max contractionbs

       complex*16,intent(in)  :: AuxMat(Nbasis,Nbasis)
       real*8,intent(in)      :: nucpos(3,Natoms),nucvel(3,Natoms)
       real*8,intent(in)      :: alpha(Ncontm,Nbasis)
       real*8,intent(in)      :: coefs(Ncontm,Nbasis)
       integer,intent(in)     :: nucof(Nbasis)

       complex*16,intent(out) :: Bmat(Nbasis,Nbasis)
       complex*16,intent(out) :: force(3,Natoms)

       complex*16 :: orbint(3)
       real*8     :: IMTX(3,4,4)
       real*8     :: posj(3),posi(3),anj,ani
       real*8     :: cij,cta,ct2
       integer    :: ii,ni,nki,jj,nj,nkj,kk
       logical    :: samenuc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       force=CMPLX(0.0d0,0.0d0)
       Bmat=CMPLX(0.0d0,0.0d0)

! Derivatives of s
!--------------------------------------------------------------------!
       do jj=1,nbs
       nkj=nucof(jj)
       posj(:)=nucpos(:,nkj)
       do nj=1,Nconts(jj)
         anj=alpha(nj,jj)


         do ii=1,nbs
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xs>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ys>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zs>

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         if (.false.) then !PUENTEO

         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

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
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

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
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo

         endif !PUENTEO

       enddo
       enddo


! Derivatives of p
!--------------------------------------------------------------------!
       if (.false.) then !PUENTEO

       do jj=nbs+1,nbs+nbp,3
       nkj=nucof(jj)
       posj(:)=nucpos(:,nkj)
       do nj=1,Nconts(jj)
         anj=alpha(nj,jj)


         do ii=1,nbs
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xpy> 
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xpz>

            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ypx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ypz>

            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zpx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zpy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)


            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xpx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xpy>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xpy> 
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xpy>

            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xpz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xpz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xpz>


            orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ypx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ypx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ypx>

            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|ypy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ypz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ypz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ypz>


            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zpx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zpx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zpx>

            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zpy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zpy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zpy>

            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|zpz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xpx>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xpx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xpx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xpx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

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

            orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,1) ! <dxx|ypy>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,1) ! <dyy|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,1) ! <dzz|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,1) ! <dxy|ypy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,1) ! <dyz|ypy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,1) ! <dzx|ypy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

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

            orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|zpz>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|zpz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|zpz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|zpz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


       enddo
       enddo


! Derivatives of d
!--------------------------------------------------------------------!
       do jj=nbs+nbp+1,Nbasis,6
       nkj=nucof(jj)
       posj(:)=nucpos(:,nkj)
       do nj=1,Nconts(jj)
         anj=alpha(nj,jj)


         do ii=1,nbs
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            ct2=2*cij
            cta=ct2*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|xdyy>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|xdzz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !

            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <s|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <s|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <s|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ydzx>

            orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zdxx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <s|zdyy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <s|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <s|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            ct2=2*cij
            cta=ct2*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xdxx>
     >                         -ct2*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|xdyy>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|xdyy>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|xdyy>

            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|xdzz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|xdzz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|xdzz>

            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xdxy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|xdyz>

            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xdzx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !


            orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ydxx>

            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <px|ydyy>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,1,1) ! <py|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,2,1) ! <pz|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <px|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <py|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <pz|ydzz>

            orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|ydxy>
     >                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <px|ydyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <py|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <pz|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ydzx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ydzx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ydzx>


            orbint(3)=orbint(3)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zdxx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zdxx>
            orbint(3)=orbint(3)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zdxx>

            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <px|zdyy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <py|zdyy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <pz|zdyy>

            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <px|zdzz>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,4) ! <py|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,4) ! <pz|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zdxy>

            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <px|zdyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <py|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <pz|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|zdzx>
     >                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo


         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          posi(:)=nucpos(:,nki)
          samenuc=.false.
          if (nkj.eq.nki) samenuc=.true.
          orbint(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)
            ani=alpha(ni,ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            ct2=2*cij
            cta=2*cij*anj
            call setim(samenuc,0,1,ani,anj,posi,posj,IMTX)

            orbint(1)=orbint(1)+cta*IMTX(1,3,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xdxx>
     >                         -ct2*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xdxx>
     >                         -ct2*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xdxx>
     >                         -ct2*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xdxx>
     >                         -ct2*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

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

            orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|xdxy>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|xdxy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|xdxy>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|xdxy>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(1)=orbint(1)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|xdyz>
            orbint(1)=orbint(1)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|xdyz>

            orbint(1)=orbint(1)+cta*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|xdzx>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|xdzx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|xdzx>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            orbint(1)=orbint(1)+cta*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|xdzx>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !


            orbint(2)=orbint(2)+cta*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ydxx>
            orbint(2)=orbint(2)+cta*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ydxx>

            orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <dxx|ydyy>
     >                         -ct2*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,4)*IMTX(3,1,1) ! <dyy|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,3,1) ! <dzz|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,4)*IMTX(3,1,1) ! <dxy|ydyy>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,2,1) ! <dyz|ydyy>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,2,1) ! <dzx|ydyy>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <dxx|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) ! <dyy|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) ! <dzz|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <dxy|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) ! <dyz|ydzz>
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <dzx|ydzz>

            orbint(2)=orbint(2)+cta*IMTX(1,3,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <dxx|ydxy>
     >                         -cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,3,3)*IMTX(3,1,1) ! <dyy|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,3,1) ! <dzz|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <dxy|ydxy>
     >                         -cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,2,1) ! <dyz|ydxy>
     >                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <dzx|ydxy>
     >                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            orbint(2)=orbint(2)+cta*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <dxx|ydyz>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,2) ! <dyy|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,2) ! <dzz|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <dxy|ydyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,2) ! <dyz|ydyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            orbint(2)=orbint(2)+cta*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <dzx|ydyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !

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

            orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <dxx|zdzz>
     >                         -ct2*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,4) ! <dyy|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,4) ! <dzz|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,4) ! <dxy|zdzz>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,4) ! <dyz|zdzz>
     >                         -ct2*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,4) ! <dzx|zdzz>
     >                         -ct2*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            orbint(3)=orbint(3)+cta*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|zdxy>
            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|zdxy>

            orbint(3)=orbint(3)+cta*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <dxx|zdyz>
     >                         -cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) ! <dyy|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) ! <dzz|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <dxy|zdyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) ! <dyz|zdyz>
     >                         -cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <dzx|zdyz>
     >                         -cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            orbint(3)=orbint(3)+cta*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|zdzx>
     >                         -cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|zdzx>
     >                         -cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|zdzx>
     >                         -cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            orbint(3)=orbint(3)+cta*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|zdzx>
     >                         -cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          do kk=1,3
           force(kk,nkj)=force(kk,nkj)+AuxMat(ii,jj)*orbint(kk)
           Bmat(ii,jj)=Bmat(ii,jj)+nucvel(kk,nkj)*orbint(kk)
          enddo
         enddo

       enddo
       enddo

       endif !PUENTEO


       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
