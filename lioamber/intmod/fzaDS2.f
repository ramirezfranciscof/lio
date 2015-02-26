!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine fzaDS2 (Natoms,Nbasis,nbs,nbp,Nconts,Ncontm,
     >                    MatT,MatD,nucpos,nucvel,alpha,coefs,nucof,
     >                    Bmat,force)
!--------------------------------------------------------------------!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       implicit none
       integer,intent(in)     :: Natoms  ! Number of atoms
       integer,intent(in)     :: Nbasis  ! Number of basis
       integer,intent(in)     :: nbs,nbp ! 
       integer,intent(in)     :: Nconts(Nbasis)
                                         ! Number of contractionbs
       integer,intent(in)     :: Ncontm  ! Max contractionbs

       complex*16,intent(in)  :: MatT(Nbasis,Nbasis)
       complex*16,intent(in)  :: MatD(Nbasis,Nbasis)
       real*8,intent(in)      :: nucpos(3,Natoms),nucvel(3,Natoms)
       real*8,intent(in)      :: alpha(Ncontm,Nbasis)
       real*8,intent(in)      :: coefs(Ncontm,Nbasis)
       integer,intent(in)     :: nucof(Nbasis)

       complex*16,intent(out) :: Bmat(Nbasis,Nbasis)
       complex*16,intent(out) :: force(3,Natoms)

       real*8     :: IMTX(3,4,4)
       complex*16 ::itrm(3)
       real*8     :: shrd,cij
       integer    :: ii,ni,nki,jj,nj,nkj
!--------------------------------------------------------------------!

       do jj=1,nbs
       nkj=nucof(jj)
       do nj=1,Nconts(jj)

         do ii=1,nbs
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xs>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ys>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zs>

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xs>

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ys>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zs>

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xs>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xs>

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ys>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ys>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zs>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zs>

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

       enddo
       enddo


       do jj=nbs+1,nbs+nbp,3
       nkj=nucof(jj)
       do nj=1,Nconts(jj)

         do ii=1,nbs
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)

          do ni=1,Nconts(ii)
            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ypx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zpx>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xpy> 
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zpy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xpz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ypz>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,1) !

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xpx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ypx>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zpx>


            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xpy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xpy> 
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xpy> 

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ypy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zpy>


            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xpz>

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ypz>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|zpz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,1) !


          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xpx>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xpx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xpx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xpx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ypx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ypx>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zpx>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zpx>


            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|xpy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|xpy> 
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|xpy> 
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|xpy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|xpy> 
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|xpy> 

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ypy>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ypy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ypy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ypy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|zpy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|zpy>


            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|xpz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|xpz>

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|ypz>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|ypz>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|zpz>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|zpz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|zpz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|zpz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

       enddo
       enddo



       do jj=nbs+nbp+1,Nbasis,6
       nkj=nucof(jj)
       do nj=1,Nconts(jj)

         do ii=1,nbs
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <s|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !

            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <s|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <s|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <s|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <s|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <s|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|ydzx>

            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <s|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <s|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <s|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <s|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <s|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <s|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,1,1) !

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+1,nbs+nbp,3
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <px|xdxx>
     >                   -2* cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,1,1) ! <py|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,2,1) ! <pz|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|xdyy>

            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|xdzz>

            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|xdxy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|xdyz>

            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|xdzx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !


            itrm(2)=itrm(2)+shrd*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <px|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <py|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <pz|ydxx>

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <px|ydyy>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,1,1) ! <py|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,2,1) ! <pz|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <px|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <py|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <pz|ydzz>

            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <px|ydxy>
     >                     - cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <py|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <pz|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <px|ydyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <py|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <pz|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|ydzx>


            itrm(3)=itrm(3)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <px|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <py|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <pz|zdxx>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <px|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <py|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <pz|zdyy>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <px|zdzz>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,4) ! <py|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,4) ! <pz|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <px|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <py|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <pz|zdxy>

            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <px|zdyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <py|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <pz|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <px|zdzx>
     >                     - cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <py|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <pz|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

         do ii=nbs+nbp+1,Nbasis,6
          nki=nucof(ii)
          itrm(:)=CMPLX(0.0d0,0.0d0)
          do ni=1,Nconts(ii)

            cij=coefs(ni,ii)*coefs(nj,jj)
            shrd=2*cij*alpha(nj,jj)
            call setim(0,1,alpha(ni,ii),alpha(nj,jj),
     >                 nucpos(:,nki),nucpos(:,nkj),IMTX)

            itrm(1)=itrm(1)+shrd*IMTX(1,3,4)*IMTX(2,1,1)*IMTX(3,1,1) ! <dxx|xdxx>
     >                   -2* cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,3,1)*IMTX(3,1,1) ! <dyy|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,1,1)*IMTX(3,3,1) ! <dzz|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,4)*IMTX(2,2,1)*IMTX(3,1,1) ! <dxy|xdxx>
     >                   -2* cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,4)*IMTX(2,2,1)*IMTX(3,2,1) ! <dyz|xdxx>
     >                   -2* cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,4)*IMTX(2,1,1)*IMTX(3,2,1) ! <dzx|xdxx>
     >                   -2* cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <dxx|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,3)*IMTX(3,1,1) ! <dyy|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,3,1) ! <dzz|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <dxy|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,2,1) ! <dyz|xdyy>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <dzx|xdyy>

            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|xdzz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|xdzz>

            itrm(1)=itrm(1)+shrd*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|xdxy>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|xdxy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|xdxy>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|xdxy>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(1)=itrm(1)+shrd*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|xdyz>
            itrm(1)=itrm(1)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|xdyz>

            itrm(1)=itrm(1)+shrd*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|xdzx>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|xdzx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|xdzx>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            itrm(1)=itrm(1)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|xdzx>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !


            itrm(2)=itrm(2)+shrd*IMTX(1,3,3)*IMTX(2,1,2)*IMTX(3,1,1) ! <dxx|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,3,2)*IMTX(3,1,1) ! <dyy|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,1,2)*IMTX(3,3,1) ! <dzz|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,3)*IMTX(2,2,2)*IMTX(3,1,1) ! <dxy|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,3)*IMTX(2,2,2)*IMTX(3,2,1) ! <dyz|ydxx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,3)*IMTX(2,1,2)*IMTX(3,2,1) ! <dzx|ydxx>

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,4)*IMTX(3,1,1) ! <dxx|ydyy>
     >                   -2* cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,4)*IMTX(3,1,1) ! <dyy|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,4)*IMTX(3,3,1) ! <dzz|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,4)*IMTX(3,1,1) ! <dxy|ydyy>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,4)*IMTX(3,2,1) ! <dyz|ydyy>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,4)*IMTX(3,2,1) ! <dzx|ydyy>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <dxx|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) ! <dyy|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) ! <dzz|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <dxy|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) ! <dyz|ydzz>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <dzx|ydzz>

            itrm(2)=itrm(2)+shrd*IMTX(1,3,2)*IMTX(2,1,3)*IMTX(3,1,1) ! <dxx|ydxy>
     >                     - cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,3,3)*IMTX(3,1,1) ! <dyy|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,3)*IMTX(3,3,1) ! <dzz|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,2,3)*IMTX(3,1,1) ! <dxy|ydxy>
     >                     - cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,3)*IMTX(3,2,1) ! <dyz|ydxy>
     >                     - cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,3)*IMTX(3,2,1) ! <dzx|ydxy>
     >                     - cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

            itrm(2)=itrm(2)+shrd*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <dxx|ydyz>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,2) ! <dyy|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,2) ! <dzz|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <dxy|ydyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,2) ! <dyz|ydyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            itrm(2)=itrm(2)+shrd*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <dzx|ydyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            itrm(2)=itrm(2)+shrd*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|ydzx>
            itrm(2)=itrm(2)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|ydzx>


            itrm(3)=itrm(3)+shrd*IMTX(1,3,3)*IMTX(2,1,1)*IMTX(3,1,2) ! <dxx|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,3,1)*IMTX(3,1,2) ! <dyy|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,1,1)*IMTX(3,3,2) ! <dzz|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,3)*IMTX(2,2,1)*IMTX(3,1,2) ! <dxy|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,3)*IMTX(2,2,1)*IMTX(3,2,2) ! <dyz|zdxx>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,3)*IMTX(2,1,1)*IMTX(3,2,2) ! <dzx|zdxx>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,3)*IMTX(3,1,2) ! <dxx|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,3)*IMTX(3,1,2) ! <dyy|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,3)*IMTX(3,3,2) ! <dzz|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,3)*IMTX(3,1,2) ! <dxy|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,3)*IMTX(3,2,2) ! <dyz|zdyy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,3)*IMTX(3,2,2) ! <dzx|zdyy>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,4) ! <dxx|zdzz>
     >                   -2* cij*IMTX(1,3,1)*IMTX(2,1,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,4) ! <dyy|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,3,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,4) ! <dzz|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,1,1)*IMTX(3,3,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,4) ! <dxy|zdzz>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,2,1)*IMTX(3,1,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,4) ! <dyz|zdzz>
     >                   -2* cij*IMTX(1,1,1)*IMTX(2,2,1)*IMTX(3,2,2) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,4) ! <dzx|zdzz>
     >                   -2* cij*IMTX(1,2,1)*IMTX(2,1,1)*IMTX(3,2,2) !

            itrm(3)=itrm(3)+shrd*IMTX(1,3,2)*IMTX(2,1,2)*IMTX(3,1,2) ! <dxx|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,3,2)*IMTX(3,1,2) ! <dyy|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,2)*IMTX(3,3,2) ! <dzz|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,2,2)*IMTX(3,1,2) ! <dxy|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,2)*IMTX(3,2,2) ! <dyz|zdxy>
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,2)*IMTX(3,2,2) ! <dzx|zdxy>

            itrm(3)=itrm(3)+shrd*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,3) ! <dxx|zdyz>
     >                     - cij*IMTX(1,3,1)*IMTX(2,1,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,3) ! <dyy|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,3,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,3) ! <dzz|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,1,2)*IMTX(3,3,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,3) ! <dxy|zdyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,2,2)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,3) ! <dyz|zdyz>
     >                     - cij*IMTX(1,1,1)*IMTX(2,2,2)*IMTX(3,2,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,3) ! <dzx|zdyz>
     >                     - cij*IMTX(1,2,1)*IMTX(2,1,2)*IMTX(3,2,1) !

            itrm(3)=itrm(3)+shrd*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,3) ! <dxx|zdzx>
     >                     - cij*IMTX(1,3,2)*IMTX(2,1,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,3) ! <dyy|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,3,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,3) ! <dzz|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,1,1)*IMTX(3,3,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,3) ! <dxy|zdzx>
     >                     - cij*IMTX(1,2,2)*IMTX(2,2,1)*IMTX(3,1,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,3) ! <dyz|zdzx>
     >                     - cij*IMTX(1,1,2)*IMTX(2,2,1)*IMTX(3,2,1) !
            itrm(3)=itrm(3)+shrd*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,3) ! <dzx|zdzx>
     >                     - cij*IMTX(1,2,2)*IMTX(2,1,1)*IMTX(3,2,1) !

          enddo
          call sumit(itrm,MatT(ii,jj),MatD(ii,jj),nucvel(:,nkj),
     >    force(:,nkj),Bmat(ii,jj))
         enddo

       enddo
       enddo


       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
