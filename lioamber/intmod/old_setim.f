!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine setim(samenuc,maxni,maxnj,alphai,alphaj,ri,rj,
     >                  IntMat)
!--------------------------------------------------------------------!
!
! IntMat(i,j) contains the integral of the product fa[i]*fb[j]
! between -infinity and +infinity, where:
!
!   fk[n](x) = (x-x[k])^n * exp( -alp[k] * (x-x[k])^2 )
!
! The matrix has an extra dimension of size 3 because there is
! one matrix for each of the 3 dimensions; the variables ai/aj
! are the same for all directions, but the x[k], y[k] and z[k]
! are different and are all introduced in ri[k],rj[k].
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       implicit none
       logical,intent(in) :: samenuc
       integer,intent(in) :: maxni,maxnj
       real*8,intent(in)  :: alphai,alphaj
       real*8,intent(in)  :: ri(3),rj(3)
       real*8,intent(out) :: IntMat(3,4,4)

       integer            :: mi,mj
       real*8             :: ai,aj
       real*8             :: posi(3),posj(3)

       logical             :: dotranspose
       integer             :: ntot,ii,jj,kk
       real*8              :: a1,b1,b2,b3,b5,b7,sqrtpi
       real*8              :: diff,fusion
       real*8              :: theta,factor(3),r0(3)
       real*8              :: BI0,BI2,BI4,BI6
       real*8,dimension(3) :: DA1,DA2,DA3,DB1,DB2,DB3
       real*8              :: term0,term2,term4

       real*8,parameter :: PI=3.14159265358979323846264338327



! If asking for a tall matrix, calculate the long matrix instead
! (and later transpose it)
!--------------------------------------------------------------------!
       if (maxnj.ge.maxni) then
         mi=maxni
         mj=maxnj
         ai=alphai
         aj=alphaj
         posi=ri
         posj=rj
         dotranspose=.false.
       else
         mi=maxnj
         mj=maxni
         ai=alphaj
         aj=alphai
         posi=rj
         posj=ri
         dotranspose=.true.
       endif


! Set initial parameters according to the matrix required
!--------------------------------------------------------------------!
       a1=ai+aj
       b2=(1/a1)
       b1=sqrt(b2)

       sqrtpi=sqrt(PI)
       BI0=sqrtpi*b1

       ntot=mi+mj
       if (ntot.ge.2) then
         b3=b2*b1
         BI2=sqrtpi*b3*(0.5)

         if (ntot.ge.4) then
           b5=b3*b2
           BI4=sqrtpi*b5*(0.75)

           if (ntot.eq.6) then
             b7=b5*b2
             BI6=sqrtpi*b7*(1.875)

           endif
         endif
       endif


       if (samenuc) then
       do kk=1,3
         DA1(kk)=0.0d0
         DB1(kk)=0.0d0
         factor(kk)=1.0d0
       enddo
       else
       do kk=1,3
         r0(kk)=ai*posi(kk)+aj*posj(kk)
         r0(kk)=r0(kk)*b2
         DA1(kk)=r0(kk)-posi(kk)
         DB1(kk)=r0(kk)-posj(kk)

         theta=posi(kk)-posj(kk)
         theta=theta**2
         theta=theta*ai*aj
         theta=theta*b2
         factor(kk)=exp(-theta)
       enddo
       endif

       if (mi.ge.2) then
       do kk=1,3
         DA2(kk)=DA1(kk)*DA1(kk)
         if (mi.eq.3) DA3(kk)=DA2(kk)*DA1(kk)
       enddo
       endif

       if (mj.ge.2) then
       do kk=1,3
         DB2(kk)=DB1(kk)*DB1(kk)
         if (mj.eq.3) DB3(kk)=DB2(kk)*DB1(kk)
       enddo
       endif



! Calculate the actual matrixes
!--------------------------------------------------------------------!
       IntMat(:,:,:)=0.0d0
       do kk=1,3
         IntMat(kk,1,1)=BI0
         IntMat(kk,1,2)=BI0*DB1(kk)
         IntMat(kk,2,1)=BI0*DA1(kk)
         if (ntot.ge.2) IntMat(2,2,kk)=BI0*DA1(kk)*DB1(kk)+BI2

         if (mj.ge.2) then
           IntMat(kk,1,3)=BI0*DB2(kk)+BI2
           IntMat(kk,2,3)=BI0*DA1(kk)*DB2(kk)+BI2*(DA1(kk)+2*DB2(kk))

           if (mi.ge.2) then
             IntMat(kk,3,1)=BI0*DA2(kk)+BI2
             IntMat(kk,3,2)=BI0*DA2(kk)*DB1(kk)+BI2*(2*DA1(kk)+DB1(kk))

             term0=BI0*DA2(kk)*DB2(kk)
             term2=BI2*(DA2(kk)+4*DA1(kk)*DB1(kk)+DB2(kk))
             IntMat(kk,3,3)=term0+term2+BI4

           endif

           if (mj.ge.3) then
             IntMat(kk,1,4)=BI0*DB3(kk)+BI2*DB1(kk)*3

             term0=BI0*DA1(kk)*DB3(kk)
             term2=BI2*3*(DA1(kk)*DB1(kk)+DB2(kk))
             IntMat(kk,2,4)=term0+term2+BI4

             if (mi.ge.2) then
               term0=BI0*DA2(kk)*DB3(kk)
               term2=BI2*(3*DA2(kk)*DB1(kk)+6*DA1(kk)*DB2(kk)+DB3(kk))
               term4=BI4*(2*DA1(kk)+3*DB1(kk))
               IntMat(kk,3,4)=term0+term2+term4

               if (mi.ge.3) then
                 IntMat(kk,4,1)=BI0*DA3(kk)+BI2*DA1(kk)*3

                 term0=BI0*DA3(kk)*DB1(kk)
                 term2=BI2*3*(DA2(kk)+DA1(kk)*DB1(kk))
                 IntMat(kk,4,2)=term0+term2+BI4

                 term0=BI0*DA3(kk)*DB2(kk)
                 term2=BI2*(DA3(kk)+6*DA2(kk)*DB1(kk)+3*DA1(kk)*DB2(kk))
                 term4=BI4*(3*DA1(kk)+2*DB1(kk))
                 IntMat(kk,4,3)=term0+term2+term4

                 term0=BI0*DA3(kk)*DB3(kk)
                 term2=0
                 term2=term2+BI2*3*DA3(kk)*DB1(kk)
                 term2=term2+BI2*9*DA2(kk)*DB2(kk)
                 term2=term2+BI2*3*DA1(kk)*DB3(kk)
                 term4=BI4*3*(DA2(kk)+3*DA1(kk)*DB1(kk)+DB2(kk))
                 IntMat(kk,4,4)=term0+term2+term4+BI6

               endif
             endif
           endif
         endif


         do jj=1,4
         do ii=1,4
           IntMat(kk,ii,jj)=IntMat(kk,ii,jj)*factor(kk)
         enddo
         enddo

       enddo


! If originally  asked for a tall matrix, transpose the long one
!--------------------------------------------------------------------!
       if (dotranspose) then
         do kk=1,3
         do ii=2,4
         do jj=1,ii-1
           term0=IntMat(kk,ii,jj)
           IntMat(kk,ii,jj)=IntMat(kk,jj,ii)
           IntMat(kk,jj,ii)=term0
         enddo
         enddo
         enddo
       endif

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
