!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine setim(ni_in,nj_in,ai_in,aj_in,ri_in,rj_in,IntMat)
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
       integer,intent(in) :: ni_in,nj_in
       real*8,intent(in)  :: ai_in,aj_in
       real*8,intent(in)  :: ri_in(3),rj_in(3)
       real*8,intent(out) :: IntMat(3,4,4)

       integer            :: ni,nj
       real*8             :: ai,aj
       real*8             :: ri(3),rj(3)

       logical             :: dotranspose
       integer             :: ntot,ii,jj,kk
       real*8              :: a1,b1,b2,b3,b5,b7,sqrtpi
       real*8              :: theta(3),factor(3),r0(3)
       real*8              :: BI0,BI2,BI4,BI6
       real*8,dimension(3) :: DA1,DA2,DA3,DB1,DB2,DB3
       real*8              :: term0,term2,term4


       real*8,parameter :: PI=3.14159265358979323846264338327



! If asking for a tall matrix, calculate the long matrix instead
! (and later transpose it)
!--------------------------------------------------------------------!
       dotranspose=.false.
       if (nj.ge.ni) then
         ni=ni_in
         nj=nj_in
         ai=ai_in
         aj=aj_in
         ri=ri_in
         rj=rj_in
         dotranspose=.false.
       else
         ni=nj_in
         nj=ni_in
         ai=aj_in
         aj=ai_in
         ri=rj_in
         rj=ri_in
         dotranspose=.true.
       endif


! Set initial parameters according to the matrix required
!--------------------------------------------------------------------!
       a1=ai+aj
       b1=sqrt(1/a1)
       sqrtpi=sqrt(PI)
       BI0=sqrtpi*b1

       ntot=ni+nj
       if (ntot.ge.2) then
         b2=b1*b1
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

       do kk=1,3
         theta(kk)=ri(kk)-rj(kk)
         theta(kk)=theta(kk)**2
         theta(kk)=theta(kk)*ai*aj
         theta(kk)=theta(kk)*b2
         factor(kk)=exp(theta(kk))

         r0(kk)=ai*ri(kk)+aj*rj(kk)
         r0(kk)=r0(kk)*b2

         DA1(kk)=r0(kk)-ri(kk)
         if (ni.ge.2) then
           DA2(kk)=DA1(kk)*DA1(kk)
           if (ni.eq.3) DA3(kk)=DA2(kk)*DA1(kk)
         endif
         
         DB1(kk)=r0(kk)-rj(kk)
         if (nj.ge.2) then
           DB2(kk)=DB1(kk)*DB1(kk)
           if (nj.eq.3) DB3(kk)=DB2(kk)*DB1(kk)
         endif
       enddo


! Calculate the actual matrixes
!--------------------------------------------------------------------!
       IntMat(:,:,:)=0.0d0
       do kk=1,3
         IntMat(kk,1,1)=BI0
         IntMat(kk,2,1)=BI0*DA1(kk)
         IntMat(kk,1,2)=BI0*DB1(kk)
         if (ntot.ge.2) IntMat(2,2,kk)=BI0*DA1(kk)*DB1(kk)+BI2



         if (nj.ge.2) then
           IntMat(kk,1,3)=BI0*DA2(kk)+BI2
           IntMat(kk,2,3)=BI0*DA1(kk)*DB2(kk)+BI2*(DA1(kk)+2*DB2(kk))

           if (ni.ge.2) then
             IntMat(kk,3,1)=BI0*DA2(kk)+BI2
             IntMat(kk,3,2)=BI0*DA2(kk)*DB1(kk)+BI2*(2*DA1(kk)+DB1(kk))

             term0=BI0*DA2(kk)*DB2(kk)
             term2=BI2*(DA2(kk)+4*DA1(kk)*DB1(kk)+DB2(kk))
             IntMat(kk,3,3)=term0+term2+BI4

           endif

           if (nj.ge.3) then
             IntMat(kk,1,4)=BI0*DB3(kk)+BI2*DB1(kk)*2

             term0=BI0*DA1(kk)*DB3(kk)
             term2=BI2*(3*DA1(kk)*DB1(kk)+3*DB2(kk))
             IntMat(kk,2,4)=term0+term2+BI4

             if (ni.ge.2) then
               term0=BI0*DA2(kk)*DB3(kk)
               term2=BI2*(3*DA2(kk)*DB1(kk)+6*DA1(kk)*DB2(kk)+DB3(kk))
               term4=BI4*(2*DA1(kk)+3*DB1(kk))
               IntMat(kk,3,4)=term0+term2+term4

               if (ni.ge.3) then
                 IntMat(kk,4,1)=BI0*DA3(kk)+BI2*DA1(kk)*2

                 term0=BI0*DA3(kk)*DB1(kk)
                 term2=BI2*(3*DA2(kk)+3*DA1(kk)*DB1(kk))
                 IntMat(kk,4,2)=term0+term2+BI4

                 term0=BI0*DA3(kk)*DB2(kk)
                 term2=BI2*(DA3(kk)+6*DA2(kk)*DB1(kk)+3*DA1(kk)*DB2(kk))
                 term4=BI4*(3*DA1(kk)+2*DB1(kk))
                 IntMat(kk,4,3)=term0+term2+term4

                 term0=BI0*DA3(kk)*DB3(kk)+BI6
                 term2=0
                 term2=term2+BI2*3*DA3(kk)*DB1(kk)
                 term2=term2+BI2*9*DA2(kk)*DB2(kk)
                 term2=term2+BI2*3*DA1(kk)*DB3(kk)
                 term4=BI4*(3*DA2(kk)+9*DA1(kk)*DB1(kk)+3*DB2(kk))
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
