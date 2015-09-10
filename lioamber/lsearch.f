      SUBROUTINE LSEARCH(N,xold,fold,G,P,X1,F,STPMAX,CHECK)
       use garcha_mod

c     !aca P en la fuerza (o el gradiente, no se bien), al menos en le primer paso, despues es el deltaX
c     ! xold son las posiciones viejas
c     ! fold es la energia vieja, F la nueva
c     ! G es la fuerza 
c
      implicit real*8 (a-h,o-z)
      REAL*8 , intent(inout)::f,fold,stpmax,g(n),p(n),x1(n),xold(n)
      real*8 dipxyz(3)

c      real*8 f,fold,stpmax,g(n),p(n),x1(n),xold(n),dipxyz(3)
      PARAMETER (ALF=1.D-03,TOLX=1.D-04)
      LOGICAL, intent(inout) :: CHECK
      logical write
c
      ntom=natom-nfrozen
c
       if (sol.and.free) then
       ntom=ntom+Nsol*natsol
       endif
       ntom=natom ! esto hay que tocarlo!!!!
c
      check=.false.
      sum=0.0D0
      do  i=1,n
       sum=sum+P(i)**2
      enddo
      sum=sqrt(sum)
c
      
      if (sum.gt.stpmax) then
       do 12 i=1,N
        p(i)=p(i)*stpmax/sum
 12    continue
      endif
c
      slope=0.0D0
      do 13 i=1,N
       slope=slope+g(i)*p(i)
 13   continue
c
      test=0.0D0
      do 14 i=1,N
        temp=abs(p(i))/dmax1(dabs(xold(i)),1.D0)
        if (temp.gt.test) test=temp ! test es el máximo Dx/X
 14   continue
c

      alamin=TOLX/test    ! 
      alam=1.D0
c 
      do kjta=1,100
c
      do 15 i=1,N
       x1(i)=xold(i)+alam*p(i) ! aca mueve las coordenadas
 15   continue
c
      k=0
		  write(*,*) 'movimiento de posiciones por lsearch (con igrid2)'      
      do i=1,ntom
       k=k+1
       r(i,1)=x1(k)
       k=k+1
       r(i,2)=x1(k)
       k=k+1
       r(i,3)=x1(k)
      enddo
       
c      do i=1,natom
c       rqm(i,:)=r(i,:)
c      enddo 
c	  	call g2g_reload_atom_positions(igrid2)
c
      GRAD=.true.
c      rqm=r
      if (OPEN) then
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,F,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      else
      call SCF(F,dipxyz)

      endif

       write(*,*) 'energias en lsearch, nueva y vieja',F,FOLD
c
      if (alam.lt.alamin) then
       write(*,*) 'alamy alamin',alam,alamin
       do 16 i=1,N
        x1(i)=xold(i)
 16    continue
c
       check=.true.
       return
c
       else if (F.le.fold+ALF*alam*slope) then
        return
       else
        if (alam.eq.1.D0) then
        tmplam=-slope/(2.D0*(f-fold-slope))
        else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold2-alam2*slope
        a1=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b1=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if (a1.eq.0.D0) then
         tmplam=-slope/(2.D0*b1)
        else
         disc=b1*b1-3.D0*a1*slope
         tmplam=(-b1+sqrt(disc))/(3.D0*a1)
        endif
c
        if (tmplam.gt.0.5D0*alam) tmplam=.5D0*alam
       endif
       endif
c
       alam2=alam
       f2=f
       fold2=fold
       alam=dmax1(tmplam,0.1D0*alam)


       enddo

c
 500  format (I3,2x,3(F12.5,2x))
c
       END
 
