c     ===========================================
c BROYDEN-FLETCHER-GOLDFARB-SHANNO SUBROUTINE,
c using search, for line minimizations, no bracketing,
c searchs for 0 of gradient
c
c----------------------------------------------------------

      SUBROUTINE dfp2(xH)
      use garcha_mod
      implicit real*8 (a-h,o-z)
c
c  STPMX 0.05 original
      real*8,intent(out) :: xH(ntatom*3,ntatom*3)

      real*8, dimension(:,:), ALLOCATABLE ::f,dxyz
      real*8, dimension(:), ALLOCATABLE :: P,XI,G,DG,HDG,PNEW
      real*8 dipxyz(3)
      PARAMETER (STPMX=0.1D0,ITMAX=200,EPS=3.E-08)
      PARAMETER (TOLX=4.D0*EPS)
c
      character*2 atname
      gtol=1D-10
c
      nt3=ntatom*3
      allocate(f(ntatom,3),dxyz(3,ntatom),P(nt3*3),XI(nt3*3)
     > ,G(nt3*3),DG(nt3),HDG(nt3),PNEW(nt3))

c
      iforce=0
c
      if (field) then
       g0=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
      endif
c
      zero=0.0D0
      Nel=2*NCO+Nunp
      write(*,*) 'NEW DAVIDSON-FLETCHER-POWELL OPTIMIZATION'
      write(*,*)
c
      ntom=natom-nfrozen
c
      if (sol.and.free)then
       ntom=ntom+Nsol*natsol
      endif

      ntom=ntatom ! ESTO HABRA QUE CAMBIARLO
c
      N=3*ntom
       
c
      k=0
      do 20 i=1,ntom
        k=k+1
        P(k)=r(i,1)
        k=k+1
        P(k)=r(i,2)
        k=k+1
        P(k)=r(i,3)
 20   continue
c
c
c
c first call , needs also E , not only gradients -----------
c

c#ifdef G2G
c			write(*,*) 'primeraraa carga de posiciones (con igrid2)'
c			call g2g_reload_atom_positions(igrid2)
c#endif
!       rqm=r
      GRAD=.false.
      if (OPEN) then
      call SCFop()
       else
      call SCF(FP,dipxyz)
      endif
c
c now gradients 
c
c
c#ifdef G2G
c			write(*,*) 'cambio de grilla para fuerza+energia (igrid)'
c      call g2g_new_grid(igrid)
c#endif
       dxyz=0
       call dft_get_qm_forces(dxyz)

       do inan=1,natom
       f(inan,1)=-dxyz(1,inan)
       f(inan,2)=-dxyz(2,inan)
       f(inan,3)=-dxyz(3,inan)

       enddo
c        write(1234,*) f
c       stop
c      call int1G(f)
c
c
c      call int3G(f,.true.)
c
c      call intSG(f)
c
c reaction field case -------
      if (field) then
        g1=g0
        call dip()
c
c        call dipg(f)
       write(*,*) 'field and geometry opt are not implemented yet'

      endif
c----------------------------
c classical solvent case ----
        if (nsol.gt.0) then
        stop 'optim with nsol/= 0 not suported yet'
c         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f)
         call intsolG(f)

        endif
c---------------------------
c      FP=FP+Exc
c ----------------------------------------------------------
      write(*,*)
c
      write(*,*) 'ENERGY GRADIENTS, IN A.U'
      k=0
      ss=0.D0
      ss1=0.0D0
      do 29 i=1,ntom
        ss=ss+f(i,1)**2+f(i,2)**2+f(i,3)**2
        ss1=ss1+r(i,1)**2+r(i,2)**2+r(i,3)**2
        k=k+1
        G(k)=f(i,1)
        XI(k)=-f(i,1)
        k=k+1
        G(k)=f(i,2)
        XI(k)=-f(i,2)
        k=k+1 
        G(k)=f(i,3)
        XI(k)=-f(i,3)
        write(*,500) i,f(i,1),f(i,2),f(i,3)
 29    continue
c
      if (ibrent.eq.2) then
      do 31 i=ntom+1,natom
        write(*,500) i,zero,zero,zero
 31   continue
      endif
c
      stpmax=STPMX*dmax1(dsqrt(ss1),dfloat(N))
      write(*,450) sqrt(ss)
      write(*,*)
      if (ss.lt.1.D-08) then
       write(*,*) 'INITIAL GEOMETRY ALREADY OPTIMIZED'
      deallocate(f,dxyz,P,XI,G,DG,HDG,PNEW)

       return
      endif
c
      DO 12 I=1,N
        DO 11 J=1,N
          xH(I,J)=0.D0
11      CONTINUE
        xH(I,I)=1.D0
12    CONTINUE
c

      DO 27 ITS=1,ITMAX
        ITER=ITS
c
        CALL LSEARCH(N,P,FP,G,XI,PNEW,  
     >               FRET,STPMAX,CHECK)
c
        write(*,*) PNEW
        FP=FRET
       open(unit=69,file='opt.xyz')
         A0=0.5291771
        write(69,*) int(N/3)
        write(69,*)
       atname='??' 
        do 121 i=1,N,3
       ii=(i/3)+1
      if(Iz(ii).eq.1) atname='H'
      if(Iz(ii).eq.6) atname='C'
      if(Iz(ii).eq.7) atname='N'
      if(Iz(ii).eq.8) atname='O'
      if(Iz(ii).eq.3) atname='Li'
       write(69,*)  Iz(ii),PNEW(i)*A0,PNEW(i+1)*A0,PNEW(i+2)*A0
 121     write(*,501) PNEW(i),PNEW(i+1),PNEW(i+2)
c
        if (ibrent.eq.2) then
        do 122 i=natom-nfrozen+1,natom
         write(*,501) r(i,1),r(i,2),r(i,3)
 122    continue
        endif
c
        DO 13 I=1,N
         XI(I)=PNEW(I)-P(I)
         P(I)=PNEW(I)
 13     CONTINUE
c
        k=0
        do 30 i=1,ntom
        k=k+1
        r(i,1)=P(k)
        k=k+1
        r(i,2)=P(k)
        k=k+1
        r(i,3)=P(k)
30      continue
c
        test=0.0D0
c
      DO 14 i=1,N
        temp=abs(XI(i))/dmax1(abs(P(i)),1.D0)
        if (temp.gt.test) test=temp
 14   CONTINUE
c
       IF (test.lt.TOLX) then
       write(*,*) 'test & tolx', test,TOLX
      deallocate(f,dxyz,P,XI,G,DG,HDG,PNEW)
       return
       endif
c
c
        DO 15 I=1,N
          DG(I)=G(I)
15      CONTINUE

#ifdef G2G
c			write(*,*) 'actualizacion de posiciones por movimiento'
c			call g2g_reload_atom_positions(igrid2)
#endif
c
      GRAD=.false.
c       rqm=r
      if (OPEN) then
      call SCFop()
       else
      call SCF(E,dipxyz)

      endif
c
c now gradients
c
c
c#ifdef G2G
c			write(*,*) 'cambio de grilla para fuerza+energia (igrid)'
c      call g2g_new_grid(igrid)
c#endif

       call dft_get_qm_forces(dxyz)

       do inan=1,natom
       f(inan,1)=dxyz(1,inan)
       f(inan,2)=dxyz(2,inan)
       f(inan,3)=dxyz(3,inan)
       write(987,*) inan,f(inan,1),f(inan,2),f(inan,3) 
       enddo

c      call int1G(f)
c

c       call int3G(f,true)
c
c      call intSG(f)
c reaction field case ------
      if (field) then
c        call dip(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,Nel,
c     >       ux,uy,uz)
c
c        g1=g0
c        call dipg(NORM,Iz,natom,r,Nuc,M,ncont,nshell,c,a,RMM,
c     >               Nel,g1,ux,uy,uz,f)
       write(*,*) 'field not posible'
        stop
      endif
c
c----------------------------
c classical solvent case ----
        if (sol) then
c         call mmsolG(natom,Nsol,natsol,Iz,pc,r,Em,Rm,f)
c         call intsolG(f)
        endif
c---------------------------
c
c      E=E+Exc
      write(*,480) E
c ----------------------------------------------------------
c
c
        k=0
        ss=0.0D0
        do 140 i=1,ntom
        ss=ss+f(i,1)**2+f(i,2)**2+f(i,3)**2
        k=k+1
        G(k)=f(i,1)
        k=k+1
        G(k)=f(i,2)
        k=k+1 
        G(k)=f(i,3)
c test
c        write(*,500) i,f(i,1),f(i,2),f(i,3)       
140   continue
c
      test=0.0D0
c

      den=dmax1(fret,1.D0)
      write(*,450) sqrt(ss)
c
        DO 16 I=1,N
         temp=abs(g(i))*dmax1(abs(P(i)),1.D0)/den
         if (temp.gt.test) test=temp
 16     CONTINUE
c
       if (test.lt.gtol) then
       write(*,*) 'AAAAAAAAAAAAAAAAA' 
      deallocate(f,dxyz,P,XI,G,DG,HDG,PNEW)
         return
       endif
c
        DO 17 i=1,N
         DG(i)=G(i)-DG(i)
 17     CONTINUE
c
        DO 19 I=1,N
          HDG(I)=0.D0
          DO 18 J=1,N
18           HDG(I)=HDG(I)+xH(I,J)*DG(J)
19        CONTINUE
        FAC2=0.D0
        FAE=0.D0
        SUMDG=0.0D0
        SUMXI=0.0D0
c
        DO 21 I=1,N
          FAC2=FAC2+DG(I)*XI(I)
          FAE=FAE+DG(I)*HDG(I)
          SUMDG=SUMDG+DG(I)**2
          SUMXI=SUMXI+XI(I)**2
21      CONTINUE
c
        if ((FAC2**2).gt.(EPS*SUMDG*SUMXI)) then
        FAC2=1.D0/FAC2
        FAD=1.D0/FAE
c
        DO 22 I=1,N
          DG(I)=FAC2*XI(I)-FAD*HDG(I)
22      CONTINUE
c
        DO 24 I=1,N
          DO 23 J=1,N
            xH(I,J)=xH(I,J)+FAC2*XI(I)*XI(J)
     *        -FAD*HDG(I)*HDG(J)+FAE*DG(I)*DG(J)
23        CONTINUE
24      CONTINUE
        ENDIF
c
        DO 26 I=1,N
          XI(I)=0.D0
          DO 25 J=1,N
            XI(I)=XI(I)-xH(I,J)*G(J)
25        CONTINUE
26      CONTINUE
27    CONTINUE
c      PAUSE 'too many iterations in DFPMIN'
      write(*,*) 'too many iterations in DFPMIN'
      deallocate(f,dxyz,P,XI,G,DG,HDG,PNEW)
      return
      
c
 450  format ('NORM OF GRADIENT ',F17.10)
 480  format ('SCF ENERGY DFP2= ',F17.10)
 500  format (I3,2x,3(F17.10,2x))
 501  format (3(F12.6,2x))
      RETURN
      END
c
c     ===========================================================
