!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module garcha_mod
!------------------------------------------------------------------------------!
       implicit real*8 (a-h,o-z)

      INCLUDE 'param.f'
      logical verbose 
      integer M,Md,natom,ntatom,NMAX,NCO,NUNP,igrid,igrid2
     >  ,Iexch,nsol,npas,npasw,idip,watermod,noconverge,
     > converge,ndiis,NGEO,nang,timedep,ntdstep,propagator,NBCH 
      integer restart_freq, energy_freq
      real*8 GOLD, TOLD, qmmmcut, dgtrig
      parameter (nng=100)
      character*65 title
      character*20 basis,whatis,stdbas
      character*40 basis_set, fitting_set
      logical int_basis
      character*4 date
      character*20 output,fcoord,fmulliken,frestart,frestartin,solv,
     > solv2
      character*4 ctype
      logical done(ntq),exists,MEMO,predcoef
      logical used,NORM,OPEN,ATRHO,DIRECT,VCINP,SHFT,DIIS
      logical TMP1,TMP2,dens,EXTR,write1,SVD,ANG,field1
      logical Coul,Scf1,Prop,GRAD,BSSE,integ,SVD1,sol,tipe
      logical exter,exter1,resp1,popf,primera,writexyz,intsoldouble
      logical OPEN1
      logical dens1,integ1,sol1,free,free1, field, extern

      logical tdrestart, writedens

      logical cubegen_only,cube_dens,cube_orb,cube_elec
      integer cube_res,cube_sel
      character*20 cube_dens_file,cube_orb_file,cube_elec_file


      dimension OCC(40),oc2(400),ATCOEF(100*ng0),ighost(ntq),
     > ighost1(ntq)
      dimension ncf(nng),lt(nng)
      real*8 e_(50,3),wang(50),e_2(116,3),wang2(116),e3(194,3), ! intg1 e intg2
     > wang3(194)                                               !
      integer Nr(0:54),Nr2(0:54)
      real*8 Fx, Fy, Fz, epsilon, a0,tdstep

      real*8, dimension (:,:), ALLOCATABLE :: r,v,rqm,d
      real*8, dimension (:), ALLOCATABLE ::  Em, Rm, pc
       integer, dimension (:), ALLOCATABLE :: Iz, nnat

      dimension isotop(54)!,Pm(nt)
      dimension Rm2(0:54), STR(880,0:21), FAC(0:16)
      dimension alpha(nss)
c Everything is dimensioned for 2 basis, normal and density
c ncf, lt,at,ct parameters for atomic basis sets
      dimension at(nng),ct(nng),nshell(0:4)
      dimension Num(0:3),nlb(ng),nld(ngd),nshelld(0:4)
       integer iconst1,idip1,ipop1,ispin1,
     > icharge1,Nsol1,natsol1,Ll(3)
      
        real*8, dimension (:), ALLOCATABLE :: af
       real*8, dimension (:,:), ALLOCATABLE :: c,a,cx,ax,cd,ad,B
       integer, dimension (:), ALLOCATABLE :: Nuc,ncont,Nucx,ncontx,Nucd
     >  ,ncontd
        integer, dimension (:), ALLOCATABLE :: indexii, indexiid

c
       
       real*8, dimension (:), ALLOCATABLE :: RMM,RMM1,RMM2,RMM3
       real*8, dimension (:), ALLOCATABLE :: rhoalpha,rhobeta
       real*8, dimension (:,:), ALLOCATABLE :: X, XX
       real*8, dimension (:), ALLOCATABLE :: old1,old2,old3 

c
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0, pi5=34.9868366552497108D0,
     >    pi52=34.9868366552497108D0)
      parameter(pis32=5.56832799683170698E0,piss=3.14159265358979312E0,
     >          rpis=1.77245385090551588E0, pis5=34.9868366552497108E0,
     >    pis52=34.9868366552497108E0)


c Angular momenta : up to f functions ( could be easily extended if
c necessary)
c
c      common /fit/ Nang,dens1,integ1,Iexch1,igridc,igrid2c
c      common /cav/ a01,epsilon1,field1,exter1,Fx1,Fy1,Fz1
c
c      common /index/ ii,iid
c
      Data Num /1,3,6,10/
      dimension jatom(2,100),coef(100),dist(100,3),distt(100)
      integer ndis,nsteps
      real*8 kjar,xini,xfinal   

      integer, dimension(:), ALLOCATABLE :: natomc,nnps,nnpp,nnpd,nns
      integer, dimension(:), ALLOCATABLE :: nnd,nnp
      real*8, dimension (:), ALLOCATABLE :: atmin
      integer, dimension(:,:), ALLOCATABLE :: jatc
      integer kknums,kknumd
      integer, dimension (:), ALLOCATABLE :: kkind,kkinds
      real*8     rmax, rmaxs
      real*8, dimension (:), ALLOCATABLE :: cool
      real*4, dimension (:), ALLOCATABLE :: cools
c      parameter rmintsol=16.0D0
!
!------------------------------------------------------------------------------!
       real*8,allocatable,dimension(:,:)     :: Smat
       real*8,allocatable,dimension(:,:)     :: RealRho
!       real*8,allocatable,dimension(:,:)     :: Gmat !DK
!       real*8,allocatable,dimension(:,:)     :: Hmat !DK
!       real*8,allocatable,dimension(:,:)     :: FockMat
!       complex*16,allocatable,dimension(:,:) :: RhoOld,RhoNew

!       real*8,allocatable,dimension(:,:) :: Lmat,Linv,Umat,Uinv
!------------------------------------------------------------------------------!
       dimension xmass(216)
*     Atomic masses (u.m.a.) of most common isotopes
      data xmass /
*      H-1             H-2             H-3
     >  1.007825037D0,  2.014101787D0,  3.016049286D0,  0.0,
*      He-4            He-3
     >  4.00260325D0,   3.016029297D0,  0.0,            0.0,
*      Li-7            Li-6
     >  7.0160045D0,    6.0151232D0,    0.0,            0.0,
*      Be-9
     >  9.0121825D0,    0.0,            0.0,            0.0,
*       B-11            B-10
     > 11.0093053D0,   10.0129380D0,    0.0,            0.0,
*       C-12            C-13
     > 12.000000000D0, 13.003354839D0,  0.0,            0.0,
*       N-14            N-15
     > 14.003074008D0, 15.000108978D0,  0.0,            0.0,
*       O-16            O-18            O-17
     > 15.99491464D0,  17.99915939D0,  16.9991306D0,    0.0,
*       F-19
     > 18.99840325D0,   0.0,            0.0,            0.0,
*      Ne-20           Ne-22           Ne-21
     > 19.9924391D0,   21.9913837D0,   20.9938453D0,    0.0,
*      Na-23
     > 22.9897697D0,    0.0,            0.0,            0.0,
*      Mg-24           Mg-26           Mg-25
     > 23.9850450D0,   25.9825954D0,   24.9858392D0,    0.0,
*      AL-27
     > 26.9815413D0,    0.0,            0.0,            0.0,
*      Si-28           Si-29           Si-30
     > 27.9769284D0,   28.9764964D0,   29.9737717D0,    0.0,
*       P-31
     > 30.9737634D0,    0.0,            0.0,            0.0,
*       S-32            S-34            S-33            S-36
     > 31.9720718D0,   33.96786774D0,  32.9714591D0,   35.9670790D0,
*      Cl-35           Cl-37
     > 34.968852729D0, 36.965902624D0,  0.0,            0.0,
*      Ar-40           Ar-36           Ar-38
     > 39.9623831D0,   35.967545605D0, 37.9627322D0,    0.0,
*       K-39            K-41            K-40
     > 38.9637079D0,   40.9618254D0,   39.9639988D0,    0.0,
*      Ca-40           Ca-44           Ca-42           Ca-48
     > 39.9625907D0,   43.9554848D0,   41.9586218D0,   47.952532D0,
*      Sc-45
     > 44.9559136D0,    0.0,            0.0,            0.0,
*      Ti-48           Ti-46           Ti-47           Ti-49
     > 47.9479467D0,   45.9526327D0,   46.9517649D0,   48.9478705D0,
*       V-51            V-50
     > 50.9439625D0,   49.9471613D0,    0.0,            0.0,
*      Cr-52           Cr-53           Cr-50           Cr-54
     > 51.9405097D0,   52.9406510D0,   49.9460463D0,   53.9388822D0,
*      Mn-55
     > 54.9380463D0,    0.0,            0.0,            0.0,
*      Fe-56           Fe-54           Fe-57           Fe-58
     > 55.9349393D0,   53.9396121D0,   56.9353957D0,   57.9332778D0,
*      Co-59
     > 58.9331978D0,    0.0,            0.0,            0.0,
*      Ni-58           Ni-60           Ni-62           Ni-61
     > 57.9353471D0,   59.9307890D0,   61.9283464D0,   60.9310586D0,
*      Cu-63           Cu-65
     > 62.9295992D0,   64.9277924D0,    0.0,            0.0,
*      Zn-64           Zn-66           Zn-68           Zn-67
     > 63.9291454D0,   65.9260352D0,   67.9248458D0,   66.9271289D0,
*      Ga-69           Ga-71
     > 68.9255809D0,   70.9247006D0,    0.0,            0.0,
*      Ge-74           Ge-72           Ge-70           Ge-73
     > 73.9211788D0,   71.9220800D0,   69.9242498D0,   72.9234639D0,
*      As-75
     > 74.9215955D0,    0.0,            0.0,            0.0,
*      Se-80           Se-78           Se-82           Se-76
     > 79.9165205D0,   77.9173040D0,   81.916709D0,    75.9192066D0,
*      Br-79           Br-81
     > 78.9183361D0,   80.916290D0,     0.0,            0.0,
*      Kr-84           Kr-86           Kr-82           Kr-83
     > 83.9115064D0,   85.910614D0,    81.913483D0,    82.914134D0,
*      Rb-85
     > 84.9117D0,      0.0,             0.0,             0.0,
*      Sr-88           Sr-84           Sr-86           Sr-87
     > 87.9056D0,      83.9134d0,      85.9094d0,      86.9089d0,
*      Y-89
     > 88.9054D0,      0.0,             0.0,             0.0,
*      Zr-90           Zr-91           Zr-92           Zr-94
     > 89.9043D0,      90.9053D0,      91.9046D0,      93.9061D0,
*      Nb-93
     > 92.9060D0,      0.0,             0.0,             0.0,
*      Mo-98           Mo-92           Mo-95           Mo-96
     > 97.9055D0,      91.9063D0,      94.90584D0,     95.9046D0,
*      Tc
     > 98.0D0,         0.0,             0.0,             0.0,
*      Ru-102          Ru-99           Ru-100          Ru-104
     > 101.9037D0,     98.9061D0,      99.9030D0,      103.9055D0,
*      Rh-103
     > 102.9048D0,     0.0,             0.0,             0.0,
*      Pd-106          Pd-104           Pd-105         Pd-108
     > 105.9032D0,     103.9036D0,      104.9046D0,    107.90389D0,
*      Ag-107          Ag-109
     > 106.90509d0,    108.9047D0,      0.0,             0.0,
*      Cd-114          Cd-110           Cd-111         Cd-112
     > 113.9036D0,     109.9030D0,      110.9042D0,    111.9028D0,
*      In-115          In-113
     > 114.9041D0,     112.9043D0,      0.0,             0.0,
*      Sn-118          Sn-116           Sn-117         Sn-119
     > 117.9018D0,     115.9021D0,      116.9031D0,    118.9034D0,
*      Sb-121          Sb-123
     > 120.9038D0,     122.9041D0,      0.0,             0.0,
*      Te-130          Te-125           Te-126         Te-128
     > 129.9067D0,     124.9044D0,      125.9032D0,    127.9047D0,
*      I-127
     > 126.9004D0,     0.0,             0.0,             0.0,
*      Xe-132          Xe-129           Xe-131         Xe-134
     > 131.9042D0,     128.9048D0,      130.9051D0,    133.9054D0/
*

       end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
