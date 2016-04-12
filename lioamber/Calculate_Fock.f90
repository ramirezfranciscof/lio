!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine Calculate_Fock(DensMao,FockMao,Energy)
  use maskrmm
  use garcha_mod, only:M,RMM,kkind,kkinds,cool,cools,igrid2 & 
                      ,nuc,natom,natomc,d,r,atmin,rmax,jatc &
                      ,nnps,nnpp,nnpd,nshell,Md

  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  real*8,intent(out)    :: FockMao(M,M)
  real*8,intent(out)    :: Energy

  real*8   :: Energy_Nuclear,Energy_1e,Energy_Coulomb
  real*8   :: Energy_SolvT,Energy_SolvF
  real*8   :: Energy_Exchange
  integer  :: kk,ii,jj,idx,idx0

  real*8   :: alf,rexp
  integer  :: igpu,zij,ti,tj
  integer  :: MM,MMd
  logical  :: MEMO


! Store DensMao in RMM
!------------------------------------------------------------------------------!
  call rmmput_dens(DensMao)


! Initializations: Copied from SCF
!------------------------------------------------------------------------------!
  do ii=1,natom
    natomc(ii)=0
    do jj=1,natom
      d(ii,jj)=0.0d0
      d(ii,jj)=d(ii,jj)+(r(ii,1)-r(jj,1))**2
      d(ii,jj)=d(ii,jj)+(r(ii,2)-r(jj,2))**2
      d(ii,jj)=d(ii,jj)+(r(ii,3)-r(jj,3))**2
      zij=atmin(ii)+atmin(jj)
      ti=atmin(ii)/zij
      tj=atmin(jj)/zij
      alf=atmin(ii)*tj
      rexp=alf*d(ii,jj)
      if (rexp.lt.rmax) then
        natomc(ii)=natomc(ii)+1
        jatc(natomc(ii),ii)=jj
      endif
    enddo
  enddo

  do ii=nshell(0),1,-1
    nnps(nuc(ii))=ii
  enddo

  do ii=nshell(0)+nshell(1),nshell(0)+1,-1
    nnpp(nuc(ii))=ii
  enddo

  do ii=M,nshell(0)+nshell(1)+1,-1
    nnpd(nuc(ii))=ii
  enddo



! Calculate S-part of fock
!------------------------------------------------------------------------------!
  call g2g_timer_sum_start('Exchange-correlation grid setup')
  call g2g_reload_atom_positions(igrid2)
  call g2g_timer_sum_stop('Exchange-correlation grid setup')

  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(Energy_Nuclear)



! Calculate fixed-parts of fock
!------------------------------------------------------------------------------!
  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

  call g2g_reload_atom_positions(igrid2)
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()

  if (igpu.le.1) then
    call intsol(Energy_SolvF,Energy_SolvT,.true.)
  else
    call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
  endif

  call int2()

  if (igpu.gt.2) then
    call aint_coulomb_init()
  endif
  if (igpu.eq.5) MEMO = .false.
  if (MEMO) then
    call g2g_timer_start('int3mem')
    call g2g_timer_sum_start('Coulomb precalc')
    call int3mem()
    call g2g_timer_stop('int3mem')
    call g2g_timer_sum_stop('Coulomb precalc')
  endif


! Calculate unfixed Fock in RMM
!------------------------------------------------------------------------------!
  call int3lu(Energy_Coulomb)
  call g2g_solve_groups(0,Energy_Exchange,0)
  MM=M*(M+1)/2
  MMd=Md*(Md+1)/2
  idx0=3*MM+2*MMd
  Energy_1e=0.0d0
  do kk=1,MM
    Energy_1e=Energy_1e+RMM(kk)*RMM(idx0+kk)
  enddo

  Energy=0.0d0
  Energy=Energy+Energy_Nuclear
  Energy=Energy+Energy_1e
  Energy=Energy+Energy_Coulomb
  Energy=Energy+Energy_SolvT
  Energy=Energy+Energy_Exchange


! Extract FockMao from RMM
!------------------------------------------------------------------------------!
  call rmmget_fock(FockMao)

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
