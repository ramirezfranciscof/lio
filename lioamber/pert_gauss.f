!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE pert_gauss(Energy)
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod,only:Fx,Fy,Fz,a0,natom,Iz,NCO,NUNP,epsilon,nstep
       implicit none
       real*8,intent(out) :: Energy

       real*8  :: dipx,dipy,dipz
       real*8  :: fldx,fldy,fldz
       real*8  :: Qc2,factor,g,term1,term2
       real*8  :: gausscurve,alpha
       integer :: mstep,kk
!
!
! SETUP
!----------------------------------------------------------!
       Fx=0.05d0
       Fy=0.0d0
       Fz=0.0d0
       g=1.0d0
       factor=2.54d0
       alpha=0.2
       mstep=50
       Qc2=0.0d0
       do kk=1,natom
         Qc2=Qc2+Iz(kk)
       enddo
       Qc2=Qc2-2*NCO+Nunp
       Qc2=Qc2**2
!
!
! GIVE GAUSSIAN SHAPE AND MODIFY FOCK
!----------------------------------------------------------!
       dipx=0.0d0
       dipy=0.0d0
       dipz=0.0d0
       call dip(dipx,dipy,dipz)
       gausscurve=exp(-alpha*(real(nstep-mstep))**2)
       fldx=fx*gausscurve
       fldy=fy*gausscurve
       fldz=fz*gausscurve
       write(*,*) fldx,fldy,fldz
!
!
! MODIFY FOCK AND CALCULATE ENERGY
!----------------------------------------------------------!
       call dip2(g,fldx,fldy,fldz)
       PRINT*
       PRINT*,Fx,Fy,Fz
       PRINT*,Fldx,Fldy,Fldz
       PRINT*,dipx,dipy,dipz
       PRINT*,g,factor,Qc2,a0,epsilon
       PRINT*
       term1=((-1.0d0)*(Fldx*dipx+Fldy*dipy+Fldz*dipz)*g)/factor
       term2=((-0.5d0)*(1.0d0-1.0d0/epsilon)*Qc2)/a0
       Energy=term1+term2
!
!
! FORMATS
!----------------------------------------------------------!
 200   FORMAT('Field shape: ',3(3X,E15.7))
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
