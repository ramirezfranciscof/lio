!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehren_verlet_n(dof,dt,mass,force,oldpos,nowpos,newpos,nowvel)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  integer,intent(in)     :: dof ! Degrees of Freedom
  real*8,intent(in)      :: dt  ! Time step

  real*8,intent(in)      :: mass(dof)
  real*8,intent(in)      :: force(dof,3)
  real*8,intent(in)      :: oldpos(dof,3)
  real*8,intent(in)      :: nowpos(dof,3)
  real*8,intent(out)     :: newpos(dof,3)
  real*8,intent(out)     :: nowvel(dof,3)

  integer :: nn,kk
  real*8  :: dt2,dtsq


  dtsq=dt*dt
  dt2=2.0d0*dt

  do kk=1,3
  do nn=1,dof
    newpos(nn,kk)=2.0d0*nowpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    newpos(nn,kk)=newpos(nn,kk)+(dtsq*force(nn,kk))/mass(nn)

    nowvel(nn,kk)=newpos(nn,kk)-oldpos(nn,kk)
    nowvel(nn,kk)=nowvel(nn,kk)/dt2
  enddo
  enddo

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
