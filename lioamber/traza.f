!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!       INTERFACE REQUIERE PONERLO EN UN MODULO
!       INTERFACE traza
!         module procedure traza_R, traza_C
!       END INTERFACE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!       SUBROUTINE traza_R(M,Matrix,Traza)
       SUBROUTINE traza(M,Matrix,Suma)
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: M
       real*8,intent(in)  :: Matrix(M,M)
       real*8,intent(out) :: Suma
       integer :: kk
!------------------------------------------------------------------------------!
       Suma=0.0d0
       do kk=1,M
         Suma=Suma+Matrix(kk,kk)
       enddo
!------------------------------------------------------------------------------!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!       SUBROUTINE traza_C(M,Matrix,Traza)
!------------------------------------------------------------------------------!
!       implicit none
!       integer,intent(in)    :: M
!       complex*8,intent(in)  :: Matrix(M,M)
!       complex*8,intent(out) :: Traza
!       integer :: kk
!------------------------------------------------------------------------------!
!       Traza=0.0d0
!       do kk=1,M
!         Traza=Traza+Matrix(kk,kk)
!       enddo
!------------------------------------------------------------------------------!
!       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
