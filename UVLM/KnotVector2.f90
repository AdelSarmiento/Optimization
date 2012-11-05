subroutine KnotVector2(E,p,C,Ui,Uf,periodic,m,U)
  implicit none
  integer(kind=4), intent(in)  :: E,p,C
  DOUBLE PRECISION, intent(in)  :: Ui,Uf
  logical(kind=4), intent(in)  :: periodic
  integer(kind=4), intent(in)  :: m
  DOUBLE PRECISION, intent(out) :: U(0:m)
  integer(kind=4) :: i, j, k
  DOUBLE PRECISION :: dU
  dU = (Uf-Ui)/E
  k = p+1
  !WRITE(*,*) E,p,C,Ui,Uf,periodic,m,U
  !WRITE(*,*) k
  do i = 1, (E-1)
  
     do j = 1, (p-C)
	    
        !U(k) = Ui + i*dU
		
        k = k+1
     end do
  end do
  
  do k = 0, p
     !U(k)   = Ui
     !U(m-k) = Uf
  end do
  
  if (periodic) then
     do k = 0, C
        !U(k)     = Ui - Uf + U(m-C-(p+1)+k)
        !U(m-C+k) = Uf - Ui + U(p+1+k)
     end do
  end if
  
  !WRITE(*,*) k
end subroutine KnotVector2

