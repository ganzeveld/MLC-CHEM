      subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none

!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
 
        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i
 
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
 	
    end subroutine solve_tridiag
	
	