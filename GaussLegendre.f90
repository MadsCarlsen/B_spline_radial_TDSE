module GaussLegendre

use utils 
implicit none 
private 

! Define working precision for the module 
integer, parameter :: wp = dp 

! Define public routines 
public gauss_leg, rescale_quad

contains 

subroutine gauss_leg(N, x_arr, w_arr, lower, upper)
    ! Get the Gauss-Legendre quad. points/weights for N'th order.
    integer, intent(in) :: N 
    real(wp), intent(out) :: x_arr(N), w_arr(N)
    integer :: i, j, k, N_iter  
    real(wp) :: xi, x_new, p0, p1, p2, pd 
    real(wp) :: eps = 3e-14_wp  ! Precision to calculate the roots
    integer :: max_iter = 100  ! Max iterations in the Newton root finder - should not be reached...
    real(wp), intent(in), optional :: lower, upper  ! Limits to rescale the quad points/weights to
    real(wp) :: a,b  ! Limits to be used if lower/upper passed

    ! Check if N is odd or even to determine number of roots to calculate 
    if (mod(N,2) == 0 ) then 
        N_iter = int(N/2) 
    else 
        N_iter = int(N/2)+1  ! To also capture the root at 0 
    end if 

    ! Obtain the roots for the N'th order legendre polynomial using Newtons's method 
    ! Only have to find the positive roots, since the polynomials are even/odd?
    do i = 1, N_iter
        ! Get guess for the i'th root (starts from 1 and goes down)
        xi = cos(pi * (4_wp*i-1_wp) / (4_wp*N + 2_wp))  ! Francesco Tricomi asymptotic formula (crazy!)

        ! Use Newton's method to refine the root 
        do j=1, max_iter 
            ! First obtain the legendre polynomial and its derivative at xi using recursion (Arfken p. 719)
            p2 = 0 
            p1 = 1_wp 
            do k = 1, N
                p0 = ((2_wp*k - 1_wp)*xi*p1 - (k-1_wp)*p2) / k
                p2 = p1 
                p1 = p0 
            end do 
            pd = N*(p2-xi*p0) / (1_wp - xi*xi)  ! The derivative in xi

            ! Then take a Newton step 
            x_new = xi - p0/pd 

            ! Check for convergence 
            if (abs(x_new - xi) < eps) exit
            xi = x_new 
        end do 

        ! Save the root and calculate the weight 
        x_arr(i) = x_new 
        w_arr(i) = 2_wp / ((1_wp - x_new*x_new) * pd*pd)
    end do 
    
    ! Add the remaining roots and weights by symmetry
    do i = N_iter+1, N 
        x_arr(i) = -x_arr(i-N_iter) 
        w_arr(i) = w_arr(i-N_iter) 
    end do

    ! Rescale if passed 
    if (present(lower) .or. present(upper)) then
        ! Set standard value of upper/lower limit if not present 
        if (present(lower)) then 
            a = lower 
        else 
            a = -1_wp 
        end if 

        if (present(upper)) then 
            b = upper 
        else 
            b = 1_wp 
        end if 
        
        ! Now rescale 
        do i=1, size(x_arr)
            x_arr(i) = a + (b-a)/2_wp * (x_arr(i)+1_wp)
            w_arr(i) = (b-a)/2_wp * w_arr(i)
        end do 
    end if 
end subroutine gauss_leg

! Subroutine to rescale quad points and weights to other interval
subroutine rescale_quad(x_scale, w_scale, x_arr, w_arr, lower, upper)
    real(wp), intent(in) :: x_arr(:), w_arr(:), lower, upper
    real(wp), intent(out) :: x_scale(size(x_arr)), w_scale(size(x_arr))
    integer :: i 

    do i=1, size(x_arr)
        x_scale(i) = lower + (upper-lower)/2._wp * (x_arr(i)+1._wp)
        w_scale(i) = (upper-lower)/2._wp * w_arr(i)
    end do
end subroutine rescale_quad 

end module GaussLegendre 
