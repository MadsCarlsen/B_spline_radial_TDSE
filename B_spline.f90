module B_spline 
use utils 
implicit none 
private 


public :: knot_squence_lin, eval_B_splines, find_index, eval_B_splines_arr, eval_function_B_basis, &
          eval_complex_function_B_basis, knot_squence_lin_para

contains 

! Function to find the rightmost knot index that is less than x 
function find_index(x, knot_list) result(index)
    real(dp), intent(in) :: x, knot_list(:)
    integer :: index
    integer :: i 

    i = 1 
    do while (knot_list(i) <= x)
        i = i + 1
    end do 
    index = i-1  ! There will always be added at least one 
end function find_index 

! Subroutine to provide a linear knot sequence, setting first and last B-spline to 0!
subroutine knot_squence_lin(knot_list, N_splines, r_list, N_points, spline_order, r_max)
    real(dp), allocatable, intent(out) :: knot_list(:)
    integer, intent(out) :: N_splines
    integer, intent(in) :: N_points, spline_order 
    real(dp), intent(in) :: r_max 
    integer :: N_k, i  ! N_k is number of knots
    real(dp), intent(out) :: r_list(N_points)

    N_k = N_points - 2 + 2*spline_order
    N_splines = N_points - 2 + spline_order 
    allocate(knot_list(N_k))

    ! Create the linear break-point list 
    call linspace(r_list, 0._dp, r_max, N_points)

    knot_list = -1._dp 
    ! First spline_order-1 points are leftmost endpoint (r=0)
    knot_list(1:spline_order) = r_list(1)
    ! Central points
    do i=2, N_points-1
        knot_list(i + spline_order-1) = r_list(i)
    end do 
    ! Last spline_order-1 points are rightmost point 
    knot_list(N_k - spline_order + 1: N_k) = r_list(N_points)
end subroutine knot_squence_lin

! Subroutine to provide a parabolic-linear, setting first and last B-spline to 0!
subroutine knot_squence_lin_para(knot_list, N_splines, r_list, N_points, spline_order, r_max, i0)
    real(dp), allocatable, intent(out) :: knot_list(:)
    integer, intent(out) :: N_splines
    integer, intent(in) :: N_points, spline_order, i0
    real(dp), intent(in) :: r_max 
    integer :: N_k, i  ! N_k is number of knots
    real(dp), intent(out) :: r_list(N_points)

    real(dp) :: r0, alpha, beta 

    N_k = N_points - 2 + 2*spline_order
    N_splines = N_points - 2 + spline_order 
    allocate(knot_list(N_k))

    ! Create break-point squence that is parabolic at center and linear after index i0
    r0 = r_max*(i0-1._dp) / (2._dp*N_points - i0 - 1._dp)
    alpha = r0 / (i0 - 1._dp)**2 
    beta = (r_max - r0) / (1._dp*N_points - i0)
    do i=1, N_points 
        if (i < i0) then 
            r_list(i) = alpha * (i-1._dp)**2 
        else 
            r_list(i) = r0 + beta * (i-1._dp*i0)
        end if 
    end do

    knot_list = -1._dp 
    ! First spline_order-1 points are leftmost endpoint (r=0)
    knot_list(1:spline_order) = r_list(1)
    ! Central points
    do i=2, N_points-1
        knot_list(i + spline_order-1) = r_list(i)
    end do 
    ! Last spline_order-1 points are rightmost point 
    knot_list(N_k - spline_order + 1: N_k) = r_list(N_points)
end subroutine knot_squence_lin_para


! Subroutine to return the spline_order number of splines on interval [t_i, t_(i+1)]
! de Boor algorithm, 'On Calculating with B-splines' p. 58.
subroutine eval_B_splines(spline_res, spline_deriv, x, knot_list, knot_init, spline_order)
    real(dp), intent(in) :: x, knot_list(:)
    integer, intent(in) :: knot_init, spline_order 
    real(dp), intent(out) :: spline_res(spline_order), spline_deriv(spline_order) ! Each spline evaluated at x_list points 

    integer :: s, r 
    real(dp) :: ta, tb, M
    real(dp) :: spline_temp(spline_order) 

    ! Outer loop over order 'column'
    spline_res = 0._dp 
    spline_temp = 0._dp 
    spline_temp(1) = 1._dp 
    do s=1, spline_order-1 
        spline_res(1) = 0._dp 

        ! Inner loop over the B-spline index 
        do r = 1, s 
            ta = knot_list(knot_init + r) - x
            tb = x - knot_list(knot_init + r - s)

            M = spline_temp(r) / (ta + tb )
            spline_res(r) = spline_res(r) + ta*M
            spline_res(r+1) = tb * M  
        end do 
        
        if (s < spline_order-1) then  ! Don't update it the last time - we use it to find derivatives 
            spline_temp = spline_res
        end if 
    end do

    ! Now calculate the derivatives of the splines 
    spline_deriv(1) = 0._dp 
    s = spline_order-1
    do r = 1, spline_order-1 
        ta = knot_list(knot_init + r) - knot_list(knot_init + r - s)
        M = spline_temp(r) / ta 
        spline_deriv(r) = spline_deriv(r) - s*M 
        spline_deriv(r+1) = s*M  
    end do
end subroutine eval_B_splines

! Evaluate B-splines on serval points. NB. entire x_list must lie within same knot interval! 
subroutine eval_B_splines_arr(spline_res, spline_deriv, x_list, knot_list, knot_init, spline_order)
    real(dp), intent(in) :: x_list(:), knot_list(:)
    integer, intent(in) :: knot_init, spline_order 
    real(dp), intent(out) :: spline_res(spline_order, size(x_list)), spline_deriv(spline_order, size(x_list)) ! Each spline evaluated at x_list points 
    real(dp) :: spline_single(spline_order), deriv_single(spline_order)
    integer :: i 

    do i=1, size(x_list) 
        ! Eval B-splines 
        call eval_B_splines(spline_single, deriv_single, x_list(i), knot_list, knot_init, spline_order)
        spline_res(:,i) = spline_single 
        spline_deriv(:,i) = deriv_single
    end do
end subroutine eval_B_splines_arr 

! Evaluate a function expanded in a B-spline basis for points in x_list.
! Assumes that the first and last B-spline is not included in the expansion! 
subroutine eval_function_B_basis(func_res, coeff_list, x_list, knot_list, spline_order, N_splines)
    integer, intent(in) :: spline_order, N_splines
    real(dp), intent(in) :: x_list(:), knot_list(:), coeff_list(N_splines)
    real(dp), intent(out) :: func_res(size(x_list))
    
    real(dp) :: spline_res(spline_order), spline_deriv(spline_order), coeff_expand(N_splines+2), res_i
    integer :: i, j, knot_init

    ! First we append a 0 in the front and back of coeff_list. Done since states do not 
    ! include the first and last B-spline, but B-spline evaluation assumes all splines 
    coeff_expand(2:N_splines+1) = coeff_list 
    coeff_expand(1) = 0._dp 
    coeff_expand(N_splines+2) = 0._dp 

    ! Determine starting index 
    knot_init = find_index(x_list(1), knot_list)

    ! Loop over x_list, bulding the linear combination for each value
    do i=1, size(x_list)
        res_i = 0._dp 

        ! Check we are performing calculations on the correct knot interval
        if (x_list(i) > knot_list(knot_init+1)) then 
            knot_init = knot_init + 1 
        end if 

        ! Evaluate the B-splines in the current point 
        call eval_B_splines(spline_res, spline_deriv, x_list(i), knot_list, knot_init, spline_order)

        ! Sum up the linear expansion, making sure to multiply with correct coefficients 
        do j=1, spline_order
            res_i = res_i + coeff_expand(knot_init - spline_order + j) * spline_res(j)
        end do
        func_res(i) = res_i 
    end do
end subroutine eval_function_B_basis

! Evaluate a function expanded in a B-spline basis for points in x_list.
! Assumes that the first and last B-spline is not included in the expansion! 
subroutine eval_complex_function_B_basis(func_res, coeff_list, x_list, knot_list, spline_order, N_splines)
    integer, intent(in) :: spline_order, N_splines
    real(dp), intent(in) :: x_list(:), knot_list(:)
    complex(dp), intent(in) :: coeff_list(N_splines)
    complex(dp), intent(out) :: func_res(size(x_list))
    
    real(dp) :: spline_res(spline_order), spline_deriv(spline_order)
    complex(dp) :: coeff_expand(N_splines+2), res_i
    integer :: i, j, knot_init

    ! First we append a 0 in the front and back of coeff_list. Done since states do not 
    ! include the first and last B-spline, but B-spline evaluation assumes all splines 
    coeff_expand(2:N_splines+1) = coeff_list 
    coeff_expand(1) = 0._dp 
    coeff_expand(N_splines+2) = 0._dp 

    ! Determine starting index 
    knot_init = find_index(x_list(1), knot_list)

    ! Loop over x_list, bulding the linear combination for each value
    do i=1, size(x_list)
        res_i = 0._dp 

        ! Check we are performing calculations on the correct knot interval
        if (x_list(i) > knot_list(knot_init+1)) then 
            knot_init = knot_init + 1 
        end if 

        ! Evaluate the B-splines in the current point 
        call eval_B_splines(spline_res, spline_deriv, x_list(i), knot_list, knot_init, spline_order)

        ! Sum up the linear expansion, making sure to multiply with correct coefficients 
        do j=1, spline_order
            res_i = res_i + coeff_expand(knot_init - spline_order + j) * spline_res(j)
        end do
        func_res(i) = res_i 
    end do
end subroutine eval_complex_function_B_basis

end module B_spline 