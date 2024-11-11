program spectrum 
use utils 
use B_spline
use GaussLegendre 
implicit none 

integer :: N_points, para_index, spline_order, N_gauss_leg, l_max, N_splines, N_splines_tot
real(dp) :: r_max 
real(dp), allocatable :: knot_list(:), r_list(:), S_mat(:,:) 
complex(dp), allocatable :: psi_prop(:,:), projection, cal_state(:)
real(dp), allocatable :: eigh_cont_arr(:,:) 
complex(dp), allocatable :: spectrum_arr(:)
integer, allocatable :: nr_states_arr(:)
integer :: i, j
character(len=100) :: load_file, save_file 


! SETUP RADIAL AND KNOT LISTS 
namelist /GRID/ r_max, N_points, spline_order, N_gauss_leg, l_max, para_index
open(file='settings.nml', unit=1)
read(nml=GRID, unit=1)
close(1)

allocate(r_list(N_points))
!call knot_squence_lin(knot_list, N_splines, r_list, N_points, spline_order, r_max)
call knot_squence_lin_para(knot_list, N_splines_tot, r_list, N_points, spline_order, r_max, para_index)
N_splines = N_splines_tot - 2  ! Set first and last B-spline to 0 to get correct boundary conditions 


! DETERMINE OVERLAP MATRIX USED IN CALCULATIONS  
allocate(S_mat(spline_order, N_splines))
call build_S_mat(S_mat, knot_list, r_list, N_points, N_splines, spline_order, N_gauss_leg)


! LOAD THE PROPAGATED STATE 
allocate(psi_prop(N_splines, l_max+1))
open(1, file='propagation_data/WF_propagated.dat', form='unformatted')
read(1) psi_prop
close(1)

! CALCULATE PARTIAL ANGULAR INTEGRATED SPECTRUMS 
allocate(nr_states_arr(l_max+1))
open(1, file='eigenstate_data/nr_continuum_states.txt')
read(1,*) nr_states_arr 
close(1)

! Main loop calculating and saving spectrum for each angular momentum 
allocate(cal_state(N_splines))
do i=1, l_max+1    
    ! First load the continuum eigenstates 
    allocate(eigh_cont_arr(nr_states_arr(i), N_splines), spectrum_arr(nr_states_arr(i)))
    write(load_file, "(A,I0,A)") 'eigenstate_data/state_continuum_', i-1, '.dat'
    open(1, file=load_file, form='unformatted')
    read(1) eigh_cont_arr 
    close(1)

    ! Perform the projections 
    call band_sym_matrix_vec_prod(S_mat, psi_prop(:,i), cal_state, N_splines, spline_order)
    do j=1, nr_states_arr(i)
        spectrum_arr(j) = sum(eigh_cont_arr(j,:) * cal_state) ! eigh state should be cc. but is real 
    end do

    ! Save the partial spectrum 
    write(save_file, "(A,I0,A)") 'spectrum_data/partial_spectrum_', i-1, '.txt'
    open(1, file=save_file)
    do j=1, nr_states_arr(i)
        write(1,*) real(spectrum_arr(j)), aimag(spectrum_arr(j))
    end do
    close(1)
 
    ! Deallocate to prepare for new iteration 
    deallocate(eigh_cont_arr, spectrum_arr)
end do

contains 

subroutine build_S_mat(S_mat, knot_list, r_list, N_points, N_splines, spline_order, N_gauss_leg)
    integer, intent(in) :: N_points, N_splines, spline_order, N_gauss_leg
    real(dp), intent(out) :: S_mat(spline_order, N_splines)
    real(dp), intent(in) :: r_list(N_points)
    real(dp), intent(in) :: knot_list(:)

    integer :: N_splines_tot, knot_init ! Total nr. of B splines on the interval 
    real(dp) :: S_temp(spline_order, N_splines+2)
    real(dp) :: spline_arr(spline_order, N_gauss_leg), deriv_arr(spline_order, N_gauss_leg)
    real(dp) :: x_quad(N_gauss_leg), x_arr(N_gauss_leg), w_quad(N_gauss_leg), w_arr(N_gauss_leg), V_arr(N_gauss_leg)
    real(dp) :: S_ij 
    integer :: i, j, n, m, k
    

    N_splines_tot = N_splines + 2
    S_temp = 0._dp 

    ! Determine rightmost index in knot sequence corresponding to r=0
    knot_init = find_index(0._dp, knot_list)

    ! Get Gauss-Legendre quad points and weights 
    call gauss_leg(N_gauss_leg, x_quad, w_quad)

    ! Loop over the different radial intervals 
    do k=1, N_points-1
        ! First rescale quad points to current interval 
        call rescale_quad(x_arr, w_arr, x_quad, w_quad, r_list(k), r_list(k+1))

        ! Evaluate B-splines
        call eval_B_splines_arr(spline_arr, deriv_arr, x_arr, knot_list, knot_init, spline_order)
        
        ! Now we need to evaluate the integral for all unique combinations of splines 
        do m=1, spline_order 
            do n=m, spline_order
                S_ij =  sum(spline_arr(m,:) * spline_arr(n,:) * w_arr)
        
                ! Now save in correct places 
                j = n + k-1  ! Index juggeling to save in band storage scheme, upper triangle, column in column, diagnoal in bottom
                i = spline_order + m - n 
                S_temp(i,j) = S_temp(i,j) + S_ij
            end do  
        end do

        knot_init = knot_init + 1  ! Update when moving along on radial grid 
    end do 

    ! Save results, excluding first and last spline 
    S_mat = S_temp(:, 2:N_splines_tot-1)
    do i=1, spline_order-1
        S_mat(spline_order-i, i) = 0._dp 
    end do
end subroutine build_S_mat

! Routine for calculating the matrix vector product, when the matrix is band symmetric
subroutine band_sym_matrix_vec_prod(matrix, vec, res, N_vec, N_band)
    integer, intent(in) :: N_vec, N_band  ! N_vec is size of diagnoal, N_band is nr. of superdiagnoals + diagonal
    real(dp), intent(in) :: matrix(N_band, N_vec) 
    complex(dp), intent(in) :: vec(N_vec)
    complex(dp), intent(out) :: res(N_vec)
    integer :: i, j, k, Nj
    complex(dp) :: sum_i 
    
    k = N_band - 1  ! Nr. of superdiagonals in band 
    res = 0._dp  ! Make sure result is initially 0, as calculating requr

    do j=1, N_vec  ! Loop over the columns first, as fast index in Fortran is rows 
        sum_i = 0._dp  ! Sum to contain contributions to res(j), from lower triangle of mat 
        Nj = N_band-j
        do i=max(1,j-k), j-1  ! Loop down over rows, not including diagonal 
            res(i) = res(i) + matrix(Nj+i, j) * vec(j)  ! Contribution to res(i) from upper triangle 
            sum_i = sum_i + matrix(Nj+i, j) * vec(i)  ! Parts for lower triagnle 
        end do 
        res(j) = res(j) + matrix(N_band, j)*vec(j) + sum_i  ! Add diagnoal + terms from lower triangle
    end do
end subroutine band_sym_matrix_vec_prod 


end program spectrum 