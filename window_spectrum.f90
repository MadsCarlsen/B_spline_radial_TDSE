program window_spectrum 
use utils 
use B_spline 
use GaussLegendre
use linalg_routines
use potential 
use spherical_harmonics
implicit none 

real(dp) :: r_max
integer :: N_points, spline_order, N_gauss_leg, l_max, para_index, N_splines, N_splines_tot, N_E, N_theta
real(dp), allocatable :: knot_list(:), r_list(:), S_mat(:,:), H0_mat(:,:), H1_mat(:,:)
real(dp), allocatable :: S_mat_lin(:,:), H0_mat_lin(:,:), H1_mat_lin(:,:)
complex(dp), allocatable :: cal_mat(:,:), psi_prop(:,:), psi_cal(:), psi_cal1(:), psi_window(:,:)
real(dp), allocatable :: E_list(:), partial_spectra(:,:), theta_list(:), PMD(:,:)
integer :: i, j, k, li, lj, LDAB, N_E_pos, E_count
real(dp) :: gamma, E_min, E_max, theta_min, theta_max 
complex(dp) :: overlap_ij, sph_harm_prod
logical :: cal_PMD 
complex(dp) :: phase
character(len=100) :: save_file 

! TEST 
complex(dp), allocatable :: psi_test(:)


! SETUP RADIAL AND KNOT LISTS 
namelist /GRID/ r_max, N_points, spline_order, N_gauss_leg, l_max, para_index
namelist /WINDOW_SPEC/ N_E, E_min, E_max, gamma, N_theta, theta_min, theta_max, cal_PMD
open(file='settings.nml', unit=1)
read(nml=GRID, unit=1)
read(nml=WINDOW_SPEC, unit=1)
close(1)

theta_min = theta_min * pi 
theta_max = theta_max * pi 

allocate(r_list(N_points))
!call knot_squence_lin(knot_list, N_splines, r_list, N_points, spline_order, r_max)
call knot_squence_lin_para(knot_list, N_splines_tot, r_list, N_points, spline_order, r_max, para_index)
N_splines = N_splines_tot - 2  ! Set first and last B-spline to 0 to get correct boundary conditions 

! Build overlap and Hamilton matrices 
allocate(H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines), S_mat(spline_order, N_splines))
call build_matrices(H0_mat, H1_mat, S_mat, knot_list, r_list, N_points, N_splines, spline_order, N_gauss_leg)

! Extend symmetric band matrices to full band matrices, needed when solving system of eqns.
LDAB = 3*(spline_order-1)+1 
allocate(S_mat_lin(LDAB, N_splines), H0_mat_lin(LDAB, N_splines), H1_mat_lin(LDAB, N_splines))

S_mat_lin = 0._dp 
H0_mat_lin = 0._dp 
H1_mat_lin = 0._dp
S_mat_lin(spline_order:2*spline_order-1, :) = S_mat
H0_mat_lin(spline_order:2*spline_order-1, :) = H0_mat
H1_mat_lin(spline_order:2*spline_order-1, :) = H1_mat

do i = 1, spline_order-1  ! Lower triangle 
    do j = 1, N_splines - i 
        S_mat_lin(2*spline_order-1 + i, j) = S_mat_lin(2*spline_order-1 - i, j+i)
        H1_mat_lin(2*spline_order-1 + i, j) = H1_mat_lin(2*spline_order-1 - i, j+i)
        H0_mat_lin(2*spline_order-1 + i, j) = H0_mat_lin(2*spline_order-1 - i, j+i)
    end do 
end do 

! Load propagated state 
allocate(psi_prop(N_splines, l_max+1))
open(1, file='propagation_data/WF_propagated.dat', form='unformatted')
read(1) psi_prop
close(1)

! Calculation of partial spectra! 
N_E = int((E_max - E_min) / (2._dp*gamma))+1
allocate(E_list(N_E), partial_spectra(l_max+1, N_E), cal_mat(LDAB, N_splines), psi_cal(N_splines), psi_cal1(N_splines))
allocate(psi_window(N_splines,l_max+1), theta_list(N_theta))

N_E_pos = 0
do i=1, N_E 
    E_list(i) = E_min + 2._dp*gamma * (i-1._dp) 
    if (E_list(i) > 0._dp) then 
        N_E_pos = N_E_pos + 1
    end if 
end do 
call linspace(theta_list, theta_min, theta_max, N_theta)

if (cal_PMD) then 
    allocate(PMD(N_theta, N_E_pos))
end if 

E_count = 1 
do j=1, N_E
    if (mod(j, 10) == 0) then 
        write(*,*) 'Calculating for energy ', E_list(j), '  ;  ', j,'/', N_E
    end if 

    do i=1, l_max+1
        ! Solve for first operator 
        phase = sqrt(cmplx(0._dp, 1._dp, dp))
        cal_mat = H0_mat_lin + (i-1._dp)*i/2._dp*H1_mat_lin - S_mat_lin*(E_list(j) - phase*gamma)
        call band_sym_matrix_vec_prod(S_mat, psi_prop(:,i), psi_cal, N_splines, spline_order)
        call solve_complex_lin_system(cal_mat, psi_cal, LDAB, N_splines, spline_order)

        ! Solve for second operator 
        cal_mat = H0_mat_lin + (i-1._dp)*i/2._dp*H1_mat_lin - S_mat_lin*(E_list(j) + phase*gamma)
        !psi_cal1 = psi_cal
        call band_sym_matrix_vec_prod(S_mat, psi_cal, psi_cal1, N_splines, spline_order)  ! DO I NEED THIS SECOND PRODUCT!?
        call solve_complex_lin_system(cal_mat, psi_cal1, LDAB, N_splines, spline_order)

        ! Save the coeff. for the energy-projected state, used for either total spectrum or PMD
        psi_window(:,i) = psi_cal1
    end do
  
    ! Now calculate total spectrum and potentially the PMD 
    do i=1, l_max+1
        call band_sym_matrix_vec_prod(S_mat, psi_window(:,i), psi_cal, N_splines, spline_order)
        partial_spectra(i, j) = real(sum(conjg(psi_window(:,i)) * psi_cal)) 
    end do

    ! Calculate PMD if wanted, and only for positive energies 
    if (cal_PMD .and. E_list(j) > 0._dp) then 
        do li=1, l_max+1
            do lj=1, li
                ! Calculate overlap integral 
                call band_sym_matrix_vec_prod(S_mat, psi_window(:,li), psi_cal, N_splines, spline_order)
                overlap_ij = sum(conjg(psi_window(:,lj)) * psi_cal)

                ! Save for the different angles 
                do k=1, N_theta
                    sph_harm_prod = sph_harm(cos(theta_list(k)), 0._dp, li-1, 0) * sph_harm(cos(theta_list(k)), 0._dp, lj-1, 0)

                    if (li == lj) then 
                        PMD(k,E_count) = PMD(k,E_count) + real(sph_harm_prod * overlap_ij) 
                    else 
                        PMD(k,E_count) = PMD(k,E_count) + real(sph_harm_prod * (overlap_ij + conjg(overlap_ij)))
                    end if 
                end do
            end do
        end do
        E_count = E_count + 1 
    end if 

end do

! Save the partial spectra 
do i=1, l_max+1 
    write(save_file, "(A,I0,A)") 'spectrum_data/window_partial_spectrum_', i-1, '.txt'
    open(1, file=save_file)
    do j=1, N_E
        write(1,*) partial_spectra(i,j)
    end do
    close(1)
end do 

open(1, file='spectrum_data/window_energies.txt')
do i=1, N_E
    write(1,*) E_list(i)
end do 

! Save PMD if calculated 
if (cal_PMD) then 
    open(1, file='spectrum_data/PMD.dat', form='unformatted')
    write(1) PMD
    close(1)

    open(2, file='spectrum_data/theta_list.txt')
    do i=1, N_theta
        write(2,*) theta_list(i)
    end do 
end if 

contains 

subroutine band_sym_complex_matrix_vec_prod(matrix, vec, res, N_vec, N_band)
    integer, intent(in) :: N_vec, N_band  ! N_vec is size of diagnoal, N_band is nr. of superdiagnoals + diagonal
    complex(dp), intent(in) :: matrix(N_band, N_vec) 
    complex(dp), intent(in) :: vec(N_vec)
    complex(dp), intent(out) :: res(N_vec)
    integer :: i, j, k, Nj
    complex(dp) :: sum_i 
    
    k = N_band - 1  ! Nr. of superdiagonals in band 
    res = 0._dp  ! Make sure result is initially 0

    do j=1, N_vec  ! Loop over the columns first, as fast index in Fortran is rows 
        sum_i = 0._dp  ! Sum to contain contributions to res(j), from lower triangle of mat 
        Nj = N_band-j
        do i=max(1,j-k), j-1  ! Loop down over rows, not including diagonal 
            res(i) = res(i) + matrix(Nj+i, j) * vec(j)  ! Contribution to res(i) from upper triangle 
            sum_i = sum_i + matrix(Nj+i, j) * vec(i)  ! Parts for lower triagnle 
        end do 
        res(j) = res(j) + matrix(N_band, j)*vec(j) + sum_i  ! Add diagnoal + terms from lower triangle
    end do
end subroutine band_sym_complex_matrix_vec_prod 


subroutine build_matrices(H0_mat, H1_mat, S_mat, knot_list, r_list, N_points, N_splines, spline_order, N_gauss_leg)
    integer, intent(in) :: N_points, N_splines, spline_order, N_gauss_leg
    real(dp), intent(out) :: H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines), S_mat(spline_order, N_splines)
    real(dp), intent(in) :: r_list(N_points)
    real(dp), intent(in) :: knot_list(:)

    integer :: N_splines_tot, knot_init ! Total nr. of B splines on the interval 
    real(dp) :: H0_temp(spline_order, N_splines+2), H1_temp(spline_order, N_splines+2), S_temp(spline_order, N_splines+2)
    real(dp) :: spline_arr(spline_order, N_gauss_leg), deriv_arr(spline_order, N_gauss_leg)
    real(dp) :: x_quad(N_gauss_leg), x_arr(N_gauss_leg), w_quad(N_gauss_leg), w_arr(N_gauss_leg), V_arr(N_gauss_leg)
    real(dp) :: H0_ij, H1_ij, S_ij 
    integer :: i, j, n, m, k
    

    N_splines_tot = N_splines + 2
    H0_temp = 0._dp 
    H1_temp = 0._dp 
    S_temp = 0._dp 

    ! Determine rightmost index in knot sequence corresponding to r=0
    knot_init = find_index(0._dp, knot_list)

    ! Get Gauss-Legendre quad points and weights 
    call gauss_leg(N_gauss_leg, x_quad, w_quad)

    ! Loop over the different radial intervals 
    do k=1, N_points-1
        ! First rescale quad points to current interval 
        call rescale_quad(x_arr, w_arr, x_quad, w_quad, r_list(k), r_list(k+1))

        ! Evaluate B-splines and potential
        call eval_B_splines_arr(spline_arr, deriv_arr, x_arr, knot_list, knot_init, spline_order)
        V_arr = V_pot(x_arr)

        ! Now we need to evaluate the integral for all unique combinations of splines 
        do m=1, spline_order 
            do n=m, spline_order
                H1_ij = sum(spline_arr(m,:)*spline_arr(n,:)/x_arr**2 * w_arr)  ! Centrifugal term
                H0_ij =  sum(spline_arr(m,:)*spline_arr(n,:)*V_arr * w_arr)  ! Potential term 
                H0_ij = H0_ij + 0.5_dp * sum(deriv_arr(n,:) * deriv_arr(m,:) * w_arr)  ! Kinetic term 
                S_ij =  sum(spline_arr(m,:) * spline_arr(n,:) * w_arr)
        
                ! Now save in correct places 
                j = n + k-1  ! Index juggeling to save in band storage scheme, upper triangle, column in column, diagnoal in bottom
                i = spline_order + m - n 
                H0_temp(i,j) = H0_temp(i,j) + H0_ij
                H1_temp(i,j) = H1_temp(i,j) + H1_ij
                S_temp(i,j) = S_temp(i,j) + S_ij
            end do  
        end do

        knot_init = knot_init + 1  ! Update when moving along on radial grid 
    end do 

    ! Save results, excluding first and last spline 
    H0_mat = H0_temp(:, 2:N_splines_tot-1)
    H1_mat = H1_temp(:, 2:N_splines_tot-1)
    S_mat = S_temp(:, 2:N_splines_tot-1)
    do i=1, spline_order-1
        H0_mat(spline_order-i, i) = 0._dp 
        H1_mat(spline_order-i, i) = 0._dp 
        S_mat(spline_order-i, i) = 0._dp 
    end do
end subroutine build_matrices

end program window_spectrum 