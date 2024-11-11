program propagate 
use utils 
use B_spline 
use GaussLegendre
use potential 
use linalg_routines
use omp_lib
implicit none 

real(dp) :: r_max, tol_residual, res_bound_low, res_bound_upper
integer :: N_splines_tot, N_splines, N_points, spline_order, N_knots, N_gauss_leg, knot_init, para_index, max_krylov
real(dp), allocatable :: knot_list(:), r_list(:)

real(dp), allocatable :: x_quad(:), w_quad(:), x_arr(:), w_arr(:), V_arr(:), psi0_real(:)
real(dp), allocatable :: spline_arr(:,:), deriv_arr(:,:)
complex(dp), allocatable :: psi0(:)

real(dp), allocatable :: H0_cal(:,:), H1_cal(:,:), D0_cal(:,:), D1_cal(:,:), S_cal(:,:)  ! Temp matrices 
real(dp) :: h0_ij, h1_ij, d0_ij, d1_ij, s_ij
real(dp), allocatable :: S_mat(:,:), H0_mat(:,:), H1_mat(:,:), D0_mat(:,:), D1_mat(:,:)  ! Mats without first/last spline
complex(dp), allocatable :: S_cholesky(:,:)  ! Complex since we use it to solve complex Sx=b with LAPACK
complex(dp), allocatable :: LU_arr(:,:,:)
integer, allocatable :: pivot_arr(:,:)
integer :: i, j, k, n, m

! Varibales for propagation 
integer :: l_max, Nt, N_krylov 
real(dp) :: dt, t_sim, time, dt_taken, dt_remain, residual
real(dp), allocatable :: alpha_list(:)
complex(dp), allocatable :: psi(:,:)

! Field params 
real(dp) :: lambda, I0, CEP, omega, E_max, Up, rtUp, laser_duration 
integer :: Nc

! TEST 
complex(dp), allocatable :: l_norms(:)
complex(dp), allocatable :: test(:), test_res(:), test_arr(:,:), test_res_arr(:,:)
complex(dp), allocatable :: norm, func_eval_arr(:)
integer :: N_steps
complex(dp), allocatable :: psi_init(:), psi_cal(:)  
real(dp), allocatable :: population(:,:), time_list(:)

! LOAD SIMULATION DATA 

! Field params 
namelist /GRID/ r_max, N_points, spline_order, N_gauss_leg, l_max, para_index
namelist /PROPAGATION/ lambda, I0, Nc, CEP, dt, t_sim, N_krylov, tol_residual, res_bound_low, res_bound_upper
open(file='settings.nml', unit=1)
read(nml=GRID, unit=1)
read(nml=PROPAGATION, unit=1)
close(1)

! Grid and B spline params 
N_gauss_leg = spline_order + N_gauss_leg 
allocate(r_list(N_points))
!call knot_squence_lin(knot_list, N_splines_tot, r_list, N_points, spline_order, r_max)
call knot_squence_lin_para(knot_list, N_splines_tot, r_list, N_points, spline_order, r_max, para_index)

N_splines = N_splines_tot - 2  ! Set first and last B-spline to 0 to get correct boundary conditions 

! LOAD THE INITIAL STATE 
allocate(psi0_real(N_splines))
open(1, file='eigenstate_data/ground_state.dat', form='unformatted')
read(1) psi0_real  
close(1)

! SET NR OF THREADS 
!call omp_set_num_threads(2)

! CALCULATE FIELD PARAMS 
! Input : lambda, I0, Nc, CEP
E_max = sqrt(I0 *1.0D14 / 3.50945D16)  ! Maximum electric field amplitude in a.u.
omega = 2._dp*pi * 137.036_dp / (lambda * 1.0D-9 / 5.29177D-11)
!omega = 0.5_dp - 0.125_dp 
Up = E_max**2 / (4._dp * omega**2)  ! Ponderomotive energy in a.u.
rtUp = sqrt(Up)
laser_duration = 2._dp*pi * Nc/omega  ! Period of the laser field 
CEP = CEP * pi  ! Input is units of pi 

t_sim = t_sim * laser_duration  ! t_sim is input in units of laser_duration 
Nt = int(t_sim / dt) + 2  ! Nr. of steps during simulation 

max_krylov = N_krylov

write(*,*) 'dt : ', dt
write(*,*) 'Up :', Up
write(*,*) 'omega : ', omega
write(*,*) 'Laser duration :', laser_duration 
write(*,*) 'CEP :', CEP  
write(*,*) 'Number of splines : ', N_splines
write(*,*) 'Propagation residual tolerance : ', tol_residual


! CALCULATE MATRICES NEEDED FOR HAMILTONIAN APPLICATION 
allocate(S_cal(spline_order, N_splines_tot), H0_cal(spline_order, N_splines_tot), H1_cal(spline_order, N_splines_tot))
allocate(D0_cal(spline_order, N_splines_tot), D1_cal(spline_order, N_splines_tot))
allocate(x_quad(N_gauss_leg), x_arr(N_gauss_leg), w_quad(N_gauss_leg), w_arr(N_gauss_leg), V_arr(N_gauss_leg))
allocate(spline_arr(spline_order, N_gauss_leg), deriv_arr(spline_order, N_gauss_leg))

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
            ! Overlap 
            s_ij =  sum(spline_arr(m,:) * spline_arr(n,:) * w_arr)
            
            ! Atomic hamiltonian 
            h0_ij = 0.5_dp * sum(deriv_arr(n,:) * deriv_arr(m,:) * w_arr)  ! Kinetic term 
            h0_ij = h0_ij + sum(spline_arr(m,:)*spline_arr(n,:)*V_arr * w_arr)  ! Potential term 
            h1_ij = sum(spline_arr(m,:)*spline_arr(n,:)/x_arr**2 * w_arr)  ! Centrifugal term
            
            ! Laser interaction
            d0_ij = sum(spline_arr(m,:) * deriv_arr(n,:) * w_arr)  
            d1_ij = sum(spline_arr(m,:)*spline_arr(n,:)/x_arr * w_arr) 

            ! Now save in correct places 
            j = n + k-1  ! Index juggeling to save in band storage scheme, upper triangle, column in column, diagnoal in bottom
            i = spline_order + m - n 
            S_cal(i,j) = S_cal(i,j) + s_ij
            H0_cal(i,j) = H0_cal(i,j) + h0_ij 
            H1_cal(i,j) = H1_cal(i,j) + h1_ij
            D0_cal(i,j) = D0_cal(i,j) + d0_ij 
            D1_cal(i,j) = D1_cal(i,j) + d1_ij 
        end do  
    end do

    knot_init = knot_init + 1  ! Update when moving along on radial grid 
end do 

! Remove the first and last spline
allocate(S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines))
allocate(D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines))
S_mat = S_cal(:, 2:N_splines_tot-1)
H0_mat = H0_cal(:, 2:N_splines_tot-1)
H1_mat = H1_cal(:, 2:N_splines_tot-1)
D0_mat = D0_cal(:, 2:N_splines_tot-1)
D1_mat = D1_cal(:, 2:N_splines_tot-1)
do i=1, spline_order-1
    S_mat(spline_order-i, i) = 0._dp
    H0_mat(spline_order-i, i) = 0._dp 
    H1_mat(spline_order-i, i) = 0._dp
    D0_mat(spline_order-i, i) = 0._dp
    D1_mat(spline_order-i, i) = 0._dp
end do
deallocate(H0_cal, H1_cal, S_cal, D0_cal, D1_cal)

! PREPARE TIME PROPAGATION 
! Make list of angular coefficients  
allocate(alpha_list(l_max+1))
do i=0, l_max
    alpha_list(i+1) = (i+1._dp) / sqrt((2._dp*i+3._dp)*(2._dp*i+1._dp))
end do 

! Setup initial state 
allocate(psi(N_splines, l_max+1))
psi = 0._dp
psi(:,1) = psi0_real 

norm = calculate_norm_complex(psi, S_mat, spline_order, N_splines)
write(*,*) 'Initial norm : ', norm 

! Obtain the Cholesky factorization of the overlap matrix S 
allocate(S_cholesky(spline_order, N_splines))
S_cholesky = S_mat 
call Cholesky_factor_band_positive_def(S_cholesky, N_splines, spline_order)

!write(*,*) 
!do i=1, spline_order 
!    !write(*,"(15(F8.3))") S_mat(i,:)
!    write(*,"(15(F10.6))") D0_mat(i,:)
!end do

! Perform LU factorizations of the preconditioner (CN matrix without laser coupling)
allocate(LU_arr(3*(spline_order-1)+1, N_splines, l_max+1))
allocate(pivot_arr(N_splines, l_max+1))
allocate(H1_cal(spline_order, N_splines))
do i=1, l_max+1
    H1_cal = H0_mat + (i-1._dp)*i/2._dp * H1_mat  ! Hamiltonian for a give angular component 
    call LU_factor_Crank_Nicolson(H1_cal, S_mat, LU_arr(:,:,i), pivot_arr(:,i), dt, N_splines, spline_order)
end do
deallocate(H1_cal)

! TIME PROPAGATION using GMRES on CN matrix scheme 
write(*,*) 'Total nr. of steps : ', Nt
allocate(l_norms(l_max+1))

allocate(population(2, int(Nt/500)), time_list(int(Nt/500)))
population(1,1) = 1._dp 
population(2,1) = 0._dp 
time_list(1) = 0._dp 

! Time propagation using GMRES  
do i=1, Nt
    ! For each step we call GMRES until converged to a solution 
    if (mod(i,500) == 0) then 
        write(*,*)
        write(*,*) i, '/', Nt
        call calculate_norm_l_components(l_norms, psi, S_mat, spline_order, N_splines)
        write(*,*) real(l_norms(1)), real(l_norms(2)), real(l_norms(3))
        write(*,*) calculate_norm_complex(psi, S_mat, spline_order, N_splines)
        write(*,*) residual
        write(*,*) N_steps, N_krylov
        population(1, int(i/500)+1) = real(l_norms(1))
        population(2, int(i/500)+1) = real(l_norms(2))
        time_list(int(i/500)+1) = (i-1)*dt
    end if
    
    !do j=1, 100
        call GMRES_precond(psi, residual, (i-1._dp)*dt, dt, N_krylov, N_splines, l_max, S_mat, H0_mat, &
                           H1_mat, D0_mat, D1_mat, alpha_list, LU_arr, pivot_arr)
        !if (j > 1) then 
        !    write(*,*) j, residual 
        !end if 
        
        ! Adjust Kyrlov space 
        if (residual > res_bound_low .and. N_krylov /= 50) then 
            N_krylov = N_krylov + 1 
        else if (residual < res_bound_upper .and. N_krylov /= 8) then 
            N_krylov = N_krylov - 1 
        end if 
        
        if (residual > tol_residual) then 
            write(*,*) 'GMRES DID NOT CONVERGE!', N_krylov         
            stop 
        end if 
        
        !if (j == 100) then 
        !    write(*,*) 'GMRES NOT CONVERGED! STOPPING!'
        !    write(*,*) 'Current time step : ', i
        !    write(*,*) 'Current residual : ', residual 
        !    stop 
        !end if 
    !end do  
end do

! Time propagation using adaptive split-Lanczos 
!do i=1, Nt 
!    if (mod(i,100) == 0) then 
!        write(*,*)
!        write(*,*) i, '/', Nt
!        call calculate_norm_l_components(l_norms, psi, S_mat, spline_order, N_splines)
!        write(*,*) real(l_norms(1)), real(l_norms(2)), real(l_norms(3))
!        write(*,*) N_steps 
!    end if 
!
!    ! Split steps: First Crank-Nicolson then Lanczos and last CN
!    call Crank_Nicolson_step(psi, dt/2._dp, S_mat, H1_mat, pivot_arr, LU_arr, N_splines, spline_order, l_max)
!    
!    !if (i==10000) then 
!    
!        ! Test norm? 
!        write(*,*) calculate_norm_complex(psi, S_mat, spline_order, N_splines)
!        do j=1, 1000
!            call Crank_Nicolson_step(psi, dt/2._dp, S_mat, H1_mat, pivot_arr, LU_arr, N_splines, spline_order, l_max)
!        end do
!        write(*,*) calculate_norm_complex(psi, S_mat, spline_order, N_splines)
!        
!        stop
!    !end if  
!
!    ! The Lanczos step automiatically splits in smaller step until convergence 
!    time = dt*(i-1._dp)
!    dt_remain = dt 
!    N_steps = 1
!    do 
!        call lanczos_step_dynamic(psi, time, dt_remain, dt_taken, N_krylov, N_splines, l_max, S_mat, S_cholesky, H0_mat, H1_mat, &
!                        D0_mat, D1_mat, alpha_list)
!        time = time + dt_taken  ! Update with the actual step taken 
!        dt_remain = dt*i - time 
!        if (dt_remain < 1.e-12_dp) then 
!            exit
!        end if 
!        N_steps = N_steps + 1
!
!        if (N_steps > 50) then 
!            write(*,*) dt_remain 
!        end if 
!    end do 
!
!    ! Last half CN step 
!    call Crank_Nicolson_step(psi, dt/2._dp, S_mat, H1_mat, pivot_arr, LU_arr, N_splines, spline_order, l_max)
!
!end do

write(*,*)
write(*,*) 'Total norm : ', calculate_norm_complex(psi, S_mat, spline_order, N_splines)

write(*,*) 'Final norms'
call calculate_norm_l_components(l_norms, psi, S_mat, spline_order, N_splines)

do i=1, l_max+1
    write(*,*) real(l_norms(i))
end do

! SAVE FINAL STATE 
open(1, file='propagation_data/WF_propagated.dat', form='unformatted')
!open(2, file='propagation_data/WF_propagated_imag.dat', form='unformatted')
write(1) psi
!write(2) aimag(psi) 
close(1)
!close(2)

! Save population test 
open(2, file='propagation_data/pop.txt')
open(3, file='propagation_data/time_list.txt')
do i=1, size(population(1,:))
    write(2,*) population(1,i), population(2,i)
    write(3,*) time_list(i)
end do 
close(2)
close(3)


! TEST OVERLAP WITH INITIAL STATE? 
allocate(psi_init(N_splines), psi_cal(N_splines))
psi_init = psi0_real

call band_sym_matrix_vec_prod(S_mat, psi(:,1), psi_cal, N_splines, spline_order)
write(*,*) 'Overlap with init state : ', abs(sum(psi_init * psi_cal))**2


contains 

! LAPACK routine for calculating the Cholesky factorization of a banded positive definite Hermitian matrix 
subroutine Cholesky_factor_band_positive_def(mat, N, N_band)
    integer, intent(in) :: N, N_band
    complex(dp), intent(inout) :: mat(N_band, N)
    integer :: INFO 

    call ZPBTRF('U', N, N_band-1, mat, N_band, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Error in Cholesky factorization!'
        stop 
    end if 
end subroutine Cholesky_factor_band_positive_def

! Simple routine to apply S to each angular component of a state vector 
subroutine apply_S(psi, psi_res, S_mat, N_splines, spline_order, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    complex(dp), intent(in) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines)
    complex(dp), intent(out) :: psi_res(N_splines, l_max+1)
    integer :: li 
    
    do li=1, l_max+1 
        call band_sym_matrix_vec_prod(S_mat, psi(:,li), psi_res(:,li), N_splines, spline_order)
    end do
end subroutine apply_S

! Routine for applying the Hamiltonian to a state. Explicitly assumes linearly polarized laser!
subroutine apply_Hamiltonian_no_centri(psi, psi_H, t, H0_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    real(dp), intent(in) :: t
    complex(dp), intent(in) :: psi(N_splines, l_max+1)
    complex(dp), intent(out) :: psi_H(N_splines, l_max+1)
    real(dp), intent(in) :: H0_mat(spline_order, N_splines), alpha_list(l_max+1)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines)

    complex(dp) :: psi_d0(N_splines, l_max+1), psi_d1(N_splines, l_max+1), A_t 
    complex(dp) :: temp_vec(N_splines)
    integer :: i, l 

    ! Precalculate dipole matrix aciting on all angular components
    do i=1, l_max+1
        call band_antisym_matrix_vec_prod(D0_mat, psi(:,i), psi_d0(:,i), N_splines, spline_order)
        call band_sym_matrix_vec_prod(D1_mat, psi(:,i), psi_d1(:,i), N_splines, spline_order)
    end do 

    ! Now apply Hamiltonian to each angular component
    A_t = cmplx(0._dp, -vector_pot(t, omega, Nc, rtUp, CEP), dp)
    
    do l=0, l_max
        ! First atomic Hamiltonian 
        call band_sym_matrix_vec_prod(H0_mat, psi(:,l+1), psi_H(:,l+1), N_splines, spline_order)
        
        ! Then laser coupling 
        temp_vec = 0._dp 
        if (l /= 0) then  ! Coupling to lower l 
            temp_vec =  alpha_list(l) * (psi_d0(:,l) - (l*1._dp)*psi_d1(:,l))  
        end if 
        if (l /= l_max) then  ! Coupling to higher l 
            temp_vec = temp_vec + alpha_list(l+1) * (psi_d0(:,l+2) + (l+1._dp)*psi_d1(:,l+2)) 
        end if 
        
        ! Add to final result 
        psi_H(:, l+1) = psi_H(:, l+1) + A_t*temp_vec
    end do 
end subroutine apply_Hamiltonian_no_centri

! Routine for obtaining LU factorization of Crank-Nicolson matrix  
! H_mat and S_mat are given in symmetric band storage 
subroutine LU_factor_Crank_Nicolson(H_mat, S_mat, res, pivot_arr, dt, N_splines, N_band)
    integer, intent(in) :: N_splines, N_band 
    real(dp), intent(in) :: dt, H_mat(N_band, N_splines), S_mat(N_band, N_splines)
    integer, intent(out) :: pivot_arr(N_splines)
    complex(dp), intent(out) :: res(3*(N_band-1)+1, N_splines)  ! Size needed for ZGBTRF
    integer :: LDAB, i, INFO

    LDAB = 3*(N_band-1)+1 

    ! First, build the full matrix S + i*dt/2 * H
    res(N_band:2*N_band-1, :) = S_mat + cmplx(0._dp, dt/2._dp, dp) * H_mat 
    do i = 1, N_band-1  ! Lower triangle 
        do j = 1, N_splines - i 
            res(2*N_band-1 + i, j) = res(2*N_band-1 - i, j+i)
        end do 
    end do 

    ! Then obtain LU factorization through ZGBTRF
    call ZGBTRF(N_splines, N_splines, N_band-1, N_band-1, res, LDAB, pivot_arr, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Problem in LU factorization for Crank-Nicolson propagator!'
        stop
    end if
end subroutine LU_factor_Crank_Nicolson 

! Routine for calculating the matrix vector product, when the matrix is band symmetric
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

! Routine for taking a Crank-Nicolson step, using the LU factorizations of CN-matrix 
! This step only invovles the centrifugal operator, does not mixing of angular components
subroutine Crank_Nicolson_step(psi, dt, S_mat, H1_mat, pivot_arr, LU_arr_centri, N_splines, spline_order, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    real(dp), intent(in) :: dt, S_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    integer, intent(in) :: pivot_arr(N_splines, l_max+1)
    complex(dp), intent(in) :: LU_arr_centri(3*(spline_order-1)+1, N_splines, l_max+1)
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    complex(dp) :: psi_cal(N_splines)
    complex(dp) :: H_cal(spline_order, N_splines)
    integer :: i, LDAB, INFO

    LDAB = 3*(spline_order-1)+1

    ! Attempt this in parallel? 

    !!$OMP PARALLEL DO PRIVATE(i, psi_cal, INFO, H_cal)
    do i=1, l_max+1
        ! First perform matrix vector multiplications 
        H_cal = S_mat - cmplx(0._dp, dt/2._dp, dp) * (i-1._dp)*i/2._dp * H1_mat
        call band_sym_complex_matrix_vec_prod(H_cal, psi(:,i), psi_cal, N_splines, spline_order)

        ! Now solve linear system of eqns. using pre-calculated LU-decomposition
        call ZGBTRS('N', N_splines, spline_order-1, spline_order-1, 1, LU_arr_centri(:,:,i), & 
                    LDAB, pivot_arr(:,i), psi_cal, N_splines, INFO)
        if (INFO /= 0) then 
            write(*,*) 'Error when solving linear systems in Crank-Nicolson propagation!'
        end if 
        psi(:,i) = psi_cal
    end do
    !!$OMP END PARALLEL DO 
end subroutine Crank_Nicolson_step

! Calculates the norm of the eigenstate expanded in the B-spline basis.  
function calculate_norm_complex(vec, S_mat, spline_order, N_splines) result(norm)
    integer, intent(in) :: N_splines, spline_order 
    real(dp), intent(in) :: S_mat(spline_order, N_splines)
    complex(dp), intent(in) :: vec(N_splines, l_max+1)
    complex(dp) :: norm, cal_vec(N_splines)
    integer :: li
    norm  = 0._dp 
    do li = 1, l_max+1
        call band_sym_matrix_vec_prod(S_mat, vec(:,li), cal_vec, N_splines, spline_order)
        norm = norm + sum(conjg(vec(:,li)) * cal_vec)
    end do
    !norm = sqrt(norm)
end function calculate_norm_complex

subroutine calculate_norm_l_components(l_norms, vec, S_mat, spline_order, N_splines)
    integer, intent(in) :: spline_order, N_splines
    real(dp), intent(in) :: S_mat(spline_order, N_splines)
    complex(dp), intent(in) :: vec(N_splines, l_max+1)
    complex(dp), intent(out) :: l_norms(l_max+1)
    integer :: li 
    complex(dp) :: cal_vec(N_splines)
    do li = 1, l_max+1
        call band_sym_matrix_vec_prod(S_mat, vec(:,li), cal_vec, N_splines, spline_order)
        l_norms(li) = sum(conjg(vec(:,li)) * cal_vec)
    end do
end subroutine calculate_norm_l_components

! Routine for calculating eigenvalues of a tridiagnoal symmetric matrix 
subroutine tridiag_symmetric_eig(diag, off_diag, eigh_mat, N)
    integer, intent(in) :: N 
    real(dp), intent(inout) :: diag(N), off_diag(N-1)
    complex(dp), intent(out) :: eigh_mat(N,N)
    real(dp) :: WORK(2*N-2)
    integer :: INFO 
    call ZSTEQR('I', N, diag, off_diag, eigh_mat, N, WORK, INFO)
    if (INFO /= 0) then 
        print *, 'Error in ZSTEQR'
        stop 
    end if
end subroutine tridiag_symmetric_eig

subroutine lanczos_step_dynamic(psi, t, dt, dt_taken, N_krylov, N_splines, l_max, S_mat, S_cholesky, & 
                                H0_mat, H1_mat, D0_mat, D1_mat, alpha_list)
    integer, intent(in) :: N_splines, l_max, N_krylov 
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines), alpha_list(l_max+1)
    complex(dp), intent(in) :: S_cholesky(spline_order, N_splines)
    real(dp), intent(in) :: t, dt
    real(dp), intent(out) :: dt_taken

    complex(dp) :: lanczos_basis(N_splines, l_max+1, N_krylov)  ! Array to contain Lanczos basis
    complex(dp) :: psi_H(N_splines, l_max+1), eigh_mat(N_krylov, N_krylov)
    complex(dp) :: diag_mat(N_krylov, N_krylov), res_mat(N_krylov, N_krylov)
    complex(dp) :: psi_cal(N_splines, l_max+1)
    real(dp) :: alpha_arr(N_krylov), beta_arr(N_krylov-1)
    integer :: INFO, i
    real(dp) :: eps_step = 1.e-9_dp 

    ! TEST VARIABLES 
    complex(dp) :: l_norms(l_max+1)
    complex(dp) :: test_mat(N_krylov, N_krylov)
    real(dp) :: beta_arr_ex(N_krylov-1), alpha_arr_ex(N_krylov)

    lanczos_basis(:,:,1) = psi 
    
    ! Apply H and determine alpha 
    call apply_Hamiltonian_no_centri(psi, psi_H, t, H0_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max) 
    alpha_arr(1) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi

    ! Solve linear system 
    call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Everybody panic!'
    end if 

    ! Determine the new, unnormalized vector 
    psi_H = psi_H - alpha_arr(1)*lanczos_basis(:,:,1) 
    
    ! Build rest of the Krylov basis 
    do i=2, N_krylov
        ! Obtain the beta coeff. by normalization requirement 
        beta_arr(i-1) = real(sqrt(calculate_norm_complex(psi_H, S_mat, spline_order, N_splines)))
           
        ! Get Lanczos vector by dividing with beta? 
        lanczos_basis(:,:,i) = psi_H / beta_arr(i-1)
        
        ! Start making new Lanczos vector and find alpha 
        call apply_Hamiltonian_no_centri(lanczos_basis(:,:,i), psi_H, t, H0_mat, D0_mat, D1_mat, &
                               N_splines, spline_order, alpha_list, l_max) 
        alpha_arr(i) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi
    
        ! Solve linear system and build next unnormalized basis element  
        if (i /= N_krylov) then 
            call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
            if (INFO /= 0) then 
                write(*,*) 'Everybody panic!'
            end if 
            psi_H = psi_H - beta_arr(i-1)*lanczos_basis(:,:,i-1) - alpha_arr(i)*lanczos_basis(:,:,i)    
        end if 
    end do

    ! Now diagnoalize the Hamilton in the Lanczos basis representation
    eigh_mat = 0._dp 
    do i=1, N_krylov
        eigh_mat(i,i) = 1._dp 
    end do 
    call tridiag_symmetric_eig(alpha_arr, beta_arr, eigh_mat, N_krylov)

    ! Make diagnoal matrix with time step on diagonal. First try full time step
    diag_mat = 0._dp
    dt_taken = dt 
    do 
        do i=1, N_krylov
            diag_mat(i,i) = exp(-cmplx(0._dp, alpha_arr(i), dp)*dt_taken)
        end do

        ! Now do matrix multiplication 
        res_mat = matmul(eigh_mat, matmul(diag_mat, conjg(transpose(eigh_mat))))

        ! Check if step if within tolerance and if so exit the propagation
        if (abs(res_mat(1, N_krylov)) < eps_step) then 
            ! Build propagated psi 
            psi = 0._dp 
            do i=1, N_krylov  
                psi = psi + res_mat(1,i) * lanczos_basis(:,:,i)
            end do
            exit
        end if 

        ! If convergence is not reached, reduce time step and try again! 
        dt_taken = dt_taken / 2._dp  
    end do
end subroutine lanczos_step_dynamic

subroutine lanczos_step(psi, t, dt, N_krylov, N_splines, l_max, S_mat, S_cholesky, H0_mat, H1_mat, D0_mat, D1_mat, alpha_list)
    integer, intent(in) :: N_splines, l_max, N_krylov 
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines), alpha_list(l_max+1)
    complex(dp), intent(in) :: S_cholesky(spline_order, N_splines)
    real(dp), intent(in) :: t, dt

    complex(dp) :: lanczos_basis(N_splines, l_max+1, N_krylov)  ! Array to contain Lanczos basis
    complex(dp) :: psi_H(N_splines, l_max+1), eigh_mat(N_krylov, N_krylov), diag_mat(N_krylov, N_krylov)
    complex(dp) :: psi_cal(N_splines, l_max+1)
    real(dp) :: alpha_arr(N_krylov), beta_arr(N_krylov-1)
    integer :: INFO, i

    ! TEST VARIABLES 
    complex(dp) :: l_norms(l_max+1)
    complex(dp) :: test_mat(N_krylov, N_krylov)
    real(dp) :: beta_arr_ex(N_krylov-1), alpha_arr_ex(N_krylov)

    lanczos_basis(:,:,1) = psi 
    
    ! Apply H and determine alpha 
    call apply_Hamiltonian_no_centri(psi, psi_H, t, H0_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max) 
    alpha_arr(1) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi

    ! Solve linear system 
    call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Everybody panic!'
    end if 

    ! Determine the new, unnormalized vector 
    psi_H = psi_H - alpha_arr(1)*lanczos_basis(:,:,1) 
    
    ! Build rest of the Krylov basis 
    do i=2, N_krylov
        ! Obtain the beta coeff. by normalization requirement 
        beta_arr(i-1) = real(sqrt(calculate_norm_complex(psi_H, S_mat, spline_order, N_splines)))
           
        ! Get Lanczos vector by dividing with beta? 
        lanczos_basis(:,:,i) = psi_H / beta_arr(i-1)
        
        ! Start making new Lanczos vector and find alpha 
        call apply_Hamiltonian_no_centri(lanczos_basis(:,:,i), psi_H, t, H0_mat, D0_mat, D1_mat, &
                               N_splines, spline_order, alpha_list, l_max) 
        alpha_arr(i) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi
    
        ! Solve linear system and build next unnormalized basis element  
        if (i /= N_krylov) then 
            call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
            if (INFO /= 0) then 
                write(*,*) 'Everybody panic!'
            end if 
            psi_H = psi_H - beta_arr(i-1)*lanczos_basis(:,:,i-1) - alpha_arr(i)*lanczos_basis(:,:,i)    
        end if 
    end do

    ! Now diagnoalize the Hamilton in the Lanczos basis representation
    eigh_mat = 0._dp 
    do i=1, N_krylov
        eigh_mat(i,i) = 1._dp 
    end do 
    call tridiag_symmetric_eig(alpha_arr, beta_arr, eigh_mat, N_krylov)

    ! Let us check that it worked? 
    !test_mat = 0._dp
    !do i=1, N_krylov 
    !    test_mat(i,i) = alpha_arr_ex(i)
    !    if (i /= 1) then 
    !        write(*,*) i,i-1
    !        test_mat(i,i-1) = beta_arr_ex(i-1)
    !    end if 
    !    if (i /= N_krylov) then 
    !        test_mat(i,i+1) = beta_arr_ex(i)
    !    end if 
    !end do 
    !eigh_mat = transpose(eigh_mat)
    !test_mat = matmul(eigh_mat, matmul(test_mat, transpose(conjg(eigh_mat))))
    !test_mat = matmul(transpose(conjg(eigh_mat)), matmul(test_mat, eigh_mat))

    !write(*,*) 
    !do i=1, N_krylov  
    !    write(*,"(15(F8.3))") real(test_mat(i,:))
    !    !write(*,"(15(F10.6))") eigh_mat(i,:)
    !end do
    !write(*,*) alpha_arr

    ! Make diagnoal matrix with time step on diagonal 
    diag_mat = 0._dp
    do i=1, N_krylov
        diag_mat(i,i) = exp(-cmplx(0._dp, alpha_arr(i), dp)*dt)
    end do

    ! Now do matrix multiplication 
    diag_mat = matmul(eigh_mat, matmul(diag_mat, conjg(transpose(eigh_mat))))

    ! Should calculate estimator here to check convergence have been reached...
    ! However cannot easily increase value of k ... 
    ! See Castro et. al 2004 for details on error estimate 

    ! First row of this matrix is coeffs. of the propagated state
    !write(*,*)
    !write(*,*) abs(diag_mat(1,1))**2, abs(diag_mat(2,1))**2 
    !write(*,*) sum(abs(diag_mat(:,1)))
    
    !if (abs(diag_mat(1,N_krylov)) > 1.e-10) then
    !    write(*,*) 'Increase Krylov basis size?'
    !    write(*,*) diag_mat(1,N_krylov)
    !end if 

    ! Build propagated psi 
    psi = 0._dp 
    do i=1, N_krylov  
        psi = psi + diag_mat(1,i) * lanczos_basis(:,:,i)
    end do
end subroutine lanczos_step

! Rotuine for propagating a state a step dt using the short iterative Lanczos propagator 
subroutine lanczos_step_old(psi, t, dt, N_krylov, N_splines, l_max, S_mat, S_cholesky, H0_mat, H1_mat, D0_mat, D1_mat, alpha_list)
    integer, intent(in) :: N_splines, l_max, N_krylov 
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines), alpha_list(l_max+1)
    complex(dp), intent(in) :: S_cholesky(spline_order, N_splines)
    real(dp), intent(in) :: t, dt

    complex(dp) :: lanczos_basis(N_splines, l_max+1, N_krylov)  ! Array to contain Lanczos basis
    complex(dp) :: psi_H(N_splines, l_max+1), eigh_mat(N_krylov, N_krylov), diag_mat(N_krylov, N_krylov)
    complex(dp) :: psi_cal(N_splines, l_max+1), q_n(N_splines, l_max+1)
    real(dp) :: alpha_arr(N_krylov), beta_arr(N_krylov-1)
    integer :: INFO, i

    ! TEST VARIABLES 
    complex(dp) :: l_norms(l_max+1)

    ! AS A TEST MAKE SURE IT PROPAGATES THE EIGENSTATE CORRECTLY. SEEMS LIKE ALPHA/BETA BECOMES 
    ! IMAGINARY AFTER THE FIRST TIME STEP! 
    ! SHOULD PROJECT OUT THE WHOLE STATE, MAYBE NEED TO INTRODUCE CHECK?

    lanczos_basis(:,:,1) = psi 
    
    ! Apply H and determine alpha 
    call apply_Hamiltonian(psi, psi_H, t, H0_mat, H1_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max) 
    alpha_arr(1) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi

    !write(*,*)
    !write(*,*) 'alpha : ', sum(conjg(psi)*psi_H)
    !write(*,*) 'alpha : ', sum(conjg(psi_H)*psi)
    !stop 

    ! Do matrix-vector products to determine q_n 
    call apply_S(psi, psi_cal, S_mat, N_splines, spline_order, l_max)
    q_n = psi_H - alpha_arr(1)*psi_cal 
    
    ! Build rest of the Krylov basis 
    do i=2, N_krylov
        ! To obtain beta coeff. we must now solve linear system, save in psi_H 
        psi_H = q_n 

        call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
        if (INFO /= 0) then 
            write(*,*) 'Everybody panic!'
        end if 
        
        beta_arr(i-1) = real(sqrt(sum(conjg(q_n)*psi_H)))
        !write(*,*) sqrt(sum(conjg(q_n)*psi_H))

        !if (abs(aimag(sqrt(sum(conjg(q_n)*psi_H)))) > 1.E-9) then 
        !    write(*,*) int(t/dt)
        !    write(*,*) 'Beta : ', sqrt(sum(conjg(q_n)*psi_H))
        !end if 
            
        ! Get Lanczos vector by dividing with beta? 
        lanczos_basis(:,:,i) = psi_H / beta_arr(i-1)

        ! Check the norm of the new Laczos vector 
        !call calculate_norm_l_components(l_norms, lanczos_basis(:,:,i), S_mat, spline_order, N_splines)
        !write(*,*) 'Norms : ', i, sum(l_norms)

        ! Start making new Lanczos vector and find alpha 
        call apply_Hamiltonian(lanczos_basis(:,:,i), psi_H, t, H0_mat, H1_mat, D0_mat, D1_mat, &
                               N_splines, spline_order, alpha_list, l_max) 
        alpha_arr(i) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi
    
        ! Do matrix-vector products to determine q_n 
        if (i /= N_krylov) then 
            q_n = psi_H - beta_arr(i-1)*psi_cal
            call apply_S(lanczos_basis(:,:,i), psi_cal, S_mat, N_splines, spline_order, l_max)
            q_n = q_n - alpha_arr(i)*psi_cal 
        end if 
    end do

    ! Now diagnoalize the Hamilton in the Lanczos basis representation
    eigh_mat = 0._dp 
    do i=1, N_krylov
        eigh_mat(i,i) = 1._dp 
    end do 
    call tridiag_symmetric_eig(alpha_arr, beta_arr, eigh_mat, N_krylov)

    !write(*,*)
    !write(*,*) alpha_arr 

    ! Make diagnoal matrix with time step on diagonal 
    diag_mat = 0._dp
    do i=1, N_krylov
        diag_mat(i,i) = exp(-cmplx(0._dp, alpha_arr(i), dp)*dt)
    end do

    ! Now do matrix multiplication 
    diag_mat = matmul(eigh_mat, matmul(diag_mat, conjg(transpose(eigh_mat))))

    ! Should calculate estimator here to check convergence have been reached...
    ! However cannot easily increase value of k ... 
    ! See Castro et. al 2004 for details on error estimate 

    ! First row of this matrix is coeffs. of the propagated state
    !write(*,*)
    !write(*,*) abs(diag_mat(1,1))**2, abs(diag_mat(2,1))**2 
    !write(*,*) sum(abs(diag_mat(:,1)))
    
    !if (abs(diag_mat(1,N_krylov)) > 1.e-10) then
    !    write(*,*) 'Increase Krylov basis size?'
    !    write(*,*) diag_mat(1,N_krylov)
    !end if 

    ! Build propagated psi 
    psi = 0._dp 
    do i=1, N_krylov  
        psi = psi + diag_mat(1,i) * lanczos_basis(:,:,i)
    end do
end subroutine lanczos_step_old

! Routine for applying the Hamiltonian to a state. Explicitly assumes linearly polarized laser!
subroutine apply_Hamiltonian(psi, psi_H, t, H0_mat, H1_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    real(dp), intent(in) :: t
    complex(dp), intent(in) :: psi(N_splines, l_max+1)
    complex(dp), intent(out) :: psi_H(N_splines, l_max+1)
    real(dp), intent(in) :: H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines), alpha_list(l_max+1)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines)

    complex(dp) :: psi_d0(N_splines, l_max+1), psi_d1(N_splines, l_max+1), A_t 
    complex(dp) :: temp_vec(N_splines)
    real(dp) :: Hl(spline_order, N_splines)
    integer :: i, l 

    ! Precalculate dipole matrix aciting on all angular components
    do i=1, l_max+1
        call band_antisym_matrix_vec_prod(D0_mat, psi(:,i), psi_d0(:,i), N_splines, spline_order)
        call band_sym_matrix_vec_prod(D1_mat, psi(:,i), psi_d1(:,i), N_splines, spline_order)
    end do 

    ! Now apply Hamiltonian to each angular component
    A_t = cmplx(0._dp, -vector_pot(t, omega, Nc, rtUp, CEP), dp)
    
    do l=0, l_max
        ! First atomic Hamiltonian 
        Hl = H0_mat + l*(l+1._dp)/2._dp * H1_mat
        call band_sym_matrix_vec_prod(Hl, psi(:,l+1), psi_H(:,l+1), N_splines, spline_order)
        
        ! Then laser coupling 
        temp_vec = 0._dp 
        if (l /= 0) then  ! Coupling to lower l 
            temp_vec =  alpha_list(l) * (psi_d0(:,l) - (l*1._dp)*psi_d1(:,l))  
        end if 
        if (l /= l_max) then  ! Coupling to higher l 
            temp_vec = temp_vec + alpha_list(l+1) * (psi_d0(:,l+2) + (l+1._dp)*psi_d1(:,l+2)) 
        end if 
        
        ! Add to final result 
        psi_H(:, l+1) = psi_H(:, l+1) + A_t*temp_vec
    end do 
end subroutine apply_Hamiltonian 

! Routine for applying the Hamiltonian to a state. Explicitly assumes linearly polarized laser!
subroutine apply_Cayley_Hamilton_form(psi, psi_H, t, dt, SH0_mat, H1_mat, D0_mat, D1_mat, & 
                                      N_splines, spline_order, alpha_list, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    real(dp), intent(in) :: t, dt
    complex(dp), intent(in) :: psi(N_splines, l_max+1)
    complex(dp), intent(out) :: psi_H(N_splines, l_max+1)
    real(dp), intent(in) :: H1_mat(spline_order, N_splines), alpha_list(l_max+1)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines)
    complex(dp), intent(in) :: SH0_mat(spline_order, N_splines)  ! S + i*dt/2 * H0  

    complex(dp) :: psi_d0(N_splines, l_max+1), psi_d1(N_splines, l_max+1), A_t 
    complex(dp) :: temp_vec(N_splines)
    complex(dp) :: Hl(spline_order, N_splines)
    integer :: i, l 

    ! Precalculate dipole matrix aciting on all angular components
    do i=1, l_max+1
        call band_antisym_matrix_vec_prod(D0_mat, psi(:,i), psi_d0(:,i), N_splines, spline_order)
        call band_sym_matrix_vec_prod(D1_mat, psi(:,i), psi_d1(:,i), N_splines, spline_order)
    end do 

    ! Now apply Hamiltonian to each angular component
    A_t = cmplx(0._dp, -vector_pot(t, omega, Nc, rtUp, CEP), dp)
    
    do l=0, l_max
        ! First atomic Hamiltonian 
        Hl = SH0_mat + cmplx(0._dp, dt/2._dp, dp)*l*(l+1._dp)/2._dp * H1_mat
        call band_sym_complex_matrix_vec_prod(Hl, psi(:,l+1), psi_H(:,l+1), N_splines, spline_order)
        
        ! Then laser coupling 
        temp_vec = 0._dp 
        if (l /= 0) then  ! Coupling to lower l 
            temp_vec =  alpha_list(l) * (psi_d0(:,l) - (l*1._dp)*psi_d1(:,l))  
        end if 
        if (l /= l_max) then  ! Coupling to higher l 
            temp_vec = temp_vec + alpha_list(l+1) * (psi_d0(:,l+2) + (l+1._dp)*psi_d1(:,l+2)) 
        end if 
        
        ! Add to final result 
        psi_H(:, l+1) = psi_H(:, l+1) + cmplx(0._dp, dt/2._dp, dp)*A_t*temp_vec
    end do 
end subroutine apply_Cayley_Hamilton_form 

! Subroutine for solving implicit Cayley-Hamilton eqns. of propagation using GMRES without preconditioner 
subroutine GMRES(psi, residual, t, dt, N_krylov, N_splines, l_max, S_mat, H0_mat, H1_mat, D0_mat, D1_mat, alpha_list)
    integer, intent(in) :: N_splines, l_max, N_krylov 
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines), alpha_list(l_max+1)
    real(dp), intent(in) :: t, dt
    real(dp), intent(out) :: residual 

    complex(dp) :: SH0_mat(spline_order, N_splines)
    complex(dp) :: krylov_basis(N_splines, l_max+1, N_krylov)
    complex(dp) :: psi_cal(N_splines, l_max+1), psi_cal1(N_splines, l_max+1)
    complex(dp) :: H_bar_krylov(N_krylov+1, N_krylov), h_ij, g_vec(N_krylov+1)
    real(dp) :: beta 
    integer :: i, j, INFO

    SH0_mat = S_mat + cmplx(0._dp, dt/2._dp, dp) * H0_mat

    ! First determine the RHS in Ax=b, i.e. (S - i*dt/2 * H) * psi(t)
    call apply_S(psi, psi_cal1, S_mat, N_splines, spline_order, l_max)
    call apply_Hamiltonian(psi, psi_cal, t, H0_mat, H1_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
    psi_cal1 = psi_cal1 - cmplx(0._dp, dt/2._dp, dp) * psi_cal 

    ! Determine initial residual and normalize
    call apply_Cayley_Hamilton_form(psi, psi_cal, t, dt, SH0_mat, H1_mat, D0_mat, D1_mat, N_splines, &
                                    spline_order, alpha_list, l_max)
    psi_cal = psi_cal1 - psi_cal
    beta = real(sqrt(sum(conjg(psi_cal)*psi_cal)))  ! Using L2 norm, not S normalized? 
    krylov_basis(:,:,1) = psi_cal / beta 

    ! Initialize H_krylov to 0? 
    ! Now perform Arnoldi using modified Gram-Schmidt 
    H_bar_krylov = 0._dp 
    do j=1, N_krylov
        ! First operator with CH matrix 
        call apply_Cayley_Hamilton_form(krylov_basis(:,:,j), psi_cal, t, dt, SH0_mat, H1_mat, & 
                                        D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
        do i=1, j
            h_ij = sum(conjg(psi_cal) * krylov_basis(:,:,i))
            psi_cal = psi_cal - h_ij * krylov_basis(:,:,i)
            H_bar_krylov(i,j) = h_ij
        end do

        ! Normalize and save new basis vector 
        h_ij = sqrt(sum(conjg(psi_cal) * psi_cal))  ! Should check if this is 0
        H_bar_krylov(j+1, j) = h_ij
        if (j /= N_krylov) then
            krylov_basis(:,:,j+1) = psi_cal / h_ij 
        end if  
    end do

    ! Test we have Hessenberg form 
    !write(*,*) 
    !do i=1, N_krylov+1
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") real(H_bar_krylov(i,:))
    !end do
    !write(*,*) 
    !do i=1, N_krylov 
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") aimag(H_bar_krylov(i,:))
    !end do
    !stop 

    ! Solve least square problem. First perform Givens rotations to turn Hessenberg (N_krylov+1, N_krylov)
    ! into (N_krylov, N_krylov) triangular system   
    g_vec = 0._dp 
    g_vec(1) = beta 
    call apply_Givens_rotations(H_bar_krylov, g_vec, N_krylov)


    ! To test Givens rotation
    !write(*,*) 
    !do i=1, N_krylov+1
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") real(H_bar_krylov(i,:))
    !end do
    !write(*,*) 
    !do i=1, N_krylov 
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") aimag(H_bar_krylov(i,:))
    !end do

    ! Solve triangular system using ZTRTRS 
    call ZTRTRS('U', 'N', 'N', N_krylov, 1, H_bar_krylov(1:N_krylov, :), N_krylov, g_vec(1:N_krylov), N_krylov, INFO)
    if (INFO /= 0) then 
        write(*,*) i, 'Problem in solving triangular least squares in GMRES!'
    end if 

    ! Now contstruct full solution, note last entry in g_vec is the risidual, not used in linear combi
    residual = abs(g_vec(N_krylov+1))
    do i=1, N_krylov
        psi = psi + g_vec(i)*krylov_basis(:,:,i)
    end do
end subroutine GMRES

! Solve linear system for x = M^(-1) * psi, with M beeing the preconditioner  
subroutine solve_precond_system(psi, res, LU_arr, pivot_arr, N_splines, spline_order, l_max)
    integer, intent(in) :: N_splines, spline_order, l_max
    complex(dp), intent(in) :: psi(N_splines, l_max+1), LU_arr(3*(spline_order-1)+1, N_splines, l_max+1)
    integer, intent(in) :: pivot_arr(N_splines, l_max+1)
    complex(dp), intent(out) :: res(N_splines, l_max+1)
    
    complex(dp) :: psi_cal(N_splines)
    integer :: LDAB, INFO, i

    LDAB = 3*(spline_order-1)+1

    !$OMP PARALLEL DO PRIVATE(psi_cal, INFO)
    do i=1, l_max+1
        psi_cal = psi(:,i)
        ! Solve linear system of eqns. using pre-calculated LU-decomposition of preconditioner
        call ZGBTRS('N', N_splines, spline_order-1, spline_order-1, 1, LU_arr(:,:,i), & 
                    LDAB, pivot_arr(:,i), psi_cal, N_splines, INFO)  ! MAKES ARRAY TEMP! 
        if (INFO /= 0) then 
            write(*,*) 'Error when solving linear system for preconditioner!'
        end if 
        res(:,i) = psi_cal
    end do
    !$OMP END PARALLEL DO
end subroutine solve_precond_system

! Subroutine for solving implicit Cayley-Hamilton eqns. of propagation using GMRES without preconditioner 
subroutine GMRES_precond(psi, residual, t, dt, N_krylov, N_splines, l_max, S_mat, H0_mat, H1_mat, D0_mat, D1_mat, alpha_list, &
                        LU_arr, pivot_arr)
    integer, intent(in) :: N_splines, l_max, N_krylov 
    complex(dp), intent(inout) :: psi(N_splines, l_max+1)
    real(dp), intent(in) :: S_mat(spline_order, N_splines), H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines)
    real(dp), intent(in) :: D0_mat(spline_order, N_splines), D1_mat(spline_order, N_splines), alpha_list(l_max+1)
    real(dp), intent(in) :: t, dt
    complex(dp), intent(in) :: LU_arr(3*(spline_order-1)+1, N_splines, l_max+1)
    integer, intent(in) :: pivot_arr(N_splines, l_max+1)
    real(dp), intent(out) :: residual 

    complex(dp) :: SH0_mat(spline_order, N_splines)
    complex(dp) :: krylov_basis(N_splines, l_max+1, N_krylov)
    complex(dp) :: psi_cal(N_splines, l_max+1), psi_cal1(N_splines, l_max+1)
    complex(dp) :: H_bar_krylov(N_krylov+1, N_krylov), h_ij, g_vec(N_krylov+1)
    real(dp) :: beta 
    integer :: i, j, INFO

    SH0_mat = S_mat + cmplx(0._dp, dt/2._dp, dp) * H0_mat

    ! First determine the RHS in Ax=b, i.e. (S - i*dt/2 * H) * psi(t)
    call apply_S(psi, psi_cal1, S_mat, N_splines, spline_order, l_max)
    call apply_Hamiltonian(psi, psi_cal, t, H0_mat, H1_mat, D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
    psi_cal1 = psi_cal1 - cmplx(0._dp, dt/2._dp, dp) * psi_cal  ! This is b in Ax=b

    ! Determine initial preconditioned residual and normalize
    call apply_Cayley_Hamilton_form(psi, psi_cal, t, dt, SH0_mat, H1_mat, D0_mat, D1_mat, N_splines, &
                                    spline_order, alpha_list, l_max)
    psi_cal = psi_cal1 - psi_cal
    call solve_precond_system(psi_cal, psi_cal1, LU_arr, pivot_arr, N_splines, spline_order, l_max)  ! Apply inverse preconditioner (M^(-1) * res)
    beta = real(sqrt(sum(conjg(psi_cal1)*psi_cal1)))  ! Using L2 norm, not S normalized? 
    krylov_basis(:,:,1) = psi_cal1 / beta 

    ! Initialize H_krylov to 0? 
    ! Now perform Arnoldi using modified Gram-Schmidt 
    H_bar_krylov = 0._dp 
    do j=1, N_krylov
        ! First operator with CH matrix 
        call apply_Cayley_Hamilton_form(krylov_basis(:,:,j), psi_cal1, t, dt, SH0_mat, H1_mat, & 
                                        D0_mat, D1_mat, N_splines, spline_order, alpha_list, l_max)
        
        ! As a temporary test apply inverse preconditioner here, do smarter if it works! 
        call solve_precond_system(psi_cal1, psi_cal, LU_arr, pivot_arr, N_splines, spline_order, l_max)
        
        do i=1, j
            h_ij = sum(conjg(psi_cal) * krylov_basis(:,:,i))
            psi_cal = psi_cal - h_ij * krylov_basis(:,:,i)
            H_bar_krylov(i,j) = h_ij
        end do

        ! Normalize and save new basis vector 
        h_ij = sqrt(sum(conjg(psi_cal) * psi_cal))  ! Should check if this is 0
        H_bar_krylov(j+1, j) = h_ij
        if (j /= N_krylov) then
            krylov_basis(:,:,j+1) = psi_cal / h_ij 
        end if  
    end do

    ! Test we have Hessenberg form 
    !write(*,*) 
    !do i=1, N_krylov+1
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") real(H_bar_krylov(i,:))
    !end do
    !write(*,*) 
    !do i=1, N_krylov 
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") aimag(H_bar_krylov(i,:))
    !end do
    !stop 

    ! Solve least square problem. First perform Givens rotations to turn Hessenberg (N_krylov+1, N_krylov)
    ! into (N_krylov, N_krylov) triangular system   
    g_vec = 0._dp 
    g_vec(1) = beta 
    call apply_Givens_rotations(H_bar_krylov, g_vec, N_krylov)


    ! To test Givens rotation
    !write(*,*) 
    !do i=1, N_krylov+1
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") real(H_bar_krylov(i,:))
    !end do
    !write(*,*) 
    !do i=1, N_krylov 
    !    !write(*,"(15(F8.3))") S_mat(i,:)
    !    write(*,"(15(F10.6))") aimag(H_bar_krylov(i,:))
    !end do

    ! Solve triangular system using ZTRTRS 
    call ZTRTRS('U', 'N', 'N', N_krylov, 1, H_bar_krylov(1:N_krylov, :), N_krylov, g_vec(1:N_krylov), N_krylov, INFO)
    if (INFO /= 0) then 
        write(*,*) i, 'Problem in solving triangular least squares in GMRES!'
    end if 

    ! Now contstruct full solution, note last entry in g_vec is the preconditioned risidual, not used in linear combi
    residual = abs(g_vec(N_krylov+1))
    do i=1, N_krylov
        psi = psi + g_vec(i)*krylov_basis(:,:,i)
    end do
end subroutine GMRES_precond 

subroutine apply_Givens_rotations(mat_hess, vec, N_krylov)
    integer, intent(in) :: N_krylov 
    complex(dp), intent(inout) :: mat_hess(N_krylov+1, N_krylov), vec(N_krylov+1)
    
    complex(dp) :: si, ci, sqrt_fac, temp
    integer :: i, j
    
    ! Loop over the diagonal elements, around which we apply the rotations (2x2 matrix products)
    do i=1, N_krylov
        ! Calculate cos/sin factors for rotation matrix 
        sqrt_fac = sqrt(abs(mat_hess(i,i))**2 + mat_hess(i+1,i)**2)
        si = mat_hess(i+1, i) / sqrt_fac 
        ci = mat_hess(i, i) / sqrt_fac 

        ! Apply the rotation, operating column wise 
        do j=i, N_krylov 
            temp = conjg(ci)*mat_hess(i,j) + si*mat_hess(i+1,j)
            mat_hess(i+1,j) = -si*mat_hess(i,j) + ci*mat_hess(i+1,j)
            mat_hess(i,j) = temp 
        end do

        ! Also apply to the RHS vector 
        temp = conjg(ci)*vec(i) + si*vec(i+1)
        vec(i+1) = -si*vec(i) + ci*vec(i+1)
        vec(i) = temp 
    end do

end subroutine apply_Givens_rotations


end program propagate 