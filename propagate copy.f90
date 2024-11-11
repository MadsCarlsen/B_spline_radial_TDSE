program propagate 
use utils 
use B_spline 
use GaussLegendre
use potential 
use linalg_routines
implicit none 

real(dp) :: r_max 
integer :: N_splines_tot, N_splines, N_points, spline_order, N_knots, N_gauss_leg, knot_init
real(dp), allocatable :: knot_list(:), r_list(:)

real(dp), allocatable :: x_quad(:), w_quad(:), x_arr(:), w_arr(:), V_arr(:), psi0_real(:)
real(dp), allocatable :: spline_arr(:,:), deriv_arr(:,:)
complex(dp), allocatable :: psi0(:)

real(dp), allocatable :: H0_cal(:,:), H1_cal(:,:), D0_cal(:,:), D1_cal(:,:), S_cal(:,:)  ! Temp matrices 
real(dp) :: h0_ij, h1_ij, d0_ij, d1_ij, s_ij
real(dp), allocatable :: S_mat(:,:), H0_mat(:,:), H1_mat(:,:), D0_mat(:,:), D1_mat(:,:)  ! Mats without first/last spline
complex(dp), allocatable :: S_cholesky(:,:)  ! Complex since we use it to solve complex Sx=b with LAPACK
integer :: i, j, k, n, m

! Varibales for propagation 
integer :: l_max, Nt, N_krylov 
real(dp) :: dt, t_sim
real(dp), allocatable :: alpha_list(:)
complex(dp), allocatable :: psi(:,:)

! Field params 
real(dp) :: lambda, I0, CEP, omega, E_max, Up, rtUp, laser_duration 
integer :: Nc

! TEST 
complex(dp), allocatable :: l_norms(:)
complex(dp), allocatable :: test(:), test_res(:), test_arr(:,:), test_res_arr(:,:)
complex(dp), allocatable :: norm, func_eval_arr(:)

! LOAD SIMULATION DATA 

! Field params 
namelist /PROPAGATION/ lambda, I0, Nc, CEP, l_max, dt, t_sim, N_krylov 
open(file='settings.nml', unit=1)
read(nml=PROPAGATION, unit=1)
close(1)

! Grid and B spline params from eigenstate data dump
open(1, file='data/sim_info.txt')
read(1,*) N_splines_tot, N_points, spline_order, N_knots, N_gauss_leg, r_max
close(1)
N_splines = N_splines_tot - 2

allocate(r_list(N_points), knot_list(N_knots), psi0_real(N_splines))
open(1, file='data/r_list.dat', form='unformatted')
open(2, file='data/knot_list.dat', form='unformatted')
open(3, file='data/init_state.dat', form='unformatted')
read(1) r_list 
read(2) knot_list 
read(3) psi0_real 
close(1)
close(2)
close(3)


! CALCULATE FIELD PARAMS 
! Input : lambda, I0, Nc, CEP
E_max = sqrt(I0 *1.0D14 / 3.50945D16)  ! Maximum electric field amplitude in a.u.
omega = 2._dp*pi * 137.036_dp / (lambda * 1.0D-9 / 5.29177D-11)
Up = E_max**2 / (4._dp * omega**2)  ! Ponderomotive energy in a.u.
rtUp = sqrt(Up)
laser_duration = 2._dp*pi * Nc/omega  ! Period of the laser field 

t_sim = t_sim * laser_duration  ! t_sim is input in units of laser_duration 
Nt = int(t_sim / dt) + 1  ! Nr. of steps during simulation 


! Test 
dt = 2._dp*pi/omega / 1000._dp
Nt = int(t_sim / dt) + 1  ! Nr. of steps during simulation 

write(*,*) 'dt : ', dt
write(*,*) 'Up :', Up
write(*,*) 'omega : ', omega
write(*,*) 'Laser duration :', laser_duration 

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

!
!write(*,*) 
!do i=1, spline_order 
!    write(*,"(15(F8.3))") real(S_cholesky(i,:))
!end do

! TEST 
!allocate(test(N_splines), test_res(N_splines))
!test = [1._dp, 2._dp, 3._dp, 4._dp]
!call band_antisym_matrix_vec_prod(D0_mat, test, test_res, N_splines, spline_order)
!write(*,*) test_res
!stop 

! TIME PROPAGATION by short iterative Lanczos progator 

! Test initial normalization 

! do loop over Lanczos steps.
! Save at end of laser pulse, project in other program? 
allocate(l_norms(l_max+1))

write(*,*) Nt
do i=1, Nt 
    !write(*,*) i
    if (mod(i, 200) == 0 ) then 
        write(*,*)
        write(*,*) i, '/', Nt
        call calculate_norm_l_components(l_norms, psi, S_mat, spline_order, N_splines)
        
        do j=1, l_max+1
            write(*,*) j, l_norms(j)
        end do
        write(*,*)
        write(*,*) sum(l_norms)
        write(*,*) calculate_norm_complex(psi, S_mat, spline_order, N_splines)
        !stop
        exit 
    end if 

    call lanczos_step(psi, dt*(i-0.5_dp), dt, N_krylov, N_splines, l_max, S_mat, S_cholesky, H0_mat, H1_mat, &
                    D0_mat, D1_mat, alpha_list)


    !if (i == 3) then
    !    stop 
    !end if  
    !    call calculate_norm_l_components(l_norms, psi, S_mat, spline_order, N_splines)
    !    write(*,*) '================================'
    !    write(*,*) l_norms(1), l_norms(2)
    ! 
    !    stop 
    !end if 
end do


! Evaluate eigenstate on radial grid and save 
deallocate(r_list)
N_points = 1000
allocate(r_list(N_points))
call linspace(r_list, 0._dp, 90._dp, N_points)
allocate(func_eval_arr(size(r_list)))
call eval_complex_function_B_basis(func_eval_arr, psi(:,3), r_list, knot_list, spline_order, N_splines)

open(file='r_list_prop.dat', unit=1)
open(file='func_list_prop_real.dat', unit=2)
open(file='func_list_prop_imag.dat', unit=3)
do i=1, N_points 
    write(1,*) r_list(i)
    write(2,*) real(func_eval_arr(i))
    write(3,*) imag(func_eval_arr(i))
end do



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

    ! Solve linear system 
    call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Everybody panic!'
    end if 

    ! TEST SOLUTION OF LINEAR SYSTEM OF EQUATIONS! 
    
    ! Determine the new, unnormalized vector 
    psi_H = psi_H - alpha_arr(1)*lanczos_basis(:,:,1) 
    
    !write(*,*) calculate_norm_complex(psi_H, S_mat, spline_order, N_splines)
    
    ! Build rest of the Krylov basis 
    do i=2, N_krylov
        ! Obtain the beta coeff. by normalization requirement 
        beta_arr(i-1) = real(sqrt(calculate_norm_complex(psi_H, S_mat, spline_order, N_splines)))
           
        ! Get Lanczos vector by dividing with beta? 
        lanczos_basis(:,:,i) = psi_H / beta_arr(i-1)
        
        ! Start making new Lanczos vector and find alpha 
        call apply_Hamiltonian(lanczos_basis(:,:,i), psi_H, t, H0_mat, H1_mat, D0_mat, D1_mat, &
                               N_splines, spline_order, alpha_list, l_max) 
        alpha_arr(i) = real(sum(conjg(psi)*psi_H))  ! Calculate alpha while we have H*psi
    
        !write(*,*) sum(conjg(psi)*psi_H)

        ! Solve linear system and build next unnormalized basis element  
        if (i /= N_krylov) then 
            call ZPBTRS('U', N_splines, spline_order-1, l_max+1, S_cholesky, spline_order, psi_H, N_splines, INFO)
            if (INFO /= 0) then 
                write(*,*) 'Everybody panic!'
            end if 
            psi_H = psi_H - beta_arr(i-1)*lanczos_basis(:,:,i-1) - alpha_arr(i)*lanczos_basis(:,:,i)    
        end if 
    end do

    !beta_arr_ex = beta_arr
    !alpha_arr_ex = alpha_arr

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

end program propagate 