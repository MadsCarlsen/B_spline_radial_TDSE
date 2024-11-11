program eigenstate 
use utils 
use B_spline
use GaussLegendre 
use potential 
use linalg_routines
implicit none 

real(dp) :: r_max, E_max, eps 
integer :: N_splines, N_points, spline_order, l_max, N_gauss_leg, max_iter, para_index
real(dp), allocatable :: knot_list(:), r_list(:) 
real(dp), allocatable :: H0_mat(:,:), H1_mat(:,:), S_mat(:,:), H_mat(:,:)

logical :: determine_bound, determine_continuum

real(dp), allocatable :: eigh_vals(:), eigh_vec(:), eigh_bound(:,:), eigh_cont(:,:)
integer :: i, j, N_bound, N_cont, li
integer, allocatable ::  N_bound_arr(:), N_cont_arr(:)
real(dp), allocatable :: eigh_val_bound(:), eigh_val_cont(:) 
integer, allocatable :: cont_index(:)
real(dp) :: DOS 
character(len=100) :: save_file

! First setup radial grid and knot sequence 
namelist /GRID/ r_max, N_points, spline_order, N_gauss_leg, l_max, para_index
namelist /EIGENSTATES/ determine_bound, determine_continuum, E_max, eps, max_iter
open(file='settings.nml', unit=1)
read(nml=GRID, unit=1)
read(nml=EIGENSTATES, unit=1)
close(1)

N_gauss_leg = spline_order + N_gauss_leg 

allocate(r_list(N_points))
!call knot_squence_lin(knot_list, N_splines, r_list, N_points, spline_order, r_max)
call knot_squence_lin_para(knot_list, N_splines, r_list, N_points, spline_order, r_max, para_index)
N_splines = N_splines - 2  ! Set first and last B-spline to 0 to get correct boundary conditions 

! DETERMINE MATRIX ELEMENTS 
allocate(H0_mat(spline_order, N_splines), H1_mat(spline_order, N_splines), S_mat(spline_order, N_splines))
call build_matrices(H0_mat, H1_mat, S_mat, knot_list, r_list, N_points, N_splines, spline_order, N_gauss_leg)

!write(*,*) 
!do i=1, spline_order 
!    write(*,"(15(F8.3))") H_mat(i,:)
!end do

! GROUND STATE CALCULATIONS (Always preform this)
write(*,*) 'Determining ground state...'
allocate(H_mat(spline_order, N_splines))
H_mat = H0_mat

! Determine eigenvalues
allocate(eigh_vals(N_splines))
call generalized_symband_eigh(eigh_vals, N_splines, spline_order, H_mat, S_mat)

! Determine ground state through inverse iteration  
allocate(eigh_vec(N_splines))
call inverse_iteration(eigh_vec, eigh_vals(1), spline_order, N_splines, H_mat, S_mat, eps, max_iter)

write(*,*) 'Ground state energy : ', eigh_vals(1)

! Save (type_l_Nr)
open(file='eigenstate_data/ground_state.dat', unit=1, form='unformatted')
write(1) eigh_vec 
close(1)


! CALCULATIONS FOR BOUND STATE AND/OR CONTINIUUM STATES 
if (determine_bound .or. determine_continuum) then 
    write(*,*) 'Determining excited states...'
    allocate(N_bound_arr(l_max+1), N_cont_arr(l_max+1))

    ! Loop over all angular momentum values 
    do li=1, l_max+1        
        ! Check if we should diagonalize
        if (li /= 1) then 
            H_mat = H0_mat + (li-1._dp)*li/2._dp * H1_mat 
            call generalized_symband_eigh(eigh_vals, N_splines, spline_order, H_mat, S_mat)
        end if 

        ! Sort the eigenvalues into bound/continuum states 
        call sort_eigh_vals(eigh_vals, E_max, N_splines, N_bound, N_cont, eigh_val_bound, eigh_val_cont, cont_index)
        N_bound_arr(li) = N_bound 
        N_cont_arr(li) = N_cont 

        ! Determine bound eigenstates 
        if (determine_bound .and. (N_bound /= 0)) then 
            allocate(eigh_bound(N_bound, N_splines))
            do i=1, N_bound
                call inverse_iteration(eigh_vec, eigh_val_bound(i), spline_order, N_splines, H_mat, S_mat, eps, max_iter)
                eigh_bound(i,:) = eigh_vec
            end do

            ! Save states 
            write(save_file, "(A,I0,A)") 'eigenstate_data/state_bound_', li-1, '.dat'
            open(file=save_file, unit=1, form='unformatted')
            write(1) eigh_bound
            close(1)

            ! Save energies
            write(save_file, "(A,I0,A)") 'eigenstate_data/energies_bound_', li-1, '.txt'
            open(file=save_file, unit=2)
            do i=1, N_bound 
                write(2,*) eigh_val_bound(i)
            end do 
            close(2)

            deallocate(eigh_bound)
        end if 

        ! Determine continuum states 
        if (determine_continuum .and. (N_cont /= 0) ) then 
            allocate(eigh_cont(N_cont, N_splines))
            do i=1, N_cont 
                call inverse_iteration(eigh_vec, eigh_val_cont(i), spline_order, N_splines, H_mat, S_mat, eps, max_iter)
                eigh_cont(i,:) = eigh_vec

                ! Normalize to energy-scale by multiplying with sqrt of density-of-states (2nd order approx for now)
                j = cont_index(i)
                if (j == 1) then  ! FIX THIS!? 
                    ! If at endpoint use forward difference  
                    DOS = 1._dp / (eigh_vals(j+1) - eigh_vals(j))
                else 
                    ! Else use central difference 
                    DOS = 2._dp / (eigh_vals(j+1) - eigh_vals(j-1)) 
                end if 
                eigh_cont(i,:) = abs(sqrt(DOS)) * eigh_cont(i,:)
            end do

            ! Save states 
            write(save_file, "(A,I0,A)") 'eigenstate_data/state_continuum_', li-1, '.dat'
            open(file=save_file, unit=1, form='unformatted')
            write(1) eigh_cont
            close(1)

            ! Save energies
            write(save_file, "(A,I0,A)") 'eigenstate_data/energies_continuum_', li-1, '.txt'
            open(file=save_file, unit=2)
            do i=1, N_cont
                write(2,*) eigh_val_cont(i)
            end do 
            close(2)

            deallocate(eigh_cont, cont_index)
        end if 

        ! Deallocate for next iteration 
        deallocate(eigh_val_bound)
        deallocate(eigh_val_cont)
    end do

    ! Last save files keeping track of nr of bound/continuum states for different l values
    open(file='eigenstate_data/nr_bound_states.txt', unit=1)
    open(file='eigenstate_data/nr_continuum_states.txt', unit=2)
    
    if (.not. determine_bound) then 
        N_bound_arr = 0
    end if 
    if (.not. determine_continuum) then 
        N_cont_arr = 0 
    end if 

    do li=1, l_max+1
        write(1,*) N_bound_arr(li)
        write(2,*) N_cont_arr(li)
    end do

end if 


! Evaluate eigenstate on radial grid and save 
!deallocate(r_list)
!allocate(r_list(N_points))
!call linspace(r_list, 0._dp, 10._dp, N_points)
!allocate(func_eval_arr(size(r_list)))
!call eval_function_B_basis(func_eval_arr, eigh_vec, r_list, knot_list, spline_order, N_splines)
!
!open(file='r_list.dat', unit=1)
!open(file='func_list.dat', unit=2)
!do i=1, N_points 
!    write(1,*) r_list(i)
!    write(2,*) func_eval_arr(i)
!end do
!
!open(file='coeff_list.dat', unit=3)
!do i=1, N_splines
!    write(3,*) eigh_vec(i)
!end do

contains 

subroutine sort_eigh_vals(eigh_vals, E_max, N_splines, N_bound, N_cont, eigh_val_bound, eigh_val_cont, cont_index)
    integer, intent(in) :: N_splines 
    real(dp), intent(in) :: eigh_vals(N_splines), E_max 
    integer, intent(out) :: N_bound, N_cont 
    real(dp), allocatable, intent(out) :: eigh_val_bound(:), eigh_val_cont(:) 
    integer, allocatable, intent(out) :: cont_index(:)  ! Can be used to calculate density of states
    integer :: i, NB, NC 

    ! First count how many 
    N_bound = 0 
    N_cont = 0 
    do i=1, N_splines
        if (eigh_vals(i) < 0._dp) then 
            N_bound = N_bound + 1 
        else if (eigh_vals(i) <= E_max) then 
            N_cont = N_cont +1 
        end if 
    end do

    ! Reserve memory and fill in the arrays 
    NB = 1
    NC = 1 
    allocate(eigh_val_bound(N_bound), eigh_val_cont(N_cont), cont_index(N_cont))
    do i=1, N_splines
        if (eigh_vals(i) < 0._dp) then 
            eigh_val_bound(NB) = eigh_vals(i)
            NB = NB + 1 
        else if (eigh_vals(i) <= E_max) then 
            cont_index(NC) = i 
            eigh_val_cont(NC) = eigh_vals(i)
            NC = NC + 1
        end if 
    end do
end subroutine sort_eigh_vals

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

end program eigenstate 
