program eigenstate 
use utils 
use B_spline
use GaussLegendre 
use potential 
use linalg_routines
implicit none 

real(dp) :: r_max 
integer :: N_splines, N_points, spline_order, knot_init
real(dp), allocatable :: knot_list(:), r_list(:)
real(dp) :: x

integer :: N_gauss_leg 
real(dp), allocatable :: x_quad(:), w_quad(:), x_arr(:), w_arr(:)
real(dp), allocatable :: spline_arr(:,:), deriv_arr(:,:)
real(dp), allocatable :: H_mat(:,:), S_mat(:,:), H_temp(:,:), S_temp(:,:)  ! Hamiltonian and overlap matrix 
real(dp) :: H_ij, S_ij 
real(dp), allocatable :: V_arr(:), eigh_list(:), eigh_vec(:), func_eval_arr(:)
integer :: i, j, k, m, n, l, N_combi, counter


! First setup radial grid and knot sequence 
r_max = 100._dp 
N_points = 100
spline_order = 9
N_gauss_leg = spline_order + 2
l = 0  ! Angular momentum 

allocate(r_list(N_points))
call knot_squence_lin(knot_list, N_splines, r_list, N_points, spline_order, r_max)
!N_splines = N_splines -   ! Remove first and last to make functions go to 0 


! MATRIX ELEMENTS 
N_combi = 1
do i=2, spline_order 
    N_combi = N_combi * i 
end do 
allocate(H_temp(spline_order, N_splines), S_temp(spline_order, N_splines))
allocate(spline_arr(spline_order, N_gauss_leg), deriv_arr(spline_order, N_gauss_leg))
allocate(x_quad(N_gauss_leg), x_arr(N_gauss_leg), w_quad(N_gauss_leg), w_arr(N_gauss_leg), V_arr(N_gauss_leg))
!allocate(H_ij_arr(N_combi), S_ij_arr(N_combi))

H_temp = 0._dp 
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
            H_ij = l*(l+1._dp)/2._dp * sum(spline_arr(m,:)*spline_arr(n,:)/x_arr**2 * w_arr)  ! Centrifugal term
            H_ij = H_ij + sum(spline_arr(m,:)*spline_arr(n,:)*V_arr * w_arr)  ! Potential term 
            H_ij = H_ij + 0.5_dp * sum(deriv_arr(n,:) * deriv_arr(m,:) * w_arr)  ! Kinetic term 
            S_ij =  sum(spline_arr(m,:) * spline_arr(n,:) * w_arr)
    
            ! Now save in correct places 
            j = n + k-1  ! Index juggeling to save in band storage scheme, upper triangle, column in column, diagnoal in bottom
            i = spline_order + m - n 
            H_temp(i,j) = H_temp(i,j) + H_ij
            S_temp(i,j) = S_temp(i,j) + S_ij
            
            !write(*,*) i, j, m, n
            !if (H_mat(i,j) == 0) then 
            !    H_mat(i,j) = 1._dp * (m+k-1) + ((n+k-1)*0.1_dp)
            !end if 
        end do  
    end do

    knot_init = knot_init + 1  ! Update when moving along on radial grid 
end do 

! We now exclude first and last spline 
allocate(H_mat(spline_order, N_splines-2), S_mat(spline_order, N_splines-2))
H_mat = H_temp(:, 2:N_splines-1)
S_mat = S_temp(:, 2:N_splines-1)
do i=1, spline_order-1
    H_mat(spline_order-i, i) = 0._dp 
    S_mat(spline_order-i, i) = 0._dp 
end do
N_splines = N_splines-2  ! Reduce for the rest of the code 

!write(*,*) 
!do i=1, spline_order 
!    write(*,"(15(F8.3))") H_mat(i,:)
!end do

! Solve for the eigenvalues 
allocate(eigh_list(N_splines))
call generalized_symband_eigh(eigh_list, N_splines, spline_order, H_mat, S_mat)

write(*,*) 
do i=1, 2
    write(*,*) eigh_list(i)
end do

! Determine selected eigenvectors using inverse iteration 
allocate(eigh_vec(N_splines))
call inverse_iteration(eigh_vec, eigh_list(1), spline_order, N_splines, H_mat, S_mat)

! Save simualtion information 
open(file='data/sim_info.txt', unit=1)
write(1,*) N_splines+2
write(1,*) N_points
write(1,*) spline_order
write(1,*) size(knot_list)
write(1,*) N_gauss_leg
write(1,*) r_max 
close(1)

open(file='data/r_list.dat', unit=1, form='unformatted')
open(file='data/knot_list.dat', unit=2, form='unformatted')
open(file='data/init_state.dat', unit=3, form='unformatted')
write(1) r_list 
write(2) knot_list
write(3) eigh_vec  
close(1)
close(2)
close(3)


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

end program eigenstate 
