module linalg_routines
use utils 
implicit none 
private 

public :: generalized_symband_eigh, inverse_iteration, band_sym_matrix_vec_prod, & 
          band_antisym_matrix_vec_prod, solve_complex_lin_system

contains 

! Routine to solve real symmetrix generalized eigenvalue problem Az = lambda Bz, 
! with A, B real symmetric and banded, and B positive definit.
! Assume A and B are given in band-stroage form, diagnoal in lowest row 
subroutine generalized_symband_eigh(eigh_list, N, N_band, A, B)
    integer, intent(in) :: N, N_band  ! N_band is nr of diagnoal + superdiagonals  
    real(dp), intent(in) :: A(N_band, N), B(N_band, N) 
    real(dp), intent(out) :: eigh_list(N)
    real(dp) :: A_cal(N_band, N), B_cal(N_band, N), X(1, N), WORK(2*N)
    real(dp) :: D(N), E(N), Q(1, N), WORK1(N)
    integer :: info, LDX, LDQ
   
    A_cal = A 
    B_cal = B 
    LDX = 1 

    ! First perform split-Cholesky factorization using DPBSTF 
    call DPBSTF('U', N, N_band-1, B_cal, N_band, info)
    if (info /= 0) then 
        write(*,*) 'Error in split-Cholesky factorization!'
        write(*,*) info 
        stop 
    end if 

    ! Now convert generalized band problem into standard band problem using DSBGST
    call DSBGST('N', 'U', N, N_band-1, N_band-1, A_cal, N_band, B_cal, N_band, X, LDX, WORK, info)
    if (info /= 0) then 
        write(*,*) 'Error in conversion to regular eigenvalue problem!'
        write(*,*) info 
        stop 
    end if 

    ! Now we can treat it like a standard real band-symmetric eigenvalue problem
    ! Start by reducing A_cal to tridiagnoal real form using DSBTRD 
    LDQ = 1 
    call DSBTRD('N', 'U', N, N_band-1, A_cal, N_band, D, E, Q, LDQ, WORK1, info)
    if (info /= 0) then 
        write(*,*) 'Error in conversion to tridiagnoal matrix!'
        write(*,*) info 
        stop 
    end if 

    ! Now we can finally obtain the eigenvalues. Several methods, for now use DSTERF 
    call DSTERF(N, D, E, info)
    if (info /= 0) then 
        write(*,*) 'Error in obtaining eigenvalues!'
        write(*,*) info 
        stop 
    end if 
    eigh_list = D
end subroutine generalized_symband_eigh

! Inverse iteration routine to determine a eigenvector using the corresponding eigenvalue.
subroutine inverse_iteration(eigh_vec, eigh_val, spline_order, N_splines, H_mat, S_mat, eps, max_iter)
    integer, intent(in) :: spline_order, N_splines, max_iter
    real(dp), intent(in) :: eigh_val, H_mat(spline_order, N_splines), S_mat(spline_order, N_splines), eps
    real(dp), intent(out) :: eigh_vec(N_splines)

    integer :: LDAB, i, j, counter 
    real(dp) :: H_full(3*(spline_order-1)+1, N_splines)  ! Size needed for DGBTRF 
    real(dp) :: S_full(3*(spline_order-1)+1, N_splines)
    real(dp) :: sub_arr(2*spline_order-1, N_splines)  ! Vector to contain S * solution vector 
    real(dp) :: norm, old_vec(N_splines), sign_new, sign_old 
    
    ! First we need to build the full banded matrices needed for the LU-factorization 
    ! Our matrix should start in (spline_order-1)+1 until 3*(spline_order-1) + 1
    LDAB = 3*(spline_order-1)+1
    H_full(spline_order:2*spline_order-1, :) = H_mat 
    S_full(spline_order:2*spline_order-1, :) = S_mat 
    
    ! Build the lower part 
    counter = 2 
    do i = 1, spline_order-1 
        do j = 1, N_splines - i 
            H_full(2*spline_order-1 + i, j) = H_full(2*spline_order-1 - i, j+i)
            S_full(2*spline_order-1 + i, j) = S_full(2*spline_order-1 - i, j+i)
        end do 
    end do 

    ! Now we perform inverse iteration 
    eigh_vec = 1._dp ! Some initial guess

    ! Normalize 
    norm = calculate_norm_real(eigh_vec, S_mat, spline_order, N_splines)
    eigh_vec = eigh_vec / norm

    ! Prepare matrices / vectors  
    H_full = H_full - eigh_val * S_full  ! Overwrite to matrix used in inverse iteration 
    sub_arr = S_full(spline_order:, :)
    old_vec = 0._dp 
    sign_old = 1._dp 

    do i=1, max_iter 
        ! Calculate matrix-vector product of S_mat on guess vec.
        call banded_matrix_vec_prod(sub_arr, eigh_vec, spline_order, N_splines)
        
        ! Solve lin system 
        call solve_lin_system(H_full, eigh_vec, LDAB, N_splines, spline_order)

        ! Normalize the solution - done by using the overlap matrix? 
        norm = calculate_norm_real(eigh_vec, S_mat, spline_order, N_splines) 
        eigh_vec = eigh_vec / norm 
       
        ! Check phase is consistent (NB: overall phase is arbitary but must agree to compare!)
        sign_new = sign(1.0_dp, eigh_vec(1))
        if (sign_new /= sign_old) then 
            old_vec = -1._dp * old_vec 
        end if 

        ! Check convergence 
        if (norm2(eigh_vec - old_vec) < eps) then 
            return 
        end if 

        old_vec = eigh_vec 
        sign_old = sign_new 
    end do 

    write(*,*) 'Inverse iteration did not converge! Returning unconverged result.'
end subroutine inverse_iteration 

function calculate_norm_real(vec, S_mat, N_band, N_splines) result(norm)
    integer, intent(in) :: N_band, N_splines 
    real(dp) :: vec(N_splines), S_mat(N_band, N_splines)
    real(dp) :: norm 
    real(dp) :: cal_vec(N_splines)

    ! First apply S-matrix to state 
    call band_sym_matrix_vec_prod_real(S_mat, vec, cal_vec, N_splines, N_band)

    ! Take inner product with state vector 
    norm = sqrt(sum(vec * cal_vec))
end function calculate_norm_real

! Calculates the norm of the eigenstate expanded in the B-spline basis. Assumes real state. 
!function calculate_norm(eigh_vec, S_mat, spline_order, N_splines) result(norm)
!    integer, intent(in) :: N_splines, spline_order 
!    real(dp), intent(in) :: eigh_vec(N_splines), S_mat(spline_order, N_splines)
!    real(dp) :: norm
!    integer :: i, j
!
!    norm  = 0._dp 
!
!    do i=1,N_splines
!        do j=i, min(N_splines, i+spline_order-1)
!            ! Handle diagonal, off-diagonal gets mult 2 
!            if (i == j) then 
!                norm = norm + eigh_vec(i) * eigh_vec(j) * S_mat(spline_order+i-j, j)
!            else 
!                norm = norm + 2._dp * eigh_vec(i) * eigh_vec(j) * S_mat(spline_order+i-j, j)
!            end if 
!        end do 
!    end do
!    norm = sqrt(norm)
!end function calculate_norm 

! Routine for calculating the matrix vector product, when the matrix is band symmetric. Real vector.
subroutine band_sym_matrix_vec_prod_real(matrix, vec, res, N_vec, N_band)
    integer, intent(in) :: N_vec, N_band  ! N_vec is size of diagnoal, N_band is nr. of superdiagnoals + diagonal
    real(dp), intent(in) :: matrix(N_band, N_vec) 
    real(dp), intent(in) :: vec(N_vec)
    real(dp), intent(out) :: res(N_vec)
    integer :: i, j, k, Nj
    real(dp) :: sum_i 
    
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
end subroutine band_sym_matrix_vec_prod_real

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

! Routine for calculating the matrix vector product, when the matrix is band symmetric
subroutine band_antisym_matrix_vec_prod(matrix, vec, res, N_vec, N_band)
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
            sum_i = sum_i - matrix(Nj+i, j) * vec(i)  ! Parts for lower triagnle, OBS minus for anti-symmetric
        end do 
        res(j) = res(j) + sum_i  ! No diagnoal, but add terms from lower triangle
    end do
end subroutine band_antisym_matrix_vec_prod 

! Calculate the matrix-vector product of a banded matrix. Not assuming symmetric mat
subroutine banded_matrix_vec_prod(mat, vec, N_band, N_vec)
    integer, intent(in) :: N_band, N_vec  ! N_band is nr of super/subdiagonals + the diagonal
    real(dp), intent(in) :: mat(2*N_band-1, N_vec)
    real(dp), intent(inout) :: vec(N_vec)
    real(dp) :: cal_vec(N_vec), sum_i

    integer :: i, j 

    ! mat_ij is stored in mat(N_band + i-j, j) for max(1,j-N_band+1) <= i <= min(m, j+N_band-1) where m is size of mat 
    ! Loop over matrix product 
    do i=1, N_vec  
        sum_i = 0._dp 
        ! Only loop over the band entires 
        do j=max(1, i-N_band+1), min(N_vec, i+N_band-1)
            sum_i = sum_i + vec(j) * mat(N_band + i-j, j)
        end do
        cal_vec(i) = sum_i  
    end do 
    vec = cal_vec 
end subroutine banded_matrix_vec_prod

! Solve linear system Ax=b through a LU-factorization.
subroutine solve_lin_system(A, vec, LDAB, N_splines, spline_order)
    integer, intent(in) :: LDAB, N_splines, spline_order
    real(dp), intent(inout) :: vec(N_splines)
    real(dp), intent(in) :: A(LDAB, N_splines)
    integer :: IPIV(N_splines), INFO 
    real(dp) :: A_cal(LDAB, N_splines), cal_vec(N_splines,1)
    integer :: NRHS 
    
    NRHS = 1 
    A_cal = A 
    cal_vec(:,1) = vec

    ! First perform LU facorization using DGBTRF
    call DGBTRF(N_splines, N_splines, spline_order-1, spline_order-1, A_cal, LDAB, IPIV, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Error when prefrom LU factorization!'
        write(*,*) INFO 
        stop 
    end if 

    ! Now solve the system using DGBTRS
    call DGBTRS('N', N_splines, spline_order-1, spline_order-1, NRHS, A_cal, LDAB, IPIV, cal_vec, N_splines, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Error when solving linear system!'
        write(*,*) INFO 
        stop 
    end if 
    vec = cal_vec(:,1)
end subroutine solve_lin_system 

! Solve linear system Ax=b through a LU-factorization.
subroutine solve_complex_lin_system(A, vec, LDAB, N_splines, spline_order)
    integer, intent(in) :: LDAB, N_splines, spline_order
    complex(dp), intent(inout) :: vec(N_splines)
    complex(dp), intent(in) :: A(LDAB, N_splines)
    integer :: IPIV(N_splines), INFO 
    complex(dp) :: A_cal(LDAB, N_splines), cal_vec(N_splines,1)
    integer :: NRHS 
    
    NRHS = 1 
    A_cal = A 
    cal_vec(:,1) = vec

    ! First perform LU facorization using DGBTRF
    call ZGBTRF(N_splines, N_splines, spline_order-1, spline_order-1, A_cal, LDAB, IPIV, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Error when prefrom LU factorization!'
        write(*,*) INFO 
        stop 
    end if 

    ! Now solve the system using DGBTRS
    call ZGBTRS('N', N_splines, spline_order-1, spline_order-1, NRHS, A_cal, LDAB, IPIV, cal_vec, N_splines, INFO)
    if (INFO /= 0) then 
        write(*,*) 'Error when solving linear system!'
        write(*,*) INFO 
        stop 
    end if 
    vec = cal_vec(:,1)
end subroutine solve_complex_lin_system 

end module linalg_routines