module utils 
implicit none 
private

! Different kinds of floating point numbers
integer, parameter :: sp = kind(1.0)
integer, parameter :: dp = kind(1.0d0) 

! Some constants 
real(dp), parameter :: pi = 3.1415926535897932_dp  ! Do I need the dp here? 

! Define the working precision for this module! 
integer, parameter :: wp = dp  

! Public constants and kinds 
public sp, dp, pi 

! Interface for overloading 
interface trapz 
    module procedure trapz_real, trapz_complex
end interface trapz 

interface trapz_2D
    module procedure trapz_2D_real, trapz_2D_complex
end interface trapz_2D

! Public routines
public linspace, trapz, trapz_2D


contains 
    ! Just a standard linspace subroutine, nothing fancy 
    subroutine linspace(arr, start, end, N)
        integer, intent(in) :: N 
        real(wp), intent(in) :: start, end 
        real(wp), intent(out) :: arr(N)
        real(wp) :: step 
        integer :: i 

        step = (end - start) / (N-1)

        do i = 1, N
            arr(i) = start + (i-1)*step 
        end do 
    end subroutine linspace


    ! A simple trapezoidal rule integration function 
    function trapz_real(func_vec, step_size)
        real(wp), intent(in) :: step_size  ! The step size of the grid used for function evaluation 
        real(wp), intent(in) :: func_vec(:)  ! The function to integrate evaluated at evenly spaced points x
        real(wp) :: trapz_real
        integer :: N  ! To hold the size of the array
        
        N = size(func_vec)
        trapz_real = 0.5_dp * step_size * (func_vec(1) + func_vec(N)) + step_size * sum(func_vec(2:N-1))
    end function trapz_real

    ! A simple trapezoidal rule integration function 
    function trapz_complex(func_vec, step_size)
        real(wp), intent(in) :: step_size  ! The step size of the grid used for function evaluation 
        complex(wp), intent(in) :: func_vec(:)  ! The function to integrate evaluated at evenly spaced points x
        complex(wp) :: trapz_complex
        integer :: N  ! To hold the size of the array
        
        N = size(func_vec)
        trapz_complex = 0.5_dp * step_size * (func_vec(1) + func_vec(N)) + step_size * sum(func_vec(2:N-1))
    end function trapz_complex

    ! Function to perform a 2D integral over an real array(x,y)
    function trapz_2D_real(arr, step_x, step_y, Nx, Ny) result(res)
        integer :: Nx, Ny
        real(wp) :: arr(Nx, Ny), step_x, step_y, res
        real(wp) :: temp_arr(Ny)
        integer :: i

        ! Perform all the x-integrals 
        do i=1, Ny
            temp_arr(i) = trapz(arr(:,i), step_x)
        end do

        ! Now perform the y-integral 
        res = trapz(temp_arr, step_y)
    end function trapz_2D_real

    !Function to perform a 2D integral over an real array(x,y)
    function trapz_2D_complex(arr, step_x, step_y, Nx, Ny) result(res)
        integer :: Nx, Ny
        complex(wp) :: arr(Nx, Ny), res
        real(wp) :: step_x, step_y
        complex(wp) :: temp_arr(Ny)
        integer :: i

        ! Perform all the x-integrals 
        do i=1, Ny
            temp_arr(i) = trapz(arr(:,i), step_x)
        end do

        ! Now perform the y-integral 
        res = trapz(temp_arr, step_y)
    end function trapz_2D_complex
    
end module utils 