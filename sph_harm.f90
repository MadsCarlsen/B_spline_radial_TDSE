module spherical_harmonics 
use utils 
implicit none 
private 

! Define working precision 
integer, parameter :: wp = dp 

public :: sph_harm, assoc_legendre

contains 

! Function to determine the spherical harmonics 
function sph_harm(cos_theta, phi, l, m)
    integer :: l,m 
    real(wp) :: cos_theta, phi 
    complex(wp) :: sph_harm 
    real(wp) :: sqrt_fac, plm 
    integer :: abs_m

    ! Calcualate legendre poly for positive m and check sign afterwards 
    abs_m = abs(m)
    sqrt_fac = sqrt((2._wp*l+1._wp)/(4._wp*pi) * gamma(l-abs_m+1._wp)/gamma(l+abs_m+1._wp))
    plm = assoc_legendre(l, abs_m, cos_theta)
    
    if (m>=0) then 
        sph_harm = sqrt_fac * plm * exp(cmplx(0._wp, m*phi, wp))
    else 
        sph_harm = (-1)**abs_m * sqrt_fac * plm * exp(cmplx(0._wp, m*phi, wp))  
    end if    
end function sph_harm 


! Function to recursively evaluate the associated legendre polynomials
function assoc_legendre(l, m, x) result(plm)
    integer :: l, m
    real(wp) :: x, plm 
    real(wp) :: pmm, pmm1, pli, fact, sqrt_fact 
    integer :: li, i
    
    ! First  calculate the Pmm legendre poly, used to start recursion 
    pmm = 1._wp  ! Value of p00 
    if(m.gt.0) then
        sqrt_fact = sqrt((1._wp-x)*(1._wp+x))
        fact=1._wp
        do i=1, m
            pmm = -pmm*fact*sqrt_fact
            fact = fact + 2._wp
        enddo 
    end if
    
    ! Check if l=m, then this is the final value 
    if (l == m) then
        plm=pmm
    else 
        ! If not calculate Pm(m+1) also needed to start recursion 
        pmm1 = x*(2._wp*m+1._wp)*pmm

        ! Check again if we are done 
        if (l == m+1) then 
            plm = pmm1 
        else 
            ! Now we really have to do recursion 
            do li = m+2, l
                pli = (x*(2._wp*li-1._wp)*pmm1 - (li+m-1._wp)*pmm) / (li-m)
                pmm = pmm1 
                pmm1 = pli 
            end do 
            plm = pli 
        end if 
    end if 
end function assoc_legendre

end module spherical_harmonics 