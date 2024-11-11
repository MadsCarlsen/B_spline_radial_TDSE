module potential 
use utils 
implicit none 
private 

public :: V_pot, vector_pot

contains 

function V_pot(r) 
    real(dp) :: r(:)
    real(dp) :: V_pot(size(r)) 
    V_pot = -1._dp / r !- exp(-(r - 5._dp)**2/0.1_dp)
end function V_pot

function vector_pot(t, omega, Nc, rtUp, CEP) 
    real(dp) vector_pot, t, omega, rtUp, CEP 
    integer Nc 

    if (t>=0. .and. t<=2._dp*pi * Nc/omega) then 
        vector_pot = 2._dp*rtUp * sin(omega*t / (2._dp*Nc))**2 * cos(omega*t + CEP)
    else 
        vector_pot = 0._dp 
    end if 
end function vector_pot 


end module potential 