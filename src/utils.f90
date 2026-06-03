MODULE UTILS
    USE CONSTANTS
    IMPLICIT NONE
    CONTAINS

    function initial_psi(x,y,z,x0,y0,z0,sigma) result(val)
        implicit none

        double precision, intent(in) :: x, y, z, x0, y0, z0, sigma
        double precision :: val

        val = 1.0d0 / ((2.0d0*pi)**(3.0d0/2.0d0) * sigma**3) * &
        exp(-((x-x0)**2 + (y-y0)**2 + (z-z0)**2) / (2.0d0*sigma**2))

    end function initial_psi


END MODULE UTILS