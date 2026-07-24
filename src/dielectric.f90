MODULE DIELECTRIC
    USE CONSTANTS
    IMPLICIT NONE
    CONTAINS

    function permitivity(eps_0,electric_field)
    implicit double precision (a-h,o-z)
    A=4.097*1e-5
    B=4.907*1e-10*(fm2au/feV2au)
    permitivity=eps_0+1.0/(A+B*abs(electric_field))

    !Inny sposob
     !B=2.55*1e4
     !E0=8.22*1e4*feV2au/fm2au
     !permitivity=1+B/(1.0+(electric_field/E0)**2)**(1.0/3.0)
    
    return
end

SUBROUTINE GET_EPSILON(electric_field, eps0, nx, ny, nz_1D, nz_3d, epsilon)
    IMPLICIT NONE
    REAL*8, INTENT (IN) :: eps0
    INTEGER*4, INTENT(IN) :: nx, ny, nz_1d, nz_3d
    REAL*8, INTENT (IN) :: electric_field(nz_1d)
    REAL*8, INTENT(OUT) :: epsilon(nx, ny, nz_3d)
    INTEGER*4 :: i, j, k,  mult_factor
    REAL*8 :: val
    mult_factor = (nz_1d-1)/(nz_3d-1)
    DO k = 1, nz_3d
        val = permitivity(eps0, electric_field(1 + (k-1)*mult_factor))
        DO i = 1, nx
            DO j = 1, ny
                epsilon(i,j,k) = val
            END DO
        END DO
    END DO
END SUBROUTINE GET_EPSILON

END MODULE DIELECTRIC