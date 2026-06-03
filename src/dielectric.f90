MODULE DIELECTRIC
    IMPLICIT NONE
    CONTAINS

    function permitivity(eps_0,electric_field)
    implicit double precision (a-h,o-z)
    fm2au=18.897261*1e9
    feV2au=0.03674932587
    A=4.097*1e-5
    B=4.907*1e-10*(fm2au/feV2au)
    permitivity=eps_0+1.0/(A+B*dabs(electric_field))

    !Inny sposob
     !B=2.55*1e4
     !E0=8.22*1e4*feV2au/fm2au
     !permitivity=1+B/(1.0+(electric_field/E0)**2)**(1.0/3.0)
    
    return
end

SUBROUTINE GET_EPSILON(electric_field, eps0, nx, ny, nz, epsilon)
    IMPLICIT NONE
    REAL*8, INTENT (IN) :: eps0
    INTEGER*4, INTENT(IN) :: nx, ny, nz
    REAL*8, INTENT (IN) :: electric_field(nz)
    REAL*8, INTENT(OUT) :: epsilon(nx, ny, nz)
    INTEGER*4 :: i, j, k
    REAL*8 :: val
    DO k=1, Nz
        val = permitivity(eps0, electric_field(k))
        DO i=1, nx
            DO j=1, ny
                epsilon(i,j,k) = val
            END DO
        END DO
    END DO
END SUBROUTINE GET_EPSILON

END MODULE DIELECTRIC