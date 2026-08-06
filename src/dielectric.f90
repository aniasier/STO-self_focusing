MODULE DIELECTRIC
    USE CONSTANTS
    USE UTILS
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

SUBROUTINE GET_EPSILON(potential_1d, eps0, &
                       nx, ny, &
                       nz_1d, dz_1d, &
                       nz_3d, dz_3d, &
                       epsilon)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx, ny
    INTEGER, INTENT(IN) :: nz_1d, nz_3d

    REAL*8, INTENT(IN) :: eps0
    REAL*8, INTENT(IN) :: dz_1d, dz_3d
    REAL*8, INTENT(IN) :: potential_1d(nz_1d)

    REAL*8, INTENT(OUT) :: epsilon(nx,ny,nz_3d)

    REAL*8 :: potential_3d(nz_3d)
    REAL*8 :: electric_field(nz_3d)

    INTEGER :: i,j,k

    !---------------------------------------------------------
    ! interpolacja potencjału
    !---------------------------------------------------------

    CALL INTERPOLATE_1D(potential_1d, nz_1d, dz_1d, &
                        potential_3d, nz_3d, dz_3d)

    !---------------------------------------------------------
    ! liczenie pola elektrycznego na siatce 3D
    !---------------------------------------------------------

    electric_field(1) = -(potential_3d(2)-potential_3d(1))/dz_3d

    DO k=2,nz_3d-1

        electric_field(k)=-(potential_3d(k+1)-potential_3d(k-1)) &
                          /(2.d0*dz_3d)

    END DO

    electric_field(nz_3d)=-(potential_3d(nz_3d)- &
                            potential_3d(nz_3d-1))/dz_3d

    !---------------------------------------------------------
    ! epsilon(E)
    !---------------------------------------------------------

    DO k=1,nz_3d

        DO i=1,nx
            DO j=1,ny

                epsilon(i,j,k)=permitivity(eps0,electric_field(k))

            END DO
        END DO

    END DO

END SUBROUTINE GET_EPSILON

END MODULE DIELECTRIC