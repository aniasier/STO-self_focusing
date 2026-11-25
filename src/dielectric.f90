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

SUBROUTINE UPDATE_EPSILON_3D(potential_el, E_trapped_3d, eps0, &
                            Nx, Ny, Nz, dx, dz, epsilon, damp)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Ny, Nz
    REAL*8,  INTENT(IN) :: potential_el(Nx, Ny, Nz)
    REAL*8,  INTENT(IN) :: E_trapped_3d(Nz)  ! pole tła od uwięzionych ładunków
    REAL*8,  INTENT(IN) :: eps0, dx, dz, damp
    REAL*8,  INTENT(INOUT) :: epsilon(Nx, Ny, Nz)

    INTEGER :: i, j, k
    REAL*8  :: Ex, Ey, Ez, E_mag, eps_calc

    DO k = 2, Nz - 1
        DO j = 2, Ny - 1
            DO i = 2, Nx - 1
                ! 1. Gradient potencjału (pole elektryczne elektronu)
                Ex = (potential_el(i+1,j,k) - potential_el(i-1,j,k)) / (2.0d0 * dx)
                Ey = (potential_el(i,j+1,k) - potential_el(i,j-1,k)) / (2.0d0 * dx)
                Ez = (potential_el(i,j,k+1) - potential_el(i,j,k-1)) / (2.0d0 * dz)

                ! 2. Dodajemy pole tła w osi Z
                Ez = Ez + E_trapped_3d(k)

                ! 3. Moduł całkowitego pola |E| w 3D
                E_mag = SQRT(Ex**2 + Ey**2 + Ez**2)

                ! 4. Wyliczenie nowej przenikalności
                eps_calc = permitivity(eps0, E_mag)

                ! 5. Miksowanie (damping), żeby pętla SCF była stabilna i nie oscylowała
                epsilon(i,j,k) = (1.0d0 - damp) * epsilon(i,j,k) + damp * eps_calc
            END DO
        END DO
    END DO

    ! Proste warunki brzegowe dla epsilon na ściankach pudła
    epsilon(1,:,:)  = epsilon(2,:,:)
    epsilon(Nx,:,:) = epsilon(Nx-1,:,:)
    epsilon(:,1,:)  = epsilon(:,2,:)
    epsilon(:,Ny,:) = epsilon(:,Ny-1,:)
    epsilon(:,:,1)  = epsilon(:,:,2)
    epsilon(:,:,Nz) = epsilon(:,:,Nz-1)

END SUBROUTINE UPDATE_EPSILON_3D



END MODULE DIELECTRIC