PROGRAM MAIN
    USE CONSTANTS
    USE Poisson_Solver_Mod
    USE DIELECTRIC 
    USE UTILS
    IMPLICIT NONE
    REAL*8 :: n0_trapped, L_trapped, eps_0, dz, m1, m2, thickness
    INTEGER*4 :: nz, nx, ny
    REAL*8, ALLOCATABLE :: charge_trapped(:)
    REAL*8, ALLOCATABLE :: potential_z(:)
    REAL*8, ALLOCATABLE :: electric_field(:)
    REAL*8, ALLOCATABLE :: electric_field_new(:)
    REAL*8, ALLOCATABLE :: epsilon(:, :, :)
    REAL*8, ALLOCATABLE :: potential(:, :, :)

    n0_trapped = 5.0*1.0e13*fne2D2au
    L_trapped=15*fnm2au
    eps_0=100
    thickness=100.0*fnm2au
    dz = 0.1*fnm2au
    nz = ceiling(thickness/dz)
    nx =100
    ny = 100
    m1=0.2
    m2=3.5

    ALLOCATE(charge_trapped(nz))
    ALLOCATE(electric_field(nz))
    ALLOCATE(electric_field_new(nz))
    ALLOCATE(potential_z(nz))
    ALLOCATE(epsilon(nx, ny,nz))
    ALLOCATE(potential(nx, ny,nz))

    potential(:,:,:) =0.0d0

    CALL POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0, nz, dz, charge_trapped, electric_field, potential_z)
    ! CALL POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)
    CALL GET_EPSILON(electric_field, eps_0, nx, ny, nz, epsilon)

    ! CALL Poisson_epsilon(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER, charge_trapped)

    DEALLOCATE(charge_trapped)
    DEALLOCATE(electric_field)
    DEALLOCATE(electric_field_new)
    DEALLOCATE(potential_z)
    DEALLOCATE(epsilon)


END PROGRAM MAIN
