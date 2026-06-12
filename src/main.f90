PROGRAM MAIN
    USE CONSTANTS
    USE Poisson_Solver_Mod
    USE DIELECTRIC 
    USE UTILS
    USE WRITERS
    USE SCHRODINGER
    IMPLICIT NONE
    REAL*8 :: n0_trapped, L_trapped, eps_0, dz, dx, m1, m2, thickness
    INTEGER*4 :: nz, nx, ny
    REAL*8, ALLOCATABLE :: charge_trapped(:)
    REAL*8, ALLOCATABLE :: charge_trapped3D(:, :, :)
    REAL*8, ALLOCATABLE :: potential_z(:)
    REAL*8, ALLOCATABLE :: electric_field(:)
    REAL*8, ALLOCATABLE :: electric_field_new(:)
    REAL*8, ALLOCATABLE :: epsilon(:, :, :)
    REAL*8, ALLOCATABLE :: potential(:, :, :)
    REAL*8, ALLOCATABLE :: potential_eps0(:, :, :)
    REAL*8, ALLOCATABLE :: density(:, :, :)
    REAL*8, ALLOCATABLE :: density_full(:, :, :)
    REAL*8 :: x0, y0, z0 ! gauss centering
    REAL*8 :: sigma ! size fo the initial gaussian
    REAL*8 :: alfa ! relaxation parameter for poisson solver
    REAL*8 :: tol
    INTEGER*4 :: MAX_ITER
    INTEGER*4 :: i, j, k, iz
    REAL*8 :: z
    tol = 1.0e-6
    MAX_ITER = 100000
    alfa = 1.5

    n0_trapped = 5.0*1.0e13*fne2D2au
    L_trapped=15*fnm2au
    eps_0=100
    ! thickness=12.0*fnm2au
    dz = 0.1*fnm2au
    dx = 0.1*fnm2au
    nz = 100!ceiling(thickness/dz)
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
    ALLOCATE(potential_eps0(nx, ny,nz))
    ALLOCATE(density(nx, ny,nz))
    ALLOCATE(density_full(nx, ny,nz))
    ALLOCATE(charge_trapped3D(nx, ny, nz))
    potential(:,:,:) =0.0d0
    potential_eps0(:,:,:) =0.0d0

    x0 = (nx-1)*dx/2.0d0
    y0 = (ny-1)*dx/2.0d0
    z0 = (nz-1)*dz/2.0d0
    sigma = 2*fnm2au
    ! stage 1: z direction
    CALL POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0, nz, dz, charge_trapped, electric_field, potential_z)
    ! CALL POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)
    ! stage 2: dielectric
    CALL GET_EPSILON(electric_field, eps_0, nx, ny, nz, epsilon)
    CALL GET_CHARGE_TRAPPED3D(charge_trapped3D, charge_trapped, nx, ny, nz)
    CALL GET_DENSITY(density, nx, ny, nz, x0, y0, z0, sigma, dx)
    CALL WRITE_DENSITY_2D_XY(density, nx, ny, nz, dx, 'data/density.dat')
    ! stage 3: poisson in 3d with changing dielectric function
    CALL Poisson_epsilon(potential, density, epsilon, alfa, nx, ny, nz, dx, tol, MAX_ITER, charge_trapped3D)
    CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz, dx, 'data/potential.dat')
    CALL WRITE_POTENTIAL_CROSS_SECTION(potential, nx, ny, nz, dx, 'data/potential_cross_section.dat')
    density_full = 0.0d0
    DO i=1, Nx
        DO j=1, Ny
            DO k=1, Nz
                density_full(i,j,k)=density(i,j,k)+ charge_trapped3D(i,j,k)
            END DO
        END DO
    END DO
    ! i need density_full because Poisson solver doesnt take charge trapped
    
    DO i=1, Nx
        DO j=1, Ny
            DO k=1, Nz
                potential(i,j,k)=potential(i,j,k)+ potential_z(k)
            END DO
        END DO
    END DO

    ! stage 4: poisson with epsilon NOT changing - epsilon 0 or epsilon R ????
    ! CALL Poisson(potential_eps0, density_full, eps_0, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER)
    ! subtracting -> only the influence of the changing eps at STO interface
    potential = potential - potential_eps0
    CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz, dx, 'data/potential_eps0.dat')

    ! state 5: imaginary time method for schrodinger equation
    potential = 0.0d0
    print*, "Expected E =", (3.14159265d0**2/(2.0d0*m1) * &
    (2.0d0/((Nx-1)*dx)**2 + 1.0d0/((Nz-1)*dz)**2))/ feV2au

    CALL IMAGINARY_TIME(potential, Nx, Ny, Nz, dx, dz, m1, m2, x0, y0, z0, sigma)


    DEALLOCATE(charge_trapped)
    DEALLOCATE(electric_field)
    DEALLOCATE(electric_field_new)
    DEALLOCATE(potential_z)
    DEALLOCATE(epsilon)
    DEALLOCATE(potential)
    DEALLOCATE(potential_eps0)
    DEALLOCATE(density)
    DEALLOCATE(density_full)
    DEALLOCATE(charge_trapped3D)


END PROGRAM MAIN
