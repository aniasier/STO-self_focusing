PROGRAM MAIN
    USE CONSTANTS
    USE INDATA
    USE Poisson_Solver
    USE DIELECTRIC 
    USE UTILS
    USE WRITERS
    USE SCHRODINGER
    IMPLICIT NONE
    REAL*8 :: eps_0
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
    REAL*8, ALLOCATABLE :: init_psi(:,:,:)
    REAL*8, ALLOCATABLE :: final_psi(:,:,:)
    REAL*8 :: x0, y0, z0 ! gauss centering
    INTEGER*4 :: i, j, k, iz, iter
    REAL*8 :: z
    REAL*8 :: energy, energy_old
    CHARACTER(LEN=50) :: filename

    eps_0=100 ! <- do wzoru na permittivity wyraz wolny
    ! thickness=12.0*fnm2au
    CALL GET_INDATA("input.nml")

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
    ALLOCATE(init_psi(nx, ny, nz))
    ALLOCATE(final_psi(nx, ny, nz))
    
    potential_eps0(:,:,:) =0.0d0

    x0 = (nx-1)*dx/2.0d0
    y0 = (ny-1)*dx/2.0d0
    z0 = (nz-1)*dz/2.0d0
    ! stage 1: z direction
    CALL POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0, nz, dz, charge_trapped, electric_field, potential_z)
    ! CALL POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)
    ! stage 2: dielectric
    CALL GET_EPSILON(electric_field, eps_0, nx, ny, nz, epsilon)
    CALL GET_CHARGE_TRAPPED3D(charge_trapped3D, charge_trapped, nx, ny, nz)
    CALL GET_INIT_PSI(init_psi, Nx, Ny, Nz, x0, y0, z0, sigma, dx, dz)
    CALL GET_DENSITY(density, init_psi, nx, ny, nz)
    CALL WRITE_DENSITY_2D_XY(density, nx, ny, nz, dx,dz, 'data/density.dat')
    ! stage 3: poisson in 3d with changing dielectric function
    energy_old = 1.d99
    ! charge_trapped3D(:,:,:) = 0.0d0
    DO iter = 1, MAX_ITER_SCF
        potential(:,:,:) =0.0d0
        PRINT*, "SCF ITERATION:", iter
        CALL Poisson_epsilon_no_charge(potential, density, epsilon, alfa, nx, ny, nz, dx, dz, tol, MAX_ITER)
        WRITE(filename, '(A,I0,A)') 'data/potential_nocharge_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz, dx, filename)
        CALL WRITE_POTENTIAL_CROSS_SECTION(potential, nx, ny, nz, dx, 'data/potential_cross_section.dat')
        ! density_full = 0.0d0
        ! DO i=1, Nx
        !     DO j=1, Ny
        !         DO k=1, Nz
        !             density_full(i,j,k)=density(i,j,k)+ charge_trapped3D(i,j,k)
        !         END DO
        !     END DO
        ! END DO
        ! i need density_full because Poisson solver doesnt take charge trapped
        
        DO i=1, Nx
            DO j=1, Ny
                DO k=1, Nz
                    potential(i,j,k)=-potential(i,j,k)-potential_z(k)
                END DO
            END DO
        END DO

        WRITE(filename, '(A,I0,A)') 'data/potential_plus_z_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz, dx, filename)

        ! stage 4: poisson with epsilon NOT changing
        CALL Poisson(potential_eps0, density, eps_0, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER)
        ! subtracting -> only the influence of the changing eps at STO interface

        WRITE(filename, '(A,I0,A)') 'data/potential_eps0', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY(potential_eps0, nx, ny, nz, dx, filename)
        potential = potential + potential_eps0
        WRITE(filename, '(A,I0,A)') 'data/potential_final_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz, dx, filename)

        ! state 5: imaginary time method for schrodinger equation
        ! potential = 0.0d0
        ! print*, "Expected E =", (3.14159265d0**2/(2.0d0*m1) * &
        ! (2.0d0/((Nx-1)*dx)**2 + 1.0d0/((Nz-1)*dz)**2))/ feV2au

        CALL IMAGINARY_TIME(potential, Nx, Ny, Nz, dx, dz, dt, MAX_TIME, m1, m2, init_psi, final_psi, energy)
        CALL GET_DENSITY(density, final_psi, nx, ny, nz)
        WRITE(filename, '(A,I0,A)') 'data/density3D_', iter, '.dat'
        CALL WRITE_DENSITY_2D_XY(density, Nx, Ny, Nz, dx, dz, filename)

        if (abs(energy-energy_old) < tol_scf) then
                print*, "Converged after", iter, "iterations"
                exit
            endif
        energy_old = energy
        init_psi=final_psi
    END DO

    PRINT*, "ENERGY (meV):", energy/feV2au


    DEALLOCATE(charge_trapped)
    DEALLOCATE(electric_field)
    DEALLOCATE(electric_field_new)
    DEALLOCATE(potential_z)
    DEALLOCATE(epsilon)
    DEALLOCATE(potential)
    DEALLOCATE(potential_eps0)
    DEALLOCATE(density)
    DEALLOCATE(init_psi)
    DEALLOCATE(final_psi)
    DEALLOCATE(density_full)
    DEALLOCATE(charge_trapped3D)


END PROGRAM MAIN
