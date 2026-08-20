PROGRAM MAIN
    USE CONSTANTS
    USE INDATA
    USE Poisson_Solver
    USE DIELECTRIC 
    USE UTILS
    USE WRITERS
    USE SCHRODINGER
    USE KP
    IMPLICIT NONE
    REAL*8 :: eps_0
    REAL*8, ALLOCATABLE :: charge_trapped(:)
    REAL*8, ALLOCATABLE :: charge_trapped3D(:, :, :)
    REAL*8, ALLOCATABLE :: potential_z(:)
    REAL*8, ALLOCATABLE :: electric_field(:)
    ! REAL*8, ALLOCATABLE :: electric_field_new(:)
    REAL*8, ALLOCATABLE :: epsilon(:, :, :)
    REAL*8, ALLOCATABLE :: potential(:, :, :)
    REAL*8, ALLOCATABLE :: potential_eps0(:, :, :)
    REAL*8, ALLOCATABLE :: density(:, :, :)
    REAL*8, ALLOCATABLE :: density_full(:, :, :)
    REAL*8, ALLOCATABLE :: init_psi(:,:,:,:)
    REAL*8, ALLOCATABLE :: final_psi(:,:,:,:)
    ! for hamiltoanian
    REAL*8, ALLOCATABLE :: Ham_z(:, :)
    REAL*8, ALLOCATABLE :: Energies_z(:)
    REAL*8, ALLOCATABLE :: Wavefunction_z(:, :)
    REAL*8, ALLOCATABLE :: dfz(:,:)
    REAL*8, ALLOCATABLE :: d2fz(:,:)


    REAL*8 :: x0, y0, z0 ! gauss centering
    INTEGER*4 :: i, j, k, iz, iter
    REAL*8 :: z, V0, sigma_v
    REAL*8 :: beta
    REAL*8 :: energy, energy_old, eps_local
    CHARACTER(LEN=50) :: filename
    integer :: clock_start, clock_end, clock_rate
    real*8 :: elapsed_time, deps
    call system_clock(clock_start, clock_rate)
    beta = 0.5d0
    eps_0=100 ! <- do wzoru na permittivity wyraz wolny
    ! thickness=12.0*fnm2au
    CALL GET_INDATA("./data/input.nml")

    ALLOCATE(charge_trapped(nz_1d))
    ALLOCATE(electric_field(nz_1d))
    ! ALLOCATE(electric_field_new(nz))
    ALLOCATE(potential_z(nz_1d))
    ALLOCATE(epsilon(nx, ny,nz_3d))
    ALLOCATE(potential(nx, ny,nz_3d))
    ALLOCATE(potential_eps0(nx, ny,nz_3d))
    ALLOCATE(density(nx, ny,nz_3d))
    ALLOCATE(density_full(nx, ny,nz_3d))
    ALLOCATE(charge_trapped3D(nx, ny, nz_3d))
    ALLOCATE(init_psi(norbital,nx, ny, nz_3d))
    ALLOCATE(final_psi(norbital,nx, ny, nz_3d))
    ALLOCATE(Ham_z(nz_3d, nz_3d))
    ALLOCATE(Energies_z(nz_3d))
    ALLOCATE(Wavefunction_z(nz_3d, nstate_1))
    ALLOCATE(dfz(Nz_3d, nstate_1))
    ALLOCATE(d2fz(Nz_3d, nstate_1))
    
    potential_eps0(:,:,:) =0.0d0
    x0 = (nx-1)*dx/2.0d0
    y0 = (ny-1)*dx/2.0d0
    z0 = (z0_indx)*dz_3d
    ! stage 1: z direction
    CALL POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0, nz_1d, dz_1d, charge_trapped, electric_field, potential_z)
    ! CALL POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)
    ! stage 2: dielectric
    CALL GET_EPSILON(potential_z, eps_0, nx, ny, nz_1d, dz_1d, nz_3d, dz_3d, epsilon)
    CALL GET_CHARGE_TRAPPED3D(charge_trapped3D, charge_trapped, nx, ny, nz_3d)
    CALL GET_INIT_PSI_KP(init_psi, Nx, Ny, Nz_3d, x0, y0, z0, sigma, dx, dz_3d)
    CALL GET_DENSITY_KP(density, init_psi, nx, ny, nz_3d)
    CALL WRITE_DENSITY_2D_XY(density, nx, ny, nz_3d, dx,dz_3d, './data/density_init.dat')
    ! CALL WRITE_DENSITY_2D_XY_SLICE(density, nx, ny, nz, dx, dz, z0_indx, './data/density_init_slice.dat')
    CALL WRITE_DENSITY_CROSS_SECTION(density, Nx, Ny, Nz_3d, dz_3d, './data/density_init_crossection.dat')
    CALL WRITE_DENSITY_CROSS_SECTION_X(density, Nx, Ny, Nz_3d, dx, z0_indx, './data/density_init_crossection_x.dat')
    CALL WRITE_DENSITY_CROSS_SECTION_Y(density, Nx, Ny, Nz_3d, dx, z0_indx, './data/density_init_crossection_y.dat')
    ! stage 3: poisson in 3d with changing dielectric function
    energy_old = 1.d99

    ! PREP FOR SCHORDINGER
    ! CALL GET_Z_HAMILTONIAN(Ham_z, nz_3d, dz_3d)
    ! CALL solve_eigenproblem(Ham_z, Energies_z, nz_3d)
    ! Wavefunction_z = Ham_z(:, 1:nstate_1)
    ! CALL WRITE_Z_WAVEFUNCTION(Ham_z, nz_3d, nstate_1, nz_3d, dz_3d, "./data/Psi_z")
    ! CALL WRITE_ENERGIES(Energies_z, nstate_1, "./data/Energies_z.dat")
    ! CALL GET_DFZ(Wavefunction_z, dfz, nstate_1, nz_3d, dz_3d)
    ! CALL GET_D2FZ(Wavefunction_z,dfz, d2fz, nstate_1, Nz_3d, dz_3d)
    ! charge_trapped3D(:,:,:) = 0.0d0
    DO iter = 1, MAX_ITER_SCF
        potential(:,:,:) =0.0d0
        ! PRINT*, "SCF ITERATION:", iter
        CALL Poisson_epsilon_no_charge(potential, density, epsilon, alfa, nx, ny, nz_3d, dx, dz_3d, tol, MAX_ITER)
        DO k=1, Nz_3d
            DO j=1, Ny
                DO i=1, Nx
                    potential(i,j,k)=potential(i,j,k)-potential_z(k)
                END DO
            END DO
        END DO
        eps_local = epsilon(ceiling((nx-1)/2.0d0), ceiling((ny-1)/2.0d0), z0_indx)  
        ! stage 4: poisson with epsilon NOT changing
        CALL Poisson(potential_eps0, density, eps_local, alfa, Nx, Ny, Nz_3d, dx, dz_3d, tol, MAX_ITER)
        potential = potential - potential_eps0
        PRINT*, "start schrodinger"
        CALL IMAGINARY_TIME_KP(potential, Nx, Ny, Nz_3d, dx,dx, dz_3d, dt, MAX_TIME, init_psi, final_psi, energy, tol)
        CALL GET_DENSITY_KP(density, final_psi, nx, ny, nz_3d)
        WRITE(filename, '(A,I0,A)') './data/density3D_', iter, '.dat'
        CALL WRITE_DENSITY_2D_XY(density, Nx, Ny, Nz_3d, dx, dz_3d, filename)

        WRITE(filename, '(A,I0,A)') './data/density_final_slice_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY_SLICE(density, nx, ny, nz_3d, dx, dz_3d, z0_indx, filename)

        WRITE(filename, '(A,I0,A)') './data/density_final_crossection_', iter, '.dat'
        CALL WRITE_DENSITY_CROSS_SECTION(density, Nx, Ny, Nz_3d, dz_3d, filename)

        WRITE(filename, '(A,I0,A)') './data/density_final_crossection_x_', iter, '.dat'
        CALL WRITE_DENSITY_CROSS_SECTION_X(density, Nx, Ny, Nz_3d, dx, z0_indx, filename)

        WRITE(filename, '(A,I0,A)') './data/density_final_crossection_y_', iter, '.dat'
        CALL WRITE_DENSITY_CROSS_SECTION_Y(density, Nx, Ny, Nz_3d, dx, z0_indx, filename)

        WRITE(filename, '(A,I0,A)') './data/potential3D_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY(potential, Nx, Ny, Nz_3d, dx, dz_3d, filename)

        WRITE(filename, '(A,I0,A)') './data/potential_final_slice_', iter, '.dat'
        CALL WRITE_POTENTIAL_2D_XY_SLICE(potential, nx, ny, nz_3d, dx, dz_3d, z0_indx, filename)

        WRITE(filename, '(A,I0,A)') './data/potential_final_crossection_', iter, '.dat'
        CALL WRITE_POTENTIAL_CROSS_SECTION(potential, Nx, Ny, Nz_3d, dz_3d, filename)

        WRITE(filename, '(A,I0,A)') './data/potential_final_crossection_x_', iter, '.dat'
        CALL WRITE_POTENTIAL_CROSS_SECTION_X(potential, Nx, Ny, Nz_3d, dx, z0_indx, filename)

        WRITE(filename, '(A,I0,A)') './data/potential_final_crossection_y_', iter, '.dat'
        CALL WRITE_POTENTIAL_CROSS_SECTION_Y(potential, Nx, Ny, Nz_3d, dx, z0_indx, filename)
        if (abs(energy-energy_old) < tol_scf) then
                print*, "Converged after", iter, "iterations"
                exit
            endif
        energy_old = energy
        init_psi=final_psi
        
        CALL POISSON_ZDIRECTION_PLUS_POTENTIAL(n0_trapped, L_trapped, eps_0,  nz_1D, nx, ny, dz_1D, dx, charge_trapped, &
        electric_field, density, nz_3d, dz_3D, potential_z, x0, y0)
        ! potential_z = beta*potential_z_new + (1.d0-beta)*potential_z_old
        CALL GET_EPSILON(potential_z, eps_0, nx, ny, nz_1d, dz_1d, nz_3d, dz_3d, epsilon)
    END DO
    CALL WRITE_POTENTIAL_2D_XY(potential, nx, ny, nz_3d, dx, dz_3d, './data/potential_final.dat')
    ! CALL WRITE_POTENTIAL_2D_XY_SLICE(potential, nx, ny, nz, dx, dz, z0_indx, './data/potential_final_slice.dat')
    CALL WRITE_POTENTIAL_CROSS_SECTION(potential, Nx, Ny, Nz_3d, dz_3d, './data/potential_final_crossection.dat')
    CALL WRITE_POTENTIAL_CROSS_SECTION_X(potential, Nx, Ny, Nz_3d, dx, z0_indx, './data/potential_final_crossection_x.dat')
    CALL WRITE_POTENTIAL_CROSS_SECTION_Y(potential, Nx, Ny, Nz_3d, dx, z0_indx, './data/potential_final_crossection_y.dat')

    CALL WRITE_DENSITY_2D_XY(density, Nx, Ny, Nz_3d, dx, dz_3d, './data/density_final.dat')
    ! CALL WRITE_DENSITY_2D_XY_SLICE(density, nx, ny, nz, dx, dz, z0_indx, './data/density_final_slice.dat')
    CALL WRITE_DENSITY_CROSS_SECTION(density, Nx, Ny, Nz_3d, dz_3d, './data/density_final_crossection.dat')
    CALL WRITE_DENSITY_CROSS_SECTION_X(density, Nx, Ny, Nz_3d, dx, z0_indx, './data/density_final_crossection_x.dat')
    CALL WRITE_DENSITY_CROSS_SECTION_Y(density, Nx, Ny, Nz_3d, dx, z0_indx, './data/density_final_crossection_y.dat')


    PRINT*, "ENERGY (meV):", energy/feV2au*1e3
    call system_clock(clock_end)
    elapsed_time = real(clock_end-clock_start,8)/real(clock_rate,8)
    print *, "Time [s] =", elapsed_time


    DEALLOCATE(charge_trapped)
    DEALLOCATE(electric_field)
    ! DEALLOCATE(electric_field_new)
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
