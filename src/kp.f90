MODULE KP
    USE UTILS
    USE CONSTANTS
    USE INDATA
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE GET_KP_HAMILTONIAN(psi, ham, potential, dx, dy, dz, Nx, Ny, Nz)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: psi(:,:,:,:)
        REAL*8, INTENT(OUT) :: ham(6,Nx, Ny, Nz)
        REAL*8, INTENT(IN) :: dx, dy, dz
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8 :: A, B, C, lx, ly, lz, mx, my, mz, n1, n2, n3
        INTEGER*4 :: i,j,k
        
       
        ham(:,:,:,:) = 0.0d0
        !!$omp parallel do collapse(3) schedule(static)
        DO k=2, Nz-1
            DO  j=2, Ny-1
                DO i=2, Nx-1
                    ham(1,i,j,k) = A * psi(1,i,j,k)&
                    -lx * (psi(1,i+1,j,k) + psi(1,i-1,j,k))&
                    -my * (psi(1,i,j+1,k) + psi(1,i,j-1,k))&
                    -mz * (psi(1,i,j,k+1) + psi(1,i,j,k-1))&
                    +n1 * (-psi(2,i+1,j+1,k)-psi(2,i-1,j-1,k)+psi(2,i+1,j-1,k)+psi(2,i-1,j+1,k))&
                    +n2 * (-psi(3,i+1,j,k+1)-psi(3,i-1,j,k-1)+psi(3,i+1,j,k-1)+psi(3,i-1,j,k+1))&
                    + potential(i,j,k) * psi(1,i,j,k)

                    ham(2,i,j,k) = B * psi(2,i,j,k)&
                    -mx * (psi(2,i+1,j,k) + psi(2,i-1,j,k))&
                    -ly * (psi(2,i,j+1,k) + psi(2,i,j-1,k))&
                    -lz * (psi(2,i,j,k+1) + psi(2,i,j,k-1))&
                    +n1 * (-psi(1,i+1,j+1,k)-psi(1,i-1,j-1,k)+psi(1,i+1,j-1,k)+psi(1,i-1,j+1,k))&
                    +n3 * (-psi(3,i,j+1,k+1)-psi(3,i,j-1,k-1)+psi(3,i,j-1,k+1)+psi(3,i,j+1,k-1))&
                    + potential(i,j,k) * psi(2,i,j,k)

                    ham(3,i,j,k) = C * psi(3,i,j,k)&
                    -mx * (psi(3,i+1,j,k) + psi(3,i-1,j,k))&
                    -my * (psi(3,i,j+1,k) + psi(3,i,j-1,k))&
                    -mz * (psi(3,i,j,k+1) + psi(3,i,j,k-1))&
                    +n2 * (-psi(1,i+1,j,k+1)-psi(1,i-1,j,k-1)+psi(1,i+1,j,k-1)+psi(1,i-1,j,k+1))&
                    +n3 * (-psi(2,i,j+1,k+1)-psi(2,i,j-1,k-1)+psi(2,i,j-1,k+1)+psi(2,i,j+1,k-1))&
                    + potential(i,j,k) * psi(3,i,j,k)


                    ham(4,i,j,k) = A * psi(4,i,j,k)&
                    -lx * (psi(4,i+1,j,k) + psi(4,i-1,j,k))&
                    -my * (psi(4,i,j+1,k) + psi(4,i,j-1,k))&
                    -mz * (psi(4,i,j,k+1) + psi(4,i,j,k-1))&
                    +n1 * (-psi(5,i+1,j+1,k)-psi(5,i-1,j-1,k)+psi(5,i+1,j-1,k)+psi(5,i-1,j+1,k))&
                    +n2 * (-psi(6,i+1,j,k+1)-psi(6,i-1,j,k-1)+psi(6,i+1,j,k-1)+psi(6,i-1,j,k+1))&
                    + potential(i,j,k) * psi(4,i,j,k)

                    ham(5,i,j,k) = B * psi(5,i,j,k)&
                    -mx * (psi(5,i+1,j,k) + psi(5,i-1,j,k))&
                    -ly * (psi(5,i,j+1,k) + psi(5,i,j-1,k))&
                    -lz * (psi(5,i,j,k+1) + psi(5,i,j,k-1))&
                    +n1 * (-psi(4,i+1,j+1,k)-psi(4,i-1,j-1,k)+psi(4,i+1,j-1,k)+psi(4,i-1,j+1,k))&
                    +n3 * (-psi(6,i,j+1,k+1)-psi(6,i,j-1,k-1)+psi(6,i,j-1,k+1)+psi(6,i,j+1,k-1))&
                    + potential(i,j,k) * psi(5,i,j,k)

                    ham(6,i,j,k) = C * psi(6,i,j,k)&
                    -mx * (psi(6,i+1,j,k) + psi(6,i-1,j,k))&
                    -my * (psi(6,i,j+1,k) + psi(6,i,j-1,k))&
                    -mz * (psi(6,i,j,k+1) + psi(6,i,j,k-1))&
                    +n2 * (-psi(4,i+1,j,k+1)-psi(4,i-1,j,k-1)+psi(4,i+1,j,k-1)+psi(4,i-1,j,k+1))&
                    +n3 * (-psi(5,i,j+1,k+1)-psi(5,i,j-1,k-1)+psi(5,i,j-1,k+1)+psi(5,i,j+1,k-1))&
                    + potential(i,j,k) * psi(6,i,j,k)
                END DO
            END DO
        END DO
        !!$omp end parallel do
        
    END SUBROUTINE GET_KP_HAMILTONIAN


    SUBROUTINE IMAGINARY_TIME_KP(potential, Nx, Ny, Nz, dx, dy, dz, dt, MAX_TIME, init_psi, final_psi, final_energy, tol)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: init_psi(:,:,:,:)
        REAL*8, INTENT(OUT) :: final_psi(norbital, Nx,Ny,Nz)
        REAL*8, INTENT(OUT) :: final_energy
        REAL*8, INTENT(IN) :: dx, dy, dz, dt, tol
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz, MAX_TIME
        REAL*8, ALLOCATABLE :: psi(:,:,:,:)
        REAL*8, ALLOCATABLE :: psi_new(:,:,:,:)
        REAL*8, ALLOCATABLE :: ham(:,:,:,:)
        INTEGER*4 :: i,j,k, iter
        REAL*8 :: energy, energy_old, norm
        REAL*8 :: kinetic_energy, potential_energy
        REAL*8 :: A, B, C, lx, ly, lz, mx, my, mz, n1, n2, n3, h1, h2, h3, h4, h5, h6
        A = 2*L/(dx**2) + M * (2/dy**2 + 2/dz**2)
        B = 2*L/(dy**2) + M * (2/dx**2 + 2/dz**2)
        C = 2*L/(dz**2) + M * (2/dx**2 + 2/dy**2)
        lx = L/(dx**2)
        ly = L/(dy**2)
        lz = L/(dz**2)
        mx = M/(dx**2)
        my = M/(dy**2)
        mz = M/(dz**2)
        n1 = N/(2*dx*dy)
        n2 = N/(2*dx*dz)
        n3 = N/(2*dz*dy)
        energy = 0.0d0
        kinetic_energy = 0.0d0
        potential_energy = 0.0d0
        print *, "Potential max (Ha):", maxval(potential)
        print *, "Potential min (Ha):", minval(potential)

        ! print*, "dt (au) =", dt
        ! print*, "dx (au) =", dx
        ! print*, "stability limit =", 1.0d0/(2.0d0*(3.0d0/(m1*dx**2)))
        if (dt > m1/(2.d0*(2.d0/dx**2 + 1.d0/dz**2))) then
            print*, "WARNING: dt too large, will blow up"
            stop
        endif

        ALLOCATE(psi(norbital, Nx, Ny, Nz))
        ALLOCATE(psi_new(norbital, Nx, Ny, Nz))
        ALLOCATE(ham(norbital, Nx, Ny, Nz))
        ham(:,:,:,:) = 0.0d0
        psi = init_psi
        energy_old = 1d100
        DO iter =1, MAX_TIME 
            ham(:,:,:,:) = 0.0d0
        !$omp parallel do collapse(3) schedule(static)
        DO k=2, Nz-1
            DO  j=2, Ny-1
                DO i=2, Nx-1
                    ham(1,i,j,k) = A * psi(1,i,j,k)&
                    -lx * (psi(1,i+1,j,k) + psi(1,i-1,j,k))&
                    -my * (psi(1,i,j+1,k) + psi(1,i,j-1,k))&
                    -mz * (psi(1,i,j,k+1) + psi(1,i,j,k-1))&
                    +n1 * (-psi(2,i+1,j+1,k)-psi(2,i-1,j-1,k)+psi(2,i+1,j-1,k)+psi(2,i-1,j+1,k))&
                    +n2 * (-psi(3,i+1,j,k+1)-psi(3,i-1,j,k-1)+psi(3,i+1,j,k-1)+psi(3,i-1,j,k+1))&
                    + potential(i,j,k) * psi(1,i,j,k)

                    ham(2,i,j,k) = B * psi(2,i,j,k)&
                    -mx * (psi(2,i+1,j,k) + psi(2,i-1,j,k))&
                    -ly * (psi(2,i,j+1,k) + psi(2,i,j-1,k))&
                    -lz * (psi(2,i,j,k+1) + psi(2,i,j,k-1))&
                    +n1 * (-psi(1,i+1,j+1,k)-psi(1,i-1,j-1,k)+psi(1,i+1,j-1,k)+psi(1,i-1,j+1,k))&
                    +n3 * (-psi(3,i,j+1,k+1)-psi(3,i,j-1,k-1)+psi(3,i,j-1,k+1)+psi(3,i,j+1,k-1))&
                    + potential(i,j,k) * psi(2,i,j,k)

                    ham(3,i,j,k) = C * psi(3,i,j,k)&
                    -mx * (psi(3,i+1,j,k) + psi(3,i-1,j,k))&
                    -my * (psi(3,i,j+1,k) + psi(3,i,j-1,k))&
                    -mz * (psi(3,i,j,k+1) + psi(3,i,j,k-1))&
                    +n2 * (-psi(1,i+1,j,k+1)-psi(1,i-1,j,k-1)+psi(1,i+1,j,k-1)+psi(1,i-1,j,k+1))&
                    +n3 * (-psi(2,i,j+1,k+1)-psi(2,i,j-1,k-1)+psi(2,i,j-1,k+1)+psi(2,i,j+1,k-1))&
                    + potential(i,j,k) * psi(3,i,j,k)


                    ham(4,i,j,k) = A * psi(4,i,j,k)&
                    -lx * (psi(4,i+1,j,k) + psi(4,i-1,j,k))&
                    -my * (psi(4,i,j+1,k) + psi(4,i,j-1,k))&
                    -mz * (psi(4,i,j,k+1) + psi(4,i,j,k-1))&
                    +n1 * (-psi(5,i+1,j+1,k)-psi(5,i-1,j-1,k)+psi(5,i+1,j-1,k)+psi(5,i-1,j+1,k))&
                    +n2 * (-psi(6,i+1,j,k+1)-psi(6,i-1,j,k-1)+psi(6,i+1,j,k-1)+psi(6,i-1,j,k+1))&
                    + potential(i,j,k) * psi(4,i,j,k)

                    ham(5,i,j,k) = B * psi(5,i,j,k)&
                    -mx * (psi(5,i+1,j,k) + psi(5,i-1,j,k))&
                    -ly * (psi(5,i,j+1,k) + psi(5,i,j-1,k))&
                    -lz * (psi(5,i,j,k+1) + psi(5,i,j,k-1))&
                    +n1 * (-psi(4,i+1,j+1,k)-psi(4,i-1,j-1,k)+psi(4,i+1,j-1,k)+psi(4,i-1,j+1,k))&
                    +n3 * (-psi(6,i,j+1,k+1)-psi(6,i,j-1,k-1)+psi(6,i,j-1,k+1)+psi(6,i,j+1,k-1))&
                    + potential(i,j,k) * psi(5,i,j,k)

                    ham(6,i,j,k) = C * psi(6,i,j,k)&
                    -mx * (psi(6,i+1,j,k) + psi(6,i-1,j,k))&
                    -my * (psi(6,i,j+1,k) + psi(6,i,j-1,k))&
                    -mz * (psi(6,i,j,k+1) + psi(6,i,j,k-1))&
                    +n2 * (-psi(4,i+1,j,k+1)-psi(4,i-1,j,k-1)+psi(4,i+1,j,k-1)+psi(4,i-1,j,k+1))&
                    +n3 * (-psi(5,i,j+1,k+1)-psi(5,i,j-1,k-1)+psi(5,i,j-1,k+1)+psi(5,i,j+1,k-1))&
                    + potential(i,j,k) * psi(6,i,j,k)
                END DO
            END DO
        END DO
        !$omp end parallel do
            psi_new = psi - dt*ham
            norm = 0.d0
            psi_new(:,1,:,:)  = 0.d0
            psi_new(:,Nx,:,:) = 0.d0
            psi_new(:,:,1,:)  = 0.d0
            psi_new(:,:,Ny,:) = 0.d0
            psi_new(:,:,:,1)  = 0.d0
            psi_new(:,:,:,Nz) = 0.d0
            energy = 0.d0
            ham(:,:,:,:) = 0.0d0
            norm = 0.d0

        !$omp parallel do collapse(3) schedule(static) reduction(+:norm)
        do k=1,Nz
            do j=1,Ny
                do i=1,Nx
                    norm = norm &
                        + psi_new(1,i,j,k)**2 &
                        + psi_new(2,i,j,k)**2 &
                        + psi_new(3,i,j,k)**2 &
                        + psi_new(4,i,j,k)**2 &
                        + psi_new(5,i,j,k)**2 &
                        + psi_new(6,i,j,k)**2
                end do
            end do
        end do
        !$omp end parallel do
            ! print*, "norm =", norm
            norm = sqrt(norm*dx*dx*dz)
            psi_new = psi_new/norm
            ! psi_new = psi_new / norm
            !$omp parallel do collapse(3) schedule(static) reduction(+:energy)&
            !$omp private(h1,h2,h3,h4,h5,h6)
            DO k=2, Nz-1
                DO  j=2, Ny-1
                    DO i=2, Nx-1
                        h1 = A * psi_new(1,i,j,k)&
                        -lx * (psi_new(1,i+1,j,k) + psi_new(1,i-1,j,k))&
                        -my * (psi_new(1,i,j+1,k) + psi_new(1,i,j-1,k))&
                        -mz * (psi_new(1,i,j,k+1) + psi_new(1,i,j,k-1))&
                        +n1 * (-psi_new(2,i+1,j+1,k)-psi_new(2,i-1,j-1,k)+psi_new(2,i+1,j-1,k)+psi_new(2,i-1,j+1,k))&
                        +n2 * (-psi_new(3,i+1,j,k+1)-psi_new(3,i-1,j,k-1)+psi_new(3,i+1,j,k-1)+psi_new(3,i-1,j,k+1))&
                        + potential(i,j,k) * psi_new(1,i,j,k)

                        h2 = B * psi_new(2,i,j,k)&
                        -mx * (psi_new(2,i+1,j,k) + psi_new(2,i-1,j,k))&
                        -ly * (psi_new(2,i,j+1,k) + psi_new(2,i,j-1,k))&
                        -lz * (psi_new(2,i,j,k+1) + psi_new(2,i,j,k-1))&
                        +n1 * (-psi_new(1,i+1,j+1,k)-psi_new(1,i-1,j-1,k)+psi_new(1,i+1,j-1,k)+psi_new(1,i-1,j+1,k))&
                        +n3 * (-psi_new(3,i,j+1,k+1)-psi_new(3,i,j-1,k-1)+psi_new(3,i,j-1,k+1)+psi_new(3,i,j+1,k-1))&
                        + potential(i,j,k) * psi_new(2,i,j,k)

                        h3 = C * psi_new(3,i,j,k)&
                        -mx * (psi_new(3,i+1,j,k) + psi_new(3,i-1,j,k))&
                        -my * (psi_new(3,i,j+1,k) + psi_new(3,i,j-1,k))&
                        -mz * (psi_new(3,i,j,k+1) + psi_new(3,i,j,k-1))&
                        +n2 * (-psi_new(1,i+1,j,k+1)-psi_new(1,i-1,j,k-1)+psi_new(1,i+1,j,k-1)+psi_new(1,i-1,j,k+1))&
                        +n3 * (-psi_new(2,i,j+1,k+1)-psi_new(2,i,j-1,k-1)+psi_new(2,i,j-1,k+1)+psi_new(2,i,j+1,k-1))&
                        + potential(i,j,k) * psi_new(3,i,j,k)


                        h4 = A * psi_new(4,i,j,k)&
                        -lx * (psi_new(4,i+1,j,k) + psi_new(4,i-1,j,k))&
                        -my * (psi_new(4,i,j+1,k) + psi_new(4,i,j-1,k))&
                        -mz * (psi_new(4,i,j,k+1) + psi_new(4,i,j,k-1))&
                        +n1 * (-psi_new(5,i+1,j+1,k)-psi_new(5,i-1,j-1,k)+psi_new(5,i+1,j-1,k)+psi_new(5,i-1,j+1,k))&
                        +n2 * (-psi_new(6,i+1,j,k+1)-psi_new(6,i-1,j,k-1)+psi_new(6,i+1,j,k-1)+psi_new(6,i-1,j,k+1))&
                        + potential(i,j,k) * psi_new(4,i,j,k)

                        h5 = B * psi_new(5,i,j,k)&
                        -mx * (psi_new(5,i+1,j,k) + psi_new(5,i-1,j,k))&
                        -ly * (psi_new(5,i,j+1,k) + psi_new(5,i,j-1,k))&
                        -lz * (psi_new(5,i,j,k+1) + psi_new(5,i,j,k-1))&
                        +n1 * (-psi_new(4,i+1,j+1,k)-psi_new(4,i-1,j-1,k)+psi_new(4,i+1,j-1,k)+psi_new(4,i-1,j+1,k))&
                        +n3 * (-psi_new(6,i,j+1,k+1)-psi_new(6,i,j-1,k-1)+psi_new(6,i,j-1,k+1)+psi_new(6,i,j+1,k-1))&
                        + potential(i,j,k) * psi_new(5,i,j,k)

                        h6 = C * psi_new(6,i,j,k)&
                        -mx * (psi_new(6,i+1,j,k) + psi_new(6,i-1,j,k))&
                        -my * (psi_new(6,i,j+1,k) + psi_new(6,i,j-1,k))&
                        -mz * (psi_new(6,i,j,k+1) + psi_new(6,i,j,k-1))&
                        +n2 * (-psi_new(4,i+1,j,k+1)-psi_new(4,i-1,j,k-1)+psi_new(4,i+1,j,k-1)+psi_new(4,i-1,j,k+1))&
                        +n3 * (-psi_new(5,i,j+1,k+1)-psi_new(5,i,j-1,k-1)+psi_new(5,i,j-1,k+1)+psi_new(5,i,j+1,k-1))&
                        + potential(i,j,k) * psi_new(6,i,j,k)

                        energy = energy + psi_new(1,i,j,k)*h1 + psi_new(2,i,j,k)*h2 &
                        + psi_new(3,i,j,k)*h3 + psi_new(4,i,j,k)*h4 &
                        + psi_new(5,i,j,k)*h5 + psi_new(6,i,j,k)*h6
                    END DO
                END DO
            END DO
            !$omp end parallel do

            energy = energy*dx*dx*dz
            if (abs((energy-energy_old)/feV2au) < tol) then
                print*, "Schrodinger converged after", iter, "iterations"
                print*, "Total energy (meV): ", energy/feV2au*1e3
                final_psi = psi_new
                final_energy = energy
                exit
            endif
            if (iter == MAX_TIME) then
                print *, "Warning: Schrodinger solver reached MAX_ITER without convergence. (max error)",&
                    abs(energy-energy_old)/feV2au
            end if
            energy_old = energy
                psi = psi_new
                ! PRINT*, "Iteration:", iter, "Energy:", energy/feV2au
        END DO
        final_psi = psi
        final_energy = energy
    END SUBROUTINE

    SUBROUTINE GET_Z_HAMILTONIAN(Ham_z, Nz, dz)
        IMPLICIT NONE
        REAL*8, INTENT(OUT) :: Ham_z(Nz,Nz)
        INTEGER*4, INTENT(IN) :: Nz
        REAL*8, INTENT(IN) :: dz
        INTEGER*4 :: i
        
        Ham_z(:,:) = 0.0d0
        DO i = 1, Nz
            Ham_z(i,i) = 2.0d0* L / (dz**2) + F_z * i * dz
            IF (i < Nz) THEN
                Ham_z(i,i+1) = -L / (dz**2)
                Ham_z(i+1,i) = -L / (dz**2)
            END IF         

        END DO

    END SUBROUTINE GET_Z_HAMILTONIAN

     SUBROUTINE GET_DFZ(Wavefunction_z, dfz, nstates, Nz, dz)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: Wavefunction_z(Nz, nstates)
        REAL*8, INTENT(OUT) :: dfz(Nz, nstates)
        INTEGER*4, INTENT(IN) :: nstates, Nz
        REAL*8, INTENT(IN) :: dz
        INTEGER*4 :: i, n

        dfz(:,:) = 0.0d0
        DO n = 1, nstates
            DO i = 2, Nz-1
                dfz(i,n) = (Wavefunction_z(i+1,n) - Wavefunction_z(i-1,n))/(2.0d0*dz)
            END DO
            dfz(1,n) = (Wavefunction_z(2,n) - Wavefunction_z(1,n))/(dz)
            dfz(Nz,n) = (Wavefunction_z(Nz,n) - Wavefunction_z(Nz-1,n))/(dz)
        END DO       
    END SUBROUTINE GET_DFZ

    SUBROUTINE GET_D2FZ(Wavefunction_z, dfz, d2fz, nstates, Nz, dz)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: Wavefunction_z(Nz, nstates), dfz(Nz, nstates)
        REAL*8, INTENT(OUT) :: d2fz(Nz, nstates)
        INTEGER*4, INTENT(IN) :: nstates, Nz
        REAL*8, INTENT(IN) :: dz
        INTEGER*4 :: i, n

        d2fz(:,:) = 0.0d0
        DO n = 1, nstates
            DO i = 2, Nz-1
                d2fz(i,n) = (Wavefunction_z(i-1,n) - 2.d0*Wavefunction_z(i,n) + Wavefunction_z(i+1,n))/(dz**2)
            END DO
            d2fz(1,n) = (dfz(2,n) - dfz(1,n))/(dz)
            d2fz(Nz,n) = (dfz(Nz,n) - dfz(Nz-1,n))/(dz)
        END DO       
    END SUBROUTINE GET_D2FZ

    SUBROUTINE GET_INIT_PSI_KP(psi, Nx, Ny, Nz, x0, y0, z0, sigma, dx, dz)

        IMPLICIT NONE

        REAL*8, INTENT(OUT) :: psi(norbital, Nx, Ny, Nz)
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: x0, y0, z0, sigma, dx, dz

        INTEGER*4 :: i, j, k
        REAL*8 :: x, y, z, norm

        psi = 0.0d0

        ! Initial Gaussian in component 1
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx

                    x = (i-1)*dx
                    y = (j-1)*dx
                    z = (k-1)*dz

                    psi(1,i,j,k) = initial_psi(x,y,z,x0,y0,z0,sigma)

                END DO
            END DO
        END DO

        ! Dirichlet boundary conditions
        psi(:,1,:,:)  = 0.0d0
        psi(:,Nx,:,:) = 0.0d0
        psi(:,:,1,:)  = 0.0d0
        psi(:,:,Ny,:) = 0.0d0
        psi(:,:,:,1)  = 0.0d0
        psi(:,:,:,Nz) = 0.0d0

        ! Norm
        norm = 0.0d0

        DO k = 2, Nz-1
            DO j = 2, Ny-1
                DO i = 2, Nx-1
                    norm = norm + DOT_PRODUCT(psi(:,i,j,k), &
                                            psi(:,i,j,k))
                END DO
            END DO
        END DO

        norm = sqrt(norm * dx * dx * dz)

        ! Normalization
        psi = psi / norm

    END SUBROUTINE GET_INIT_PSI_KP

    SUBROUTINE GET_DENSITY_KP(density, psi, Nx, Ny, Nz)
        IMPLICIT NONE

        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: psi(:,:,:,:)
        REAL*8, INTENT(OUT) :: density(Nx, Ny, Nz)

        INTEGER*4 :: i, j, k, norb

        density = 0.d0

        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx
                    DO norb = 1, 6
                        density(i,j,k) = density(i,j,k) &
                            + abs(psi(norb,i,j,k))**2
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE GET_DENSITY_KP


END MODULE