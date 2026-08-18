MODULE SCHRODINGER
    USE UTILS
    USE CONSTANTS
    USE INDATA
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE IMAGINARY_TIME(potential, Nx, Ny, Nz, dx, dz, dt, MAX_TIME, m1, m2, init_psi, final_psi, final_energy, tol)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: init_psi(:,:,:)
        REAL*8, INTENT(OUT) :: final_psi(Nx,Ny,Nz)
        REAL*8, INTENT(OUT) :: final_energy
        REAL*8, INTENT(IN) :: dx, dz, m1, m2, dt, tol
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz, MAX_TIME
        REAL*8, ALLOCATABLE :: psi(:,:,:)
        REAL*8, ALLOCATABLE :: psi_new(:,:,:)
        REAL*8, ALLOCATABLE :: ham(:,:,:)
        INTEGER*4 :: i,j,k, iter
        REAL*8 :: energy, energy_old, norm
        REAL*8 :: kinetic_energy, potential_energy
        REAL*8 :: ham_temp

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

        ALLOCATE(psi(Nx, Ny, Nz))
        ALLOCATE(psi_new(Nx, Ny, Nz))
        ALLOCATE(ham(Nx, Ny, Nz))
        ham(:,:,:) = 0.0d0
        psi = init_psi
       
        energy_old = 1.d99
        ! tol = 1.d-6
        DO iter =1, MAX_TIME
            ham(:,:,:) = 0.0d0
            !!$omp parallel do collapse(3) schedule(static)
            DO k=2, Nz-1
                DO  j=2, Ny-1
                    DO i=2, Nx-1
                        ham(i,j,k) = -1/(2.0d0*m1)*(( psi(i+1,j,k)+psi(i-1,j,k) + psi(i,j+1,k)+psi(i,j-1,k)-&
                        4.d0*psi(i,j,k))/(dx**2) + (psi(i,j,k+1)+psi(i,j,k-1) - 2.d0*psi(i,j,k) )/dz**2) &
                        + potential(i,j,k) * psi(i,j,k)
                    END DO
                END DO
            END DO
            !!$omp end parallel do
            psi_new = psi - dt*ham
            norm = 0.d0

            do k=1,Nz
                do j=1,Ny
                    do i=1,Nx

                    norm = norm + (psi_new(i,j,k))**2

                    end do
                end do
            end do

            psi_new(1,:,:)  = 0.d0
            psi_new(Nx,:,:) = 0.d0
            psi_new(:,1,:)  = 0.d0
            psi_new(:,Ny,:) = 0.d0
            psi_new(:,:,1)  = 0.d0
            psi_new(:,:,Nz) = 0.d0
            energy = 0.d0
            ham(:,:,:) = 0.0d0

            ! print*, "norm =", norm
            norm = sqrt(norm*dx*dx*dz)

            psi_new = psi_new/norm
            ! psi_new = psi_new / norm

            DO i=2, Nx-1
                DO  j=2, Ny-1
                    DO k=2, Nz-1
                        ham_temp = -1/(2.0d0*m1)*((psi_new(i+1,j,k)+psi_new(i-1,j,k) + psi_new(i,j+1,k)+psi_new(i,j-1,k)-&
                        4.d0*psi_new(i,j,k))/(dx**2) + (psi_new(i,j,k+1)+psi_new(i,j,k-1) - 2.d0*psi_new(i,j,k) )/dz**2) &
                        + potential(i,j,k) * psi_new(i,j,k)
                        energy = energy + psi_new(i,j,k)*ham_temp
                    END DO
                END DO
            END DO


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

    SUBROUTINE GET_KP_HAMILTONIAN()
        IMPLICIT NONE

        ! (i, j, k)
        fijk(1,1) = -2*L/(dx**2) + M * (-2/dy**2 + -2/dz**2)

        fijk(2,2) = -2*L/(dy**2) + M * (-2/dx**2 + -2/dz**2)

        fijk(3,3) = -2*L/(dz**2) + M * (-2/dx**2 + -2/dy**2)
        
        ! (i+1, j, k)
        fiplusjk(1,1) = L/(dx**2)

        fiplusjk(2,2) = M/(dx**2)

        fiplusjk(3,3) = M/(dx**2)

        ! (i-1, j, k)
        fiminusjk(1,1) = L/(dx**2)

        fiminusjk(2,2) = M/(dx**2)

        fiminusjk(3,3) = M/(dx**2)

        ! (i, j+1, k)
        fijplusk(1,1) = M/(dy**2)

        fijplusk(2,2) = L/(dy**2)

        fijplusk(3,3) = M/(dy**2)

        ! (i, j-1, k)
        fijminusk(1,1) = M/(dy**2)

        fijminusk(2,2) = L/(dy**2)

        fijminusk(3,3) = M/(dy**2)

        ! (i, j, k+1)
        fijkplus(1,1) = M/(dz**2)
        
        fijkplus(2,2) = M/(dz**2)

        fijkplus(3,3) = L/(dz**2)

        ! (i, j, k-1)
        fijkminus(1,1) = M/(dz**2)

        fijkminus(2,2) = M/(dz**2)

        fijkminus(3,3) = L/(dz**2)

        

END MODULE SCHRODINGER