MODULE SCHRODINGER
    USE UTILS
    USE CONSTANTS
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE IMAGINARY_TIME(potential, Nx, Ny, Nz, dx, dz, dt, MAX_TIME, m1, m2, init_psi, final_psi, final_energy)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: init_psi(:,:,:)
        REAL*8, INTENT(OUT) :: final_psi(Nx,Ny,Nz)
        REAL*8, INTENT(OUT) :: final_energy
        REAL*8, INTENT(IN) :: dx, dz, m1, m2, dt
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz, MAX_TIME
        REAL*8, ALLOCATABLE :: psi(:,:,:)
        REAL*8, ALLOCATABLE :: psi_new(:,:,:)
        REAL*8, ALLOCATABLE :: ham(:,:,:)
        INTEGER*4 :: i,j,k, iter
        REAL*8 :: energy, energy_old, tol, norm

        energy = 0.0d0

        ! print*, "dt (au) =", dt
        ! print*, "dx (au) =", dx
        ! print*, "stability limit =", 1.0d0/(2.0d0*(3.0d0/(m1*dx**2)))
        if (dt > 1.0d0/(2.0d0*(3.0d0/(m1*dx**2)))) then
            print*, "WARNING: dt too large, will blow up"
            stop
        endif

        ALLOCATE(psi(Nx, Ny, Nz))
        ALLOCATE(psi_new(Nx, Ny, Nz))
        ALLOCATE(ham(Nx, Ny, Nz))
        ham(:,:,:) = 0.0d0
        psi = init_psi
       
        energy_old = 1.d99
        tol = 1.d-6
        DO iter =1, MAX_TIME
            ham(:,:,:) = 0.0d0
            DO i=2, Nx-1
                DO  j=2, Ny-1
                    DO k=2, Nz-1
                        ham(i,j,k) = -1/(2.0d0*m1)*(( psi(i+1,j,k)+psi(i-1,j,k) + psi(i,j+1,k)+psi(i,j-1,k)-&
                        4.d0*psi(i,j,k))/(dx**2) + (psi(i,j,k+1)+psi(i,j,k-1) - 2.d0*psi(i,j,k) )/dz**2) &
                        + potential(i,j,k) * psi(i,j,k)
                    END DO
                END DO
            END DO
            psi_new = psi - dt*ham
            norm = 0.d0

            do i=1,Nx
                do j=1,Ny
                    do k=1,Nz

                    norm = norm + abs(psi_new(i,j,k))**2

                    end do
                end do
            end do
            ! print*, "norm =", norm
            norm = sqrt(norm*dx*dx*dz)

            psi_new = psi_new/norm
            ! psi_new = psi_new / norm

            psi_new(1,:,:)  = 0.d0
            psi_new(Nx,:,:) = 0.d0
            psi_new(:,1,:)  = 0.d0
            psi_new(:,Ny,:) = 0.d0
            psi_new(:,:,1)  = 0.d0
            psi_new(:,:,Nz) = 0.d0
            energy = 0.d0
            ham(:,:,:) = 0.0d0

            DO i=2, Nx-1
                DO  j=2, Ny-1
                    DO k=2, Nz-1
                        ham(i,j,k) = -1/(2.0d0*m1)*((psi_new(i+1,j,k)+psi_new(i-1,j,k) + psi_new(i,j+1,k)+psi_new(i,j-1,k)-&
                        4.d0*psi_new(i,j,k))/(dx**2) + (psi_new(i,j,k+1)+psi_new(i,j,k-1) - 2.d0*psi_new(i,j,k) )/dz**2) &
                        + potential(i,j,k) * psi_new(i,j,k)
                    END DO
                END DO
            END DO

            do i=2,Nx-1
            do j=2,Ny-1
            do k=2,Nz-1

            energy = energy + psi_new(i,j,k)*ham(i,j,k)

            end do
            end do
            end do

            energy = energy*dx*dx*dz
            if (abs(energy-energy_old) < tol) then
                print*, "Converged after", iter, "iterations"
                print*, "Energy (eV): ", energy/feV2au
                final_psi = psi_new
                final_energy = energy
                exit
            endif

            energy_old = energy
            psi = psi_new
            ! PRINT*, "Iteration:", iter, "Energy:", energy/feV2au
            if (iter == MAX_TIME) then
                print *, "Warning: Schrodinger solver reached MAX_ITER without convergence. (max error)", abs(energy-energy_old)
            end if
        END DO
        final_psi = psi
        final_energy = energy
    END SUBROUTINE


END MODULE SCHRODINGER