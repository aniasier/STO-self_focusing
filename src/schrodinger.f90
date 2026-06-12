MODULE SCHRODINGER
    USE UTILS
    USE CONSTANTS
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE IMAGINARY_TIME(potential, Nx, Ny, Nz, dx, dz, m1, m2, x0, y0, z0, sigma)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: dx, dz, m1, m2, x0, y0, z0, sigma
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, ALLOCATABLE :: psi(:,:,:)
        REAL*8, ALLOCATABLE :: psi_new(:,:,:)
        REAL*8, ALLOCATABLE :: ham(:,:,:)
        REAL*8 :: val
        INTEGER*4 :: i,j,k, iter, MAX_TIME
        REAL*8 :: x, y, z, energy, energy_old, tol, norm, dt

        MAX_TIME = 1000
        dt = 0.1*fns2au
        energy = 0.0d0


        ALLOCATE(psi(Nx, Ny, Nz))
        ALLOCATE(psi_new(Nx, Ny, Nz))
        ALLOCATE(ham(Nx, Ny, Nz))
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx
                    x=(i-1)*dx
                    y=(j-1)*dx
                    z=(k-1)*dx
                    psi(i,j,k)=initial_psi(x,y,z,x0,y0,z0,sigma)
                enddo
            enddo
        enddo

        ! calculating the norm
        val=0.0
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx         
                    val = val+abs(psi(i,j,k))**2*dx*dx*dx
                enddo
            enddo
        enddo
        ! normalization
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx         
                    psi(i,j,k) = psi(i,j,k)/sqrt(val)
                enddo
            enddo
        enddo

        
        energy_old = 1.d99
        tol = 1.d-10
        DO iter =1, MAX_TIME
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

                    norm = norm + abs(psi_new(i,j,k)**2)

                    end do
                end do
            end do

            norm = sqrt(norm*dx*dx*dz)

            psi_new = psi_new/norm
            energy = 0.d0

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
                exit
            endif

            energy_old = energy
            psi = psi_new
            PRINT*, "Iteration:", iter, "Energy:", energy
        END DO

    END SUBROUTINE


END MODULE SCHRODINGER