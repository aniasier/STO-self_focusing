MODULE SCHRODINGER
    USE UTILS
    USE CONSTANTS
    USE INDATA
    IMPLICIT NONE
    INTEGER, PARAMETER :: norb = 6
    TYPE :: KP_STENCIL

        REAL*8 :: fijk(norb,norb)

        REAL*8 :: fiplusjk(norb,norb)
        REAL*8 :: fiminusjk(norb,norb)

        REAL*8 :: fijplusk(norb,norb)
        REAL*8 :: fijminusk(norb,norb)

        REAL*8 :: fijkplus(norb,norb)
        REAL*8 :: fijkminus(norb,norb)

        REAL*8 :: fiplusjplusk(norb,norb)
        REAL*8 :: fiplusjminusk(norb,norb)
        REAL*8 :: fiminusjplusk(norb,norb)
        REAL*8 :: fiminusjminusk(norb,norb)

        REAL*8 :: fiplusjkplus(norb,norb)
        REAL*8 :: fiplusjkminus(norb,norb)
        REAL*8 :: fiminusjkplus(norb,norb)
        REAL*8 :: fiminusjkminus(norb,norb)

        REAL*8 :: fijpluskplus(norb,norb)
        REAL*8 :: fijpluskminus(norb,norb)
        REAL*8 :: fijminuskplus(norb,norb)
        REAL*8 :: fijminuskminus(norb,norb)

    END TYPE KP_STENCIL
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

    SUBROUTINE GET_KP_HAMILTONIAN(kp, dx, dy, dz)
        IMPLICIT NONE
        TYPE(KP_STENCIL), INTENT(OUT) :: kp
        REAL*8, INTENT(IN) :: dx, dy, dz

        kp%fijk = 0.0d0
        kp%fiplusjk = 0.0d0
        kp%fiminusjk = 0.0d0
        kp%fijplusk = 0.0d0
        kp%fijminusk = 0.0d0
        kp%fijkplus = 0.0d0
        kp%fijkminus = 0.0d0
        kp%fiplusjplusk = 0.0d0
        kp%fiplusjminusk = 0.0d0
        kp%fiminusjplusk = 0.0d0
        kp%fiminusjminusk = 0.0d0
        kp%fiplusjkplus = 0.0d0
        kp%fiplusjkminus = 0.0d0
        kp%fiminusjkplus = 0.0d0
        kp%fiminusjkminus = 0.0d0
        kp%fijpluskplus = 0.0d0
        kp%fijpluskminus = 0.0d0
        kp%fijminuskplus = 0.0d0
        kp%fijminuskminus = 0.0d0


        ! (i, j, k)
        kp%fijk(1,1) = 2*L/(dx**2) + M * (2/dy**2 + 2/dz**2)
        kp%fijk(4,4) = kp%fijk(1,1)
        kp%fijk(2,2) = 2*L/(dy**2) + M * (2/dx**2 + 2/dz**2)
        kp%fijk(5,5) = kp%fijk(2,2)
        kp%fijk(3,3) = 2*L/(dz**2) + M * (2/dx**2 + 2/dy**2)
        kp%fijk(6,6) = kp%fijk(3,3) 
        
        ! (i+1, j, k)
        kp%fiplusjk(1,1) = -L/(dx**2)
        kp%fiplusjk(4,4) = kp%fiplusjk(1,1)
        kp%fiplusjk(2,2) = -M/(dx**2)
        kp%fiplusjk(5,5) = kp%fiplusjk(2,2)
        kp%fiplusjk(3,3) = -M/(dx**2)
        kp%fiplusjk(6,6) = kp%fiplusjk(3,3) 

        ! (i-1, j, k)
        kp%fiminusjk(1,1) = -L/(dx**2)
        kp%fiminusjk(4,4) = kp%fiminusjk(1,1)
        kp%fiminusjk(2,2) = -M/(dx**2)
        kp%fiminusjk(5,5) = kp%fiminusjk(2,2)
        kp%fiminusjk(3,3) = -M/(dx**2)
        kp%fiminusjk(6,6) = kp%fiminusjk(3,3)

        ! (i, j+1, k)
        kp%fijplusk(1,1) = -M/(dy**2)
        kp%fijplusk(4,4) = kp%fijplusk(1,1)
        kp%fijplusk(2,2) = -L/(dy**2)
        kp%fijplusk(5,5) = kp%fijplusk(2,2)
        kp%fijplusk(3,3) = -M/(dy**2)
        kp%fijplusk(6,6) = kp%fijplusk(3,3)

        ! (i, j-1, k)
        kp%fijminusk(1,1) = -M/(dy**2)
        kp%fijminusk(4,4) = kp%fijminusk(1,1)
        kp%fijminusk(2,2) = -L/(dy**2)
        kp%fijminusk(5,5) = kp%fijminusk(2,2)
        kp%fijminusk(3,3) = -M/(dy**2)
        kp%fijminusk(6,6) = kp%fijminusk(3,3)

        ! (i, j, k+1)
        kp%fijkplus(1,1) = -M/(dz**2)
        kp%fijkplus(4,4) = kp%fijkplus(1,1)
        kp%fijkplus(2,2) = -M/(dz**2)
        kp%fijkplus(5,5) = kp%fijkplus(2,2)
        kp%fijkplus(3,3) = -L/(dz**2)
        kp%fijkplus(6,6) = kp%fijkplus(3,3) 

        ! (i, j, k-1)
        kp%fijkminus(1,1) = -M/(dz**2)
        kp%fijkminus(4,4) = kp%fijkminus(1,1)
        kp%fijkminus(2,2) = -M/(dz**2)
        kp%fijkminus(5,5) = kp%fijkminus(2,2)
        kp%fijkminus(3,3) = -L/(dz**2)
        kp%fijkminus(6,6) = kp%fijkminus(3,3)

        ! (i+1, j+1, k)
        kp%fiplusjplusk(1,2) = -N/(2*dx*dy)
        kp%fiplusjplusk(2,1) = kp%fiplusjplusk(1,2)
        kp%fiplusjplusk(4,5) = kp%fiplusjplusk(1,2)
        kp%fiplusjplusk(5,4) = kp%fiplusjplusk(1,2)
        ! (i+1, j-1, k)
        kp%fiplusjminusk(1,2) = N/(2*dx*dy)
        kp%fiplusjminusk(2,1) = kp%fiplusjminusk(1,2)
        kp%fiplusjminusk(4,5) = kp%fiplusjminusk(1,2)
        kp%fiplusjminusk(5,4) = kp%fiplusjminusk(1,2)
        ! (i-1, j+1, k)
        kp%fiminusjplusk(1,2) = N/(2*dx*dy)
        kp%fiminusjplusk(2,1) = kp%fiminusjplusk(1,2)
        kp%fiminusjplusk(4,5) = kp%fiminusjplusk(1,2)
        kp%fiminusjplusk(5,4) = kp%fiminusjplusk(1,2)
        ! (i-1, j-1, k)
        kp%fiminusjminusk(1,2) = -N/(2*dx*dy)
        kp%fiminusjminusk(2,1) = kp%fiminusjminusk(1,2)
        kp%fiminusjminusk(4,5) = kp%fiminusjminusk(1,2)
        kp%fiminusjminusk(5,4) = kp%fiminusjminusk(1,2)

        ! (i+1, j, k+1)
        kp%fiplusjkplus(1,3) = -N/(2*dx*dz)
        kp%fiplusjkplus(3,1) = kp%fiplusjkplus(1,3)
        kp%fiplusjkplus(4,6) = kp%fiplusjkplus(1,3)
        kp%fiplusjkplus(6,4) = kp%fiplusjkplus(1,3)
        ! (i+1, j, k-1)
        kp%fiplusjkminus(1,3) = N/(2*dx*dz)
        kp%fiplusjkminus(3,1) = kp%fiplusjkminus(1,3)
        kp%fiplusjkminus(4,6) = kp%fiplusjkminus(1,3)
        kp%fiplusjkminus(6,4) = kp%fiplusjkminus(1,3)
        ! (i-1, j, k+1)
        kp%fiminusjkplus(1,3) = N/(2*dx*dz)
        kp%fiminusjkplus(3,1) = kp%fiminusjkplus(1,3)
        kp%fiminusjkplus(4,6) = kp%fiminusjkplus(1,3)
        kp%fiminusjkplus(6,4) = kp%fiminusjkplus(1,3)
        ! (i-1, j, k-1)
        kp%fiminusjkminus(1,3) = -N/(2*dx*dz)
        kp%fiminusjkminus(3,1) = kp%fiminusjkminus(1,3)
        kp%fiminusjkminus(4,6) = kp%fiminusjkminus(1,3)
        kp%fiminusjkminus(6,4) = kp%fiminusjkminus(1,3)

        ! (i, j+1, k+1)
        kp%fijpluskplus(2,3) = -N/(2*dx*dz)
        kp%fijpluskplus(3,2) = kp%fijpluskplus(2,3)
        kp%fijpluskplus(5,6) = kp%fijpluskplus(2,3)
        kp%fijpluskplus(6,5) = kp%fijpluskplus(2,3)
        ! (i, j+1, k-1)
        kp%fijpluskminus(2,3) = N/(2*dx*dz)
        kp%fijpluskminus(3,2) = kp%fijpluskminus(2,3)
        kp%fijpluskminus(5,6) = kp%fijpluskminus(2,3)
        kp%fijpluskminus(6,5) = kp%fijpluskminus(2,3)
        ! (i, j-1, k+1)
        kp%fijminuskplus(2,3) = N/(2*dx*dz)
        kp%fijminuskplus(3,2) = kp%fijminuskplus(2,3)
        kp%fijminuskplus(5,6) = kp%fijminuskplus(2,3)
        kp%fijminuskplus(6,5) = kp%fijminuskplus(2,3)
        ! (i, j-1, k-1)
        kp%fijminuskminus(2,3) = -N/(2*dx*dz)
        kp%fijminuskminus(3,2) = kp%fijminuskminus(2,3)
        kp%fijminuskminus(5,6) = kp%fijminuskminus(2,3)
        kp%fijminuskminus(6,5) = kp%fijminuskminus(2,3)

    END SUBROUTINE GET_KP_HAMILTONIAN


    SUBROUTINE IMAGINARY_TIME_KP(potential, Nx, Ny, Nz, dx, dy, dz, dt, MAX_TIME, kp, init_psi, final_psi, final_energy, tol)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: potential(:,:,:)
        REAL*8, INTENT(IN) :: init_psi(:,:,:,:)
        REAL*8, INTENT(OUT) :: final_psi(norb, Nx,Ny,Nz)
        REAL*8, INTENT(OUT) :: final_energy
        REAL*8, INTENT(IN) :: dx, dy, dz, dt, tol
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz, MAX_TIME
        TYPE(KP_STENCIL), INTENT(IN) :: kp
        REAL*8, ALLOCATABLE :: psi(:,:,:,:)
        REAL*8, ALLOCATABLE :: psi_new(:,:,:,:)
        REAL*8, ALLOCATABLE :: ham(:,:,:,:)
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

        ALLOCATE(psi(norb, Nx, Ny, Nz))
        ALLOCATE(psi_new(norb, Nx, Ny, Nz))
        ALLOCATE(ham(norb, Nx, Ny, Nz))
        ham(:,:,:,:) = 0.0d0
        psi = init_psi
       
        energy_old = 1.d99
        ! tol = 1.d-6
        DO iter =1, MAX_TIME
            ham(:,:,:,:) = 0.0d0
            !!$omp parallel do collapse(3) schedule(static)
            DO k=2, Nz-1
                DO  j=2, Ny-1
                    DO i=2, Nx-1
                        ham(:,i,j,k) = &
                        MATMUL(kp%fijk, psi(:,i,j,k)) +&
                        MATMUL(kp%fiplusjk, psi(:,i+1,j,k)) +&
                        MATMUL(kp%fiminusjk, psi(:,i-1,j,k)) +&
                        MATMUL(kp%fijplusk, psi(:,i,j+1,k)) +&
                        MATMUL(kp%fijminusk, psi(:,i,j-1,k)) +&
                        MATMUL(kp%fijkplus, psi(:,i,j,k+1)) +&
                        MATMUL(kp%fijkminus, psi(:,i,j,k-1)) +&
                        MATMUL(kp%fiplusjplusk, psi(:,i+1,j+1,k)) +&
                        MATMUL(kp%fiplusjminusk, psi(:,i+1,j-1,k)) +&
                        MATMUL(kp%fiminusjplusk, psi(:,i-1,j+1,k)) +&
                        MATMUL(kp%fiminusjminusk, psi(:,i-1,j-1,k)) +&
                        MATMUL(kp%fiplusjkplus, psi(:,i+1,j,k+1)) +&
                        MATMUL(kp%fiplusjkminus, psi(:,i+1,j,k-1)) +&
                        MATMUL(kp%fiminusjkplus, psi(:,i-1,j,k+1)) +&
                        MATMUL(kp%fiminusjkminus, psi(:,i-1,j,k-1)) +&
                        MATMUL(kp%fijpluskplus, psi(:,i,j+1,k+1)) +&
                        MATMUL(kp%fijpluskminus, psi(:,i,j+1,k-1)) +&
                        MATMUL(kp%fijminuskplus, psi(:,i,j-1,k+1)) +&
                        MATMUL(kp%fijminuskminus, psi(:,i,j-1,k-1)) +&
                        potential(i,j,k) * psi(:,i,j,k)

                        !-1/(2.0d0*m1)*(( psi(i+1,j,k)+psi(i-1,j,k) + psi(i,j+1,k)+psi(i,j-1,k)-&
                        !4.d0*psi(i,j,k))/(dx**2) + (psi(i,j,k+1)+psi(i,j,k-1) - 2.d0*psi(i,j,k) )/dz**2) &
                        !+ potential(i,j,k) * psi(i,j,k)
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
        
END MODULE SCHRODINGER