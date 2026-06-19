!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Schrodinger Poisson method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Poisson_Solver
    USE CONSTANTS
    USE DIELECTRIC
    IMPLICIT NONE

CONTAINS

    ! =====================================================================
    ! MAIN MULTILEVEL WRAPPER
    ! =====================================================================
    SUBROUTINE Poisson_Multilevel(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, n_levels, tol, MAX_ITER)
        INTEGER, INTENT(IN)           :: Nx, Ny, Nz
        REAL*8,  INTENT(INOUT)        :: potential(Nx, Ny, Nz)
        REAL*8,  INTENT(IN)           :: density(Nx, Ny, Nz)
        REAL*8,  INTENT(IN)           :: epsilon(Nx, Ny, Nz)
        REAL*8,  INTENT(IN)           :: alfa, dx
        INTEGER, INTENT(IN)           :: n_levels ! Power of 2 to reduce by
        REAL*8, INTENT(IN)            :: tol
        INTEGER, INTENT(IN)           :: MAX_ITER           

        REAL*8, ALLOCATABLE :: pot_c(:,:,:), dens_c(:,:,:), eps_c(:,:,:)
        REAL*8, ALLOCATABLE :: pot_f(:,:,:), dens_f(:,:,:), eps_f(:,:,:)
        INTEGER :: L, Ncx, Ncy, Ncz, Nfx, Nfy, Nfz, stride, i, j, k
        REAL*8  :: dx_curr

        ! 1. START AT THE COARSEST LEVEL
        stride = 2**n_levels
        Ncx = (Nx - 1) / stride + 1
        Ncy = (Ny - 1) / stride + 1
        Ncz = (Nz - 1) / stride + 1
        dx_curr = dx * stride

        ALLOCATE(pot_c(Ncx, Ncy, Ncz), dens_c(Ncx, Ncy, Ncz), eps_c(Ncx, Ncy, Ncz))
        
        ! Injection: Sample fine grid to coarse grid
        DO k = 1, Ncz
        DO j = 1, Ncy
        DO i = 1, Ncx
            eps_c(i,j,k)  = epsilon((i-1)*stride+1, (j-1)*stride+1, (k-1)*stride+1)
            dens_c(i,j,k) = density((i-1)*stride+1, (j-1)*stride+1, (k-1)*stride+1)
            pot_c(i,j,k)  = 0.0d0
        END DO
        END DO
        END DO

        PRINT *, "Solving Level 0 (Coarsest): ", Ncx, "x", Ncy, "x", Ncz
        CALL Poisson_Kernel(pot_c, dens_c, eps_c, alfa, Ncx, Ncy, Ncz, dx_curr, tol, MAX_ITER) ! initial

        ! 2. PROGRESSIVE REFINEMENT LOOP
        DO L = n_levels-1, 0, -1
            stride = 2**L
            Nfx = (Nx - 1) / stride + 1
            Nfy = (Ny - 1) / stride + 1
            Nfz = (Nz - 1) / stride + 1
            dx_curr = dx * stride

            ALLOCATE(pot_f(Nfx, Nfy, Nfz), dens_f(Nfx, Nfy, Nfz), eps_f(Nfx, Nfy, Nfz))
            
            ! Prolongate (Interpolate) previous solution to this level
            CALL Prolongate_3D(pot_c, pot_f, Ncx, Ncy, Ncz, Nfx, Nfy, Nfz)
            
            ! Sample physics for this level
            DO k = 1, Nfz
            DO j = 1, Nfy
            DO i = 1, Nfx
                eps_f(i,j,k)  = epsilon((i-1)*stride+1, (j-1)*stride+1, (k-1)*stride+1)
                dens_f(i,j,k) = density((i-1)*stride+1, (j-1)*stride+1, (k-1)*stride+1)
            END DO
            END DO
            END DO

            PRINT *, "Solving Level: ", n_levels - L, " Size: ", Nfx
            CALL Poisson_Kernel(pot_f, dens_f, eps_f, alfa, Nfx, Nfy, Nfz, dx_curr, tol, MAX_ITER) ! iteracyjnie

            ! Prepare for next level
            DEALLOCATE(pot_c, dens_c, eps_c)
            IF (L > 0) THEN
                ALLOCATE(pot_c(Nfx, Nfy, Nfz), dens_c(Nfx, Nfy, Nfz), eps_c(Nfx, Nfy, Nfz))
                pot_c = pot_f
                dens_c = dens_f
                eps_c = eps_f
                Ncx = Nfx; Ncy = Nfy; Ncz = Nfz
            ELSE
                potential = pot_f ! Final fine grid
            END IF
            DEALLOCATE(pot_f, dens_f, eps_f)
        END DO

    END SUBROUTINE Poisson_Multilevel

    ! =====================================================================
    ! GAUSS-SEIDEL / SOR KERNEL (Inhomogeneous)
    ! =====================================================================
    SUBROUTINE Poisson_Kernel(phi, rho, eps, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER)
        INTEGER, INTENT(IN)    :: Nx, Ny, Nz, MAX_ITER
        REAL*8,  INTENT(IN)    :: dx, alfa, tol, rho(Nx,Ny,Nz), eps(Nx,Ny,Nz)
        REAL*8,  INTENT(INOUT) :: phi(Nx,Ny,Nz)
        INTEGER :: i, j, k, iter
        REAL*8  :: e_ip, e_im, e_jp, e_jm, e_kp, e_km, sum_e, source, res, max_err

        DO iter = 1, MAX_ITER
            max_err = 0.0d0
            DO k = 2, Nz-1
            DO j = 2, Ny-1
            DO i = 2, Nx-1
                ! Averaged epsilons at cell faces
                e_ip = 0.5d0 * (eps(i+1,j,k) + eps(i,j,k))
                e_im = 0.5d0 * (eps(i-1,j,k) + eps(i,j,k))
                e_jp = 0.5d0 * (eps(i,j+1,k) + eps(i,j,k))
                e_jm = 0.5d0 * (eps(i,j-1,k) + eps(i,j,k))
                e_kp = 0.5d0 * (eps(i,j,k+1) + eps(i,j,k))
                e_km = 0.5d0 * (eps(i,j,k-1) + eps(i,j,k))
                
                sum_e = e_ip + e_im + e_jp + e_jm + e_kp + e_km
                source = (rho(i,j,k) / epsilon0) * dx**2
                
                res = (e_ip*phi(i+1,j,k) + e_im*phi(i-1,j,k) + &
                       e_jp*phi(i,j+1,k) + e_jm*phi(i,j-1,k) + &
                       e_kp*phi(i,j,k+1) + e_km*phi(i,j,k-1) + source) / sum_e
                
                max_err = MAX(max_err, ABS(phi(i,j,k) - res))
                phi(i,j,k) = (1.0d0 - alfa) * phi(i,j,k) + alfa * res
            END DO
            END DO
            END DO
            IF (max_err < tol) EXIT
        END DO
    END SUBROUTINE Poisson_Kernel

    ! =====================================================================
    ! PROLONGATION (Linear Interpolation)
    ! =====================================================================
    SUBROUTINE Prolongate_3D(P_c, P_f, Ncx, Ncy, Ncz, Nfx, Nfy, Nfz)
        INTEGER, INTENT(IN) :: Ncx, Ncy, Ncz, Nfx, Nfy, Nfz
        REAL*8,  INTENT(IN) :: P_c(Ncx, Ncy, Ncz)
        REAL*8,  INTENT(OUT):: P_f(Nfx, Nfy, Nfz)
        INTEGER :: i, j, k, fi, fj, fk

        P_f = 0.0d0
        DO k = 1, Ncz - 1
        DO j = 1, Ncy - 1
        DO i = 1, Ncx - 1
            fi = 2*i-1; fj = 2*j-1; fk = 2*k-1
            
            ! Direct mapping
            P_f(fi,fj,fk) = P_c(i,j,k)
            
            ! Linear midpoints
            P_f(fi+1, fj, fk) = 0.5d0 * (P_c(i,j,k) + P_c(i+1,j,k))
            P_f(fi, fj+1, fk) = 0.5d0 * (P_c(i,j,k) + P_c(i,j+1,k))
            P_f(fi, fj, fk+1) = 0.5d0 * (P_c(i,j,k) + P_c(i,j,k+1))
            
            ! Face/Cell centers (Averages)
            P_f(fi+1,fj+1,fk)   = 0.25d0 * (P_c(i,j,k)+P_c(i+1,j,k)+P_c(i,j+1,k)+P_c(i+1,j+1,k))
            P_f(fi+1,fj+1,fk+1) = 0.125d0 * (P_c(i,j,k) + P_c(i+1,j,k) + P_c(i,j+1,k) + &
                                             P_c(i,j,k+1) + P_c(i+1,j+1,k) + P_c(i+1,j,k+1) + &
                                             P_c(i,j+1,k+1) + P_c(i+1,j+1,k+1))
        END DO
        END DO
        END DO
        ! Boundary faces (simplified)
        P_f(Nfx,:,:) = 0.0d0; P_f(:,Nfy,:) = 0.0d0; P_f(:,:,Nfz) = 0.0d0 
    END SUBROUTINE Prolongate_3D


    !!!!!!!!!!!! normal Gauss-Seidl method for a spatially changing dielectric constant !!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Poisson_epsilon(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, dz, tol, MAX_ITER, charge_trapped)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)           :: Nx, Ny, Nz, MAX_ITER
        REAL*8,  INTENT(INOUT)        :: potential(:,:,:)
        REAL*8,  INTENT(IN)           :: density(:,:,:)
        REAL*8,  INTENT(IN)           :: epsilon(:,:,:)
        REAL*8,  INTENT(IN)           :: charge_trapped(:,:,:)
        REAL*8,  INTENT(IN)           :: alfa, dx, tol, dz
        
        REAL*8, PARAMETER :: pi       = 3.141592653589793d0
        REAL*8, PARAMETER :: epsilon0 = 1.0d0/(4.0d0 * pi)
        
        INTEGER           :: i, j, k, iter
        REAL*8            :: e_ip, e_im, e_jp, e_jm, e_kp, e_km
        REAL*8            :: sum_e, res, old_val, max_err, source, val

        REAL*8, ALLOCATABLE :: sigma_2d(:,:)
        
        allocate(sigma_2d(Nx,Ny))

        !determine sigma_2d
        do i = 1, Nx
            do j = 1, Ny
            val=0.0
            do k = 1, Nz
                val=val+(density(i,j,k)+charge_trapped(i, j, k))*dz
            enddo
            sigma_2d(i,j)=val
            ENDDO
        enddo


        do iter = 1, MAX_ITER
            max_err = 0.0d0

            ! --- 1. Update the Mapped Boundary at z=0 (k=1) ---
            k = 1
            do j = 2, Ny - 1
                do i = 2, Nx - 1
                    old_val = potential(i,j,k)
                    
                    res=potential(i,j,k+1) +  dx * sigma_2d(i,j) /epsilon(i,j,k)/epsilon0 ! 

                    potential(i,j,k) = (1.0d0 - alfa) * old_val + alfa * res
                    max_err = max(max_err, abs(potential(i,j,k) - old_val))
                end do
            end do
            
            ! ---------2. Neumann BC at x and y boundaries (zero normal derivative) ---------
            i=1
            do k = 2, Nz - 1
            do j = 1, Ny 
                potential(i,j,k) = potential(i+1,j,k)
            enddo
            enddo

            i=Nx
            do k = 2, Nz - 1
            do j = 1, Ny 
                potential(i,j,k) = potential(i-1,j,k)
            enddo
            enddo

            j=1
            do k = 2, Nz - 1
            do i = 1, Nx
                potential(i,j,k) = potential(i,j+1,k)
            enddo
            enddo

            j=Ny
            do k = 2, Nz - 1
            do i = 1, Nx
                potential(i,j,k) = potential(i,j-1,k)
            enddo
            enddo


            do k = 2, Nz - 1
                do j = 2, Ny - 1
                    do i = 2, Nx - 1
                        
                        old_val = potential(i,j,k)
                        
                        ! 1. Calculate epsilon at mid-points (half-steps)
                        e_ip = (epsilon(i+1,j,k) + epsilon(i,j,k)) / 2.0d0
                        e_im = (epsilon(i-1,j,k) + epsilon(i,j,k)) / 2.0d0
                        e_jp = (epsilon(i,j+1,k) + epsilon(i,j,k)) / 2.0d0
                        e_jm = (epsilon(i,j-1,k) + epsilon(i,j,k)) / 2.0d0
                        e_kp = (epsilon(i,j,k+1) + epsilon(i,j,k)) / 2.0d0
                        e_km = (epsilon(i,j,k-1) + epsilon(i,j,k)) / 2.0d0
                        
                        sum_e = e_ip + e_im + e_jp + e_jm + e_kp + e_km
                        
                        ! 2. The source term (Right Hand Side)
                        source = -(density(i,j,k) + charge_trapped(i,j,k)) / epsilon0 * dx * dx
                        
                        ! 3. Calculate Gauss-Seidel result
                        ! This is now a weighted average of neighbors
                        res = (e_ip*potential(i+1,j,k) + e_im*potential(i-1,j,k) + &
                            e_jp*potential(i,j+1,k) + e_jm*potential(i,j-1,k) + &
                            e_kp*potential(i,j,k+1) + e_km*potential(i,j,k-1) + &
                            source) / sum_e
                        
                        ! 4. SOR Mixing
                        potential(i,j,k) = (1.0d0 - alfa) * old_val + alfa * res
                        
                        max_err = max(max_err, abs(potential(i,j,k) - old_val))
                     end do
                end do
            end do

            !BC at z=0

            print *, iter, max_err
            if (max_err < TOL) exit
        end do
        print *, "Converged in ", iter, " iterations. Final error: ", max_err
        DEALLOCATE(sigma_2d)

    END SUBROUTINE Poisson_epsilon



     !!!!!!!!!!!! normal Gauss-Seidl method for a spatially changing dielectric constant !!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Poisson_epsilon_no_charge(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, dz, tol, MAX_ITER)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)           :: Nx, Ny, Nz, MAX_ITER
        REAL*8,  INTENT(INOUT)        :: potential(:,:,:)
        REAL*8,  INTENT(IN)           :: density(:,:,:)
        REAL*8,  INTENT(IN)           :: epsilon(:,:,:)
        REAL*8,  INTENT(IN)           :: alfa, dx, tol, dz
        
        INTEGER           :: i, j, k, iter
        REAL*8            :: e_ip, e_im, e_jp, e_jm, e_kp, e_km
        REAL*8            :: sum_e, res, old_val, max_err, source, val
        print*, maxval(density)
        print*, 1/(2*pi*(5.0*fnm2au)**2)**(1.5)

        do iter = 1, MAX_ITER
            max_err = 0.0d0

            ! --- 1. Update the Mapped Boundary at z=0 (k=1) ---
            k = 1
            do j = 2, Ny - 1
                do i = 2, Nx - 1
                    old_val = potential(i,j,k)
                    
                    res=potential(i,j,k+1) 

                    potential(i,j,k) = (1.0d0 - alfa) * old_val + alfa * res
                    max_err = max(max_err, abs(potential(i,j,k) - old_val))
                end do
            end do
            
            ! ---------2. Neumann BC at x and y boundaries (zero normal derivative) ---------
            ! update: Dirichlet boundary conditions
            i=1
            do k = 2, Nz - 1
            do j = 1, Ny 
                potential(i,j,k) = 0.0!potential(i+1,j,k)
            enddo
            enddo

            i=Nx
            do k = 2, Nz - 1
            do j = 1, Ny 
                potential(i,j,k) = 0.0 !potential(i-1,j,k)
            enddo
            enddo

            j=1
            do k = 2, Nz - 1
            do i = 1, Nx
                potential(i,j,k) = 0.0! potential(i,j+1,k)
            enddo
            enddo

            j=Ny
            do k = 2, Nz - 1
            do i = 1, Nx
                potential(i,j,k) = 0.0 !potential(i,j-1,k)
            enddo
            enddo


            do k = 2, Nz - 1
                do j = 2, Ny - 1
                    do i = 2, Nx - 1
                        
                        old_val = potential(i,j,k)
                        
                        ! 1. Calculate epsilon at mid-points (half-steps)
                        e_ip = (epsilon(i+1,j,k) + epsilon(i,j,k)) / 2.0d0
                        e_im = (epsilon(i-1,j,k) + epsilon(i,j,k)) / 2.0d0
                        e_jp = (epsilon(i,j+1,k) + epsilon(i,j,k)) / 2.0d0
                        e_jm = (epsilon(i,j-1,k) + epsilon(i,j,k)) / 2.0d0
                        e_kp = (epsilon(i,j,k+1) + epsilon(i,j,k)) / 2.0d0
                        e_km = (epsilon(i,j,k-1) + epsilon(i,j,k)) / 2.0d0
                        
                        sum_e = e_ip + e_im + e_jp + e_jm + e_kp + e_km
                        
                        ! 2. The source term (Right Hand Side)
                        source = -(density(i,j,k)) / epsilon0 * dx * dx
                        
                        ! 3. Calculate Gauss-Seidel result
                        ! This is now a weighted average of neighbors
                        res = (e_ip*potential(i+1,j,k) + e_im*potential(i-1,j,k) + &
                            e_jp*potential(i,j+1,k) + e_jm*potential(i,j-1,k) + &
                            e_kp*potential(i,j,k+1) + e_km*potential(i,j,k-1) + &
                            source) / sum_e
                        
                        ! 4. SOR Mixing
                        potential(i,j,k) = (1.0d0 - alfa) * old_val + alfa * res
                        
                        max_err = max(max_err, abs(potential(i,j,k) - old_val))
                     end do
                end do
            end do

            !BC at z=0

            print *, iter, max_err
            if (max_err < TOL) exit
        end do
        print *, "Converged in ", iter, " iterations. Final error: ", max_err

    END SUBROUTINE Poisson_epsilon_no_charge


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! normal Gauss-Seidl method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Poisson(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER)
        IMPLICIT NONE
        
        ! Arguments
        INTEGER, INTENT(IN)           :: Nx, Ny, Nz, MAX_ITER
        REAL*8,  INTENT(INOUT)        :: potential(:,:,:)
        REAL*8,  INTENT(IN)           :: density(:,:,:)
        REAL*8,  INTENT(IN)           :: alfa, epsilon
        REAL*8,  INTENT(IN)           :: dx, tol
        
        ! Internal Solver Parameters
        INTEGER           :: i, j, k, iter
        REAL*8            :: source, res, old_val, max_err
        print *, "Entered Poisson"
        print *, shape(potential)
        print *, shape(density)

        ! Iterative Loop
        do iter = 1, MAX_ITER
            max_err = 0.0d0
            

            ! Perform Gauss-Seidel update on interior nodes
            ! Boundaries (1 and N) are treated as Dirichlet (fixed)
            do k = 2, Nz - 1
                do j = 2, Ny - 1
                    do i = 2, Nx - 1
                        
                        source = density(i,j,k)/epsilon0/epsilon
                        old_val = potential(i,j,k)
                        
                        ! 2. Calculate the updated value based on neighbors
                        ! Standard 7-point finite difference stencil
                        res = (potential(i+1,j,k) + potential(i-1,j,k) + &
                            potential(i,j+1,k) + potential(i,j-1,k) + &
                            potential(i,j,k+1) + potential(i,j,k-1) - &
                            source * dx * dx) / 6.0d0
                        
                        ! 3. Apply mixing (Successive Over-Relaxation)
                        potential(i,j,k) = (1.0d0 - alfa) * old_val + alfa * res
                        ! Track maximum change for convergence
                        max_err = max(max_err, abs(potential(i,j,k) - old_val))
                        
                    end do
                end do
            end do

            print *, "Poisson at iteration: iter, max_err ", iter, max_err
            ! Convergence Check
            if (max_err < tol) then
                print *, "Poisson converged at iteration: ", iter, max_err
                exit
            end if
            
            if (iter == MAX_ITER) then
                print *, "Warning: Poisson solver reached MAX_ITER without convergence."
            end if
        end do

    END SUBROUTINE Poisson


    SUBROUTINE POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0,  nz, dz, charge_trapped, electric_field, pot_hartree)
        IMPLICIT NONE

        REAL*8, intent(in) :: eps_0
        REAL*8, intent(in) :: L_trapped
        REAL*8, intent(in) :: n0_trapped
        REAL*8, intent(in) :: dz
        INTEGER, intent(in) :: nz
        REAL*8, INTENT(OUT) :: charge_trapped(nz)
        REAL*8, INTENT(OUT) :: electric_field(nz)
        REAL*8, INTENT(OUT) :: pot_hartree(nz)

        !zmienne pomocnicze
        INTEGER :: i, iz
        REAL*8 :: charge_bc
        REAL*8 :: z

        !PARDISO       
        INTEGER, ALLOCATABLE :: nmat(:)
        INTEGER :: pt_prd(64)= 0.
        INTEGER :: maxfct_prd, mnum_prd, mtype_prd, phase_prd
        INTEGER :: n_prd, nrhs_prd, error_prd
        INTEGER, ALLOCATABLE :: ia_prd(:)
        INTEGER, ALLOCATABLE :: ja_prd(:)
        INTEGER, ALLOCATABLE ::  perm_prd(:)
        INTEGER :: iparm_prd(64), msglvl_prd
        REAL*8, ALLOCATABLE :: a_prd(:)
        REAL*8, ALLOCATABLE :: b_prd(:,:)
        REAL*8, ALLOCATABLE :: x_prd(:,:)
        INTEGER :: maxnonzeroprd, nelem

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!! tablica charge_trapped !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !trapped electrons in STO (note that there are eletrons !!! )
        do iz=1,nz
            z=(iz-1)*dz
            charge_trapped(iz)=(n0_trapped/L_trapped)*dexp(-z/L_trapped)
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!! rozwiazanie rownania Poissona dla pierwszej iteracji z charge_trapped !!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nrhs_prd=1 ! number of rhs
        n_prd=nz ! size of the matrix
        maxnonzeroprd=(n_prd-2)*3+2+1 ! max non-zero

        ALLOCATE( perm_prd(n_prd) )
        ALLOCATE( ia_prd(n_prd+1) )
        ALLOCATE( ja_prd(maxnonzeroprd) )
        ALLOCATE( a_prd(maxnonzeroprd) )
        ALLOCATE( b_prd(n_prd,nrhs_prd) )
        ALLOCATE( x_prd(n_prd,nrhs_prd) )

        pt_prd(:)= 0 ! default solver handle
        maxfct_prd= 1 ! max number of factorzation matrix
        mnum_prd= 1 ! which factorization
        mtype_prd= 11   ! real nonsymmetric  
        phase_prd= 13 ! analysis, numerical factorization and solve
        msglvl_prd= 0 ! no verobse
        iparm_prd(:)= 0 ! default params

        perm_prd(:)= 0 ! default permutation
        ia_prd(:)= 0 
        ja_prd(:)= 0
        a_prd(:)= 0. 
        b_prd(:,:)= 0. ! right hand side
        x_prd(:,:)= 0. ! solution vector


        charge_bc=0.0
        do iz=1,nz
            charge_bc=charge_bc+charge_trapped(iz)*dz
        enddo

        nelem=1
        DO i=1,nz

            IF (i.eq.1) THEN
                ia_prd(i)= nelem

                ja_prd( nelem )=  i
                a_prd( nelem )= -1.0
                nelem= nelem + 1

                ja_prd( nelem )=  i+1
                a_prd( nelem )= 1.0
                nelem= nelem + 1
                
            ELSE IF (i.eq.(nz)) THEN
                ia_prd(i)= nelem
                
                ja_prd( nelem )=  i
                a_prd( nelem )= 1.0
                nelem= nelem + 1

            ELSE
                ia_prd(i)= nelem

                ja_prd(nelem)=i-1
                a_prd(nelem)= 1.0
                nelem= nelem + 1

                ja_prd( nelem )=  i
                a_prd( nelem )= -2.0
                nelem= nelem + 1

                ja_prd( nelem )=  i+1
                a_prd( nelem )= 1.0
                nelem= nelem + 1
            ENDIF

        enddo

        ia_prd( n_prd+1 )= nelem
        IF (nelem-1 /= maxnonzeroprd) STOP "PARDISO:  nelem /= maxnonzeroprd"

        !wektor wyrazow wolnych
        DO i=2,nz-1
        b_prd(i, 1)= -(-charge_trapped(i))*dz*dz/epsilon0/eps_0
        ENDDO
        !warunki brzegowe
        b_prd( 1, 1 )= -charge_bc*dz/epsilon0/eps_0
        b_prd( nz, 1 )= 0.0

        CALL pardiso (pt_prd, maxfct_prd, mnum_prd, mtype_prd, phase_prd,      &
                &        n_prd, a_prd, ia_prd, ja_prd, perm_prd, nrhs_prd,     &
                &        iparm_prd, msglvl_prd, b_prd, x_prd, error_prd)
        IF (error_prd /= 0)  THEN
            PRINT*, "pardiso: error_prd =", error_prd
            STOP
        END IF

        do iz=1,nz
            pot_hartree(iz)=x_prd(iz,1)
        enddo

        do iz=1,nz-1
            electric_field(iz)=(-pot_hartree(iz+1)+pot_hartree(iz))/dz
        enddo
        electric_field(nz)=electric_field(nz-1)

        !!!!!!!!!!!!! zapis do pliku !!!!!!!!!!!!!!!!!!
        OPEN(1, FILE="data/charge_trapped.dat")
        do iz=1,nz
            z=(iz-1)*dz
            write(1, '(200e20.12)') z/fnm2au, charge_trapped(iz)/fne2au   
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/potential_trapped.dat")
        do iz=1,nz
        z=(iz-1)*dz
        write(1, '(200e20.12)') z/fnm2au, -pot_hartree(iz)/feV2au 
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/electric_field_trapped.dat")
        do iz=1,nz
        z=(iz-1)*dz
        write(1, '(200e20.12)') z/fnm2au, electric_field(iz)*fnm2au/feV2au 
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/epsilon_trapped.dat")
        do iz=1,nz
            z=(iz-1)*dz
            write(1, '(200e20.12)') z/fnm2au, permitivity(eps_0,electric_field(iz))
        enddo
        CLOSE(1)
        DEALLOCATE(perm_prd)
        DEALLOCATE(ia_prd)
        DEALLOCATE(ja_prd)
        DEALLOCATE(a_prd)
        DEALLOCATE(b_prd)
        DEALLOCATE(x_prd)

        
    END SUBROUTINE POISSON_ZDIRECTION_INIT


    SUBROUTINE POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)
        ! solving in z direction but with changing epsilon
        IMPLICIT NONE
        REAL*8, INTENT (IN) :: electric_field(nz)
        REAL*8, intent(in) :: eps_0
        REAL*8, intent(in) :: charge_trapped(nz)
        REAL*8, intent(in) :: dz
        INTEGER, intent(in) :: nz
        REAL*8, INTENT (OUT) :: electric_field_new(nz)

        INTEGER :: i, iz
        REAL*8 :: charge_bc, z
        !PARDISO       
        INTEGER :: pt_prd(64)= 0.
        INTEGER :: maxfct_prd, mnum_prd, mtype_prd, phase_prd
        INTEGER :: n_prd, nrhs_prd, error_prd
        INTEGER, ALLOCATABLE :: ia_prd(:)
        INTEGER, ALLOCATABLE :: ja_prd(:)
        INTEGER, ALLOCATABLE ::  perm_prd(:)
        INTEGER :: iparm_prd(64), msglvl_prd
        REAL*8, ALLOCATABLE :: a_prd(:)
        REAL*8, ALLOCATABLE :: b_prd(:,:)
        REAL*8, ALLOCATABLE :: x_prd(:,:)
        INTEGER :: maxnonzeroprd, nelem

        REAL*8, ALLOCATABLE :: pot_hartree(:)

        allocate(pot_hartree(nz))
        nrhs_prd=1
        n_prd=nz
        maxnonzeroprd=(n_prd-2)*3+2+1

        ALLOCATE( perm_prd(n_prd) )
        ALLOCATE( ia_prd(n_prd+1) )
        ALLOCATE( ja_prd(maxnonzeroprd) )
        ALLOCATE( a_prd(maxnonzeroprd) )
        ALLOCATE( b_prd(n_prd,nrhs_prd) )
        ALLOCATE( x_prd(n_prd,nrhs_prd) )
        
        charge_bc=0.0
        do iz=1,nz
            charge_bc=charge_bc+(charge_trapped(iz))*dz
        enddo
	
        perm_prd(:)= 0
        ia_prd(:)= 0
        ja_prd(:)= 0
        a_prd(:)= 0.
        b_prd(:,:)= 0.
        x_prd(:,:)= 0.
        pot_hartree(:)= 0.

        !uzupelnianie macierzy
            nelem=1
            DO i=1,nz
            IF (i.eq.1) THEN
                ia_prd(i)= nelem

                ja_prd( nelem )=  i
                a_prd( nelem )= -1.0
                nelem= nelem + 1

                ja_prd( nelem )=  i+1
                a_prd( nelem )= 1.0
                nelem= nelem + 1
            
            ELSE IF (i.eq.(nz)) THEN
                ia_prd(i)= nelem
                
                ja_prd( nelem )=  i
                a_prd( nelem )= 1.0
                nelem= nelem + 1

            ELSE
                ia_prd(i)= nelem

                ja_prd(nelem)=i-1
                a_prd(nelem)= 0.5*(permitivity(eps_0,electric_field(i))+permitivity(eps_0,electric_field(i-1))) 
                nelem= nelem + 1

                ja_prd( nelem )=  i
                a_prd( nelem )= -0.5*(permitivity(eps_0,electric_field(i+1))+ &
                            & 2.0*permitivity(eps_0,electric_field(i))+permitivity(eps_0,electric_field(i-1))) 
                nelem= nelem + 1

                ja_prd( nelem )=i+1
                a_prd( nelem )= 0.5*(permitivity(eps_0,electric_field(i))+permitivity(eps_0,electric_field(i+1))) 
                nelem= nelem + 1
            ENDIF
            enddo


        ia_prd( n_prd+1 )= nelem
        IF (nelem-1 /= maxnonzeroprd) STOP "PARDISO:  nelem /= maxnonzeroprd"

        !wektor wyrazow wolnych
        DO i=2,nz-1
               b_prd(i, 1)= -(-charge_trapped(i))*dz*dz/epsilon0
        ENDDO
        
        !warunki brzegowe
        b_prd( 1, 1 )= -charge_bc*dz/permitivity(eps_0,electric_field(1))/epsilon0
        b_prd( nz, 1 )= 0.0


        CALL pardiso (pt_prd, maxfct_prd, mnum_prd, mtype_prd, phase_prd,      &
                &        n_prd, a_prd, ia_prd, ja_prd, perm_prd, nrhs_prd,     &
                &        iparm_prd, msglvl_prd, b_prd, x_prd, error_prd)
        IF (error_prd /= 0)  THEN
            PRINT*, "pardiso: error_prd =", error_prd
            STOP
        END IF

        do iz=1,nz
            pot_hartree(iz)=x_prd(iz,1)
        enddo

        do iz=1,nz-1
            electric_field_new(iz)=(-pot_hartree(iz+1)+pot_hartree(iz))/dz
        enddo
        electric_field_new(nz)=electric_field_new(nz-1)

        !!!!!!!!!!!!! zapis do pliku !!!!!!!!!!!!!!!!!!
        OPEN(1, FILE="data/charge.dat")
        do iz=1,nz
            z=(iz-1)*dz
            write(1, '(200e20.12)') z/fnm2au, charge_trapped(iz)/fne2au   
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/potential.dat")
        do iz=1,nz
        z=(iz-1)*dz
        write(1, '(200e20.12)') z/fnm2au, -pot_hartree(iz)/feV2au 
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/electric_field.dat")
        do iz=1,nz
        z=(iz-1)*dz
        write(1, '(200e20.12)') z/fnm2au, electric_field(iz)*fnm2au/feV2au 
        enddo
        CLOSE(1)

        OPEN(1, FILE="data/epsilon.dat")
        do iz=1,nz
            z=(iz-1)*dz
            write(1, '(200e20.12)') z/fnm2au, permitivity(eps_0,electric_field(iz))
        enddo
        CLOSE(1)

        DEALLOCATE(pot_hartree)
        DEALLOCATE(ia_prd)
        DEALLOCATE(ja_prd)
        DEALLOCATE(a_prd)
        DEALLOCATE(b_prd)
        DEALLOCATE(x_prd)
    END SUBROUTINE POISSON_ZDIRECTION


END MODULE Poisson_Solver

