MODULE UTILS
    USE CONSTANTS
    IMPLICIT NONE
    ! interface do dsyev !!!!
    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info) bind(C,name="dsyev_")
            character(len=1), intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real*8, intent(inout) :: a(lda,*)
            real*8, intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine
    end interface

    CONTAINS

    function initial_psi(x,y,z,x0,y0,z0,sigma) result(val)
        implicit none

        double precision, intent(in) :: x, y, z, x0, y0, z0, sigma
        double precision :: val

        val = 1.0d0 / ((2.0d0*pi)**(3.0d0/2.0d0) * sigma**3) * &
        exp(-((x-x0)**2 + (y-y0)**2 + (z-z0)**2) / (2.0d0*sigma**2))

    end function initial_psi

    SUBROUTINE GET_DENSITY(density, psi, Nx, Ny, Nz)
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) ::Nx, Ny, Nz
        REAL*8, INTENT(IN) :: psi(:,:,:)
        REAL*8, INTENT(OUT) ::density(Nx, Ny, Nz)
        INTEGER*4 :: i, j, k

        !!!!! inital electron concentration
        do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx         
                density(i,j,k) = abs(psi(i,j,k))**2
            enddo
        enddo
        enddo

    END SUBROUTINE GET_DENSITY

    SUBROUTINE GET_NORM(density, Nx, Ny, Nz, dx, dz, norm)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: density(:,:,:)
        REAL*8, INTENT(IN) :: dx, dz
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(OUT) :: norm
        INTEGER*4 :: i,j,k
        norm = 0.d0
        do i=1,Nx
            do j=1,Ny
                do k=1,Nz
                    norm = norm + density(i,j,k)
                enddo
            enddo
        enddo

        norm = norm * dx * dx * dz

    END SUBROUTINE GET_NORM

    SUBROUTINE GET_CHARGE_TRAPPED3D(charge_trapped3D, charge_trapped, nx, ny, nz)
        IMPLICIT NONE
        REAL*8, INTENT(OUT) :: charge_trapped3D(nx, ny, nz)
        REAL*8, INTENT(IN) :: charge_trapped(nz)
        INTEGER*4, INTENT(IN) :: nx, ny, nz
        INTEGER*4 :: i, j, k

        DO k=1, Nz
            DO i=1, nx
                DO j=1, ny
                    charge_trapped3D(i,j,k) = charge_trapped(k)
                END DO
            END DO
        END DO
    END SUBROUTINE GET_CHARGE_TRAPPED3D

    SUBROUTINE GET_INIT_PSI(psi, Nx, Ny, Nz, x0, y0, z0, sigma, dx, dz)
        IMPLICIT NONE
        REAL*8, INTENT(OUT) :: psi(Nx, Ny, Nz)
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: x0, y0, z0, sigma, dx, dz
        INTEGER*4 :: i, j, k
        REAL*8 :: x, y, z, val

        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx
                    x=(i-1)*dx
                    y=(j-1)*dx
                    z=(k-1)*dz
                    psi(i,j,k)=initial_psi(x,y,z,x0,y0,z0,sigma)
                enddo
            enddo
        enddo

        ! calculating the norm
        val=0.0
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx         
                    val = val+abs(psi(i,j,k))**2*dx*dx*dz
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
    END SUBROUTINE GET_INIT_PSI

    SUBROUTINE ADD_POTENTIAL(V0, potential, Nx, Ny, Nz, x0, y0, z0, sigma, dx, dz)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: V0
        REAL*8, INTENT(INOUT) :: potential(Nx, Ny, Nz)
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: x0, y0, z0, sigma, dx, dz
        INTEGER*4 :: i,j,k
        REAL*8 :: r2, x, y, z

        DO k=1,Nz
            z = (k-1)*dz
            DO j=1,Ny
                y = (j-1)*dx
                DO i=1,Nx
                    x = (i-1)*dx
                    r2 = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
                    potential(i,j,k) = -V0 * EXP(-r2/(2.0d0*sigma**2))
                END DO
            END DO
        END DO

    END SUBROUTINE ADD_POTENTIAL

    SUBROUTINE INTERPOLATE_1D(data_coarse, nz_coarse, dz_coarse, &
                          data_fine,   nz_fine,   dz_fine)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nz_coarse, nz_fine
    REAL*8, INTENT(IN) :: dz_coarse, dz_fine
    REAL*8, INTENT(IN) :: data_coarse(nz_coarse)
    REAL*8, INTENT(OUT) :: data_fine(nz_fine)

    INTEGER :: kf, kc
    REAL*8 :: z, alpha

    do kf = 1, nz_fine

        z = (kf-1)*dz_fine

        kc = int(z/dz_coarse) + 1

        if (kc >= nz_coarse) then
            data_fine(kf) = data_coarse(nz_coarse)

        else

            alpha = (z-(kc-1)*dz_coarse)/dz_coarse

            data_fine(kf) = (1.d0-alpha)*data_coarse(kc) + &
                             alpha*data_coarse(kc+1)

        endif

    enddo

    END SUBROUTINE INTERPOLATE_1D

    SUBROUTINE GET_ELECTRON_POTENTIAL_Z(density, nx, ny, nz3d, &
                                    dx, dz3D, z_fine, nz_fine, &
                                    x_center, y_center, potential_e)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nx, ny, nz3d, nz_fine
    REAL*8, INTENT(IN) :: density(nx,ny,nz3d)
    REAL*8, INTENT(IN) :: dx, dz3D
    REAL*8, INTENT(IN) :: z_fine(nz_fine)
    REAL*8, INTENT(IN) :: x_center, y_center

    REAL*8, INTENT(OUT) :: potential_e(nz_fine)

    INTEGER :: i, j, k, iz
    REAL*8 :: x, y, z
    REAL*8 :: rx, ry, rz, r
    REAL*8 :: dV
    REAL*8 :: eps_reg

    dV = dx*dx*dz3D

    ! skala regularizacji dla komórki źródłowej
    eps_reg = dx/2.d0

    potential_e(:) = 0.d0

    DO iz = 1, nz_fine

        DO k = 1, nz3d

            z = (k-1)*dz3D

            DO j = 1, ny

                y = (j-1)*dx

                DO i = 1, nx

                    x = (i-1)*dx

                    rx = x - x_center
                    ry = y - y_center
                    rz = z - z_fine(iz)

                    r = sqrt(rx*rx + ry*ry + rz*rz)

                    ! regularizacja osobliwości
                    r = sqrt(r*r + eps_reg*eps_reg)

                    potential_e(iz) = potential_e(iz) &
                        - density(i,j,k)*dV/r

                END DO
            END DO
        END DO

    END DO

    END SUBROUTINE GET_ELECTRON_POTENTIAL_Z


    subroutine solve_eigenproblem(H, w, num)
        real*8, intent(inout) :: H(:,:)
        real*8, intent(out) :: w(:)
        integer, intent(in) :: num
        real*8, allocatable :: work(:)
        integer :: lwork, info

        lwork = -1
        allocate(work(1))
        call dsyev('V', 'U', num, H, num, w, work, lwork, info)

        if (info /= 0) then
            print *, "Błąd przy query workspace, info=", info
            stop
        end if

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call dsyev('V', 'U', num, H, num, w, work, lwork, info)

        if ( info /= 0) then
            print *, 'Błąd, info=', info
        end if
    end subroutine solve_eigenproblem

END MODULE UTILS