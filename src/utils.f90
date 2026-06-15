MODULE UTILS
    USE CONSTANTS
    IMPLICIT NONE
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

END MODULE UTILS