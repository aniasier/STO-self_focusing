MODULE WRITERS
    USE CONSTANTS
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE WRITE_DENSITY(density, Nx, Ny, Nz, dx, filename)
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: density(Nx, Ny, Nz), dx
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER*4 :: i, j, k, unit
        REAL*8 :: x, y, z
        
        ! Open file for writing
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        
        ! Write header comment for clarity
        WRITE(unit, '(A)') '# 3D Density Data'
        WRITE(unit, '(A)') '# x y z density'
        
        ! Write density data
        ! Format: separate each z-slice with blank line for gnuplot compatibility
        DO k = 1, Nz
            DO j = 1, Ny
                DO i = 1, Nx
                    x = (i-1) * dx/fnm2au
                    y = (j-1) * dx/fnm2au
                    z = (k-1) * dx/fnm2au
                    WRITE(unit, '(4E15.7)') x, y, z, density(i,j,k)
                END DO
                WRITE(unit, '(A)') ''  ! Blank line for gnuplot
            END DO
            WRITE(unit, '(A)') ''  ! Extra blank line between z-slices
        END DO
        
        CLOSE(unit)
        
    END SUBROUTINE WRITE_DENSITY

    SUBROUTINE WRITE_DENSITY_2D_XY(density, Nx, Ny, Nz, dx, dz, filename)
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: density(Nx, Ny, Nz), dx, dz
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER*4 :: i, j, k, unit
        REAL*8 :: x, y, val
        
        ! Open file for writing
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        
        ! Write header comment for clarity
        WRITE(unit, '(A)') '# 2D Density Projection onto X-Y plane'
        WRITE(unit, '(A)') '# x y integrated_density'
        
        ! Write 2D projection (integrate over z)
        DO i = 1, Nx
            DO j = 1, Ny
                val = 0.0d0
                ! Integrate over z dimension
                DO k = 1, Nz
                    val = val + density(i,j,k) * dz
                END DO
                x = (i-1) * dx / fnm2au
                y = (j-1) * dx / fnm2au
                WRITE(unit, '(3G25.16)') x, y, val
            END DO
            WRITE(unit, '(A)') ''  ! Blank line for gnuplot
        END DO
        
        CLOSE(unit)
        
    END SUBROUTINE WRITE_DENSITY_2D_XY

    SUBROUTINE WRITE_POTENTIAL_2D_XY(potential, Nx, Ny, Nz, dx,dz, filename)
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: potential(Nx, Ny, Nz), dx, dz
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER*4 :: i, j, k, unit
        REAL*8 :: x, y, val
        
        ! Open file for writing
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        
        ! Write header comment for clarity
        WRITE(unit, '(A)') '# 2D Potential Projection onto X-Y plane'
        WRITE(unit, '(A)') '# x y integrated_potential'
        
        ! Write 2D projection (integrate over z)
        DO i = 1, Nx
            DO j = 1, Ny
                val = 0.0d0
                ! Integrate over z dimension
                DO k = 1, Nz
                    val = val + potential(i,j,k) * dz
                END DO
                x = (i-1) * dx / fnm2au
                y = (j-1) * dx / fnm2au
                WRITE(unit, '(200e20.12)') x, y, val / feV2au
            END DO
            WRITE(unit, '(200e20.12)')
        END DO
        
        CLOSE(unit)
        
    END SUBROUTINE WRITE_POTENTIAL_2D_XY

    SUBROUTINE WRITE_POTENTIAL_CROSS_SECTION(potential, Nx, Ny, Nz, dx, filename)
        IMPLICIT NONE
        INTEGER*4, INTENT(IN) :: Nx, Ny, Nz
        REAL*8, INTENT(IN) :: potential(Nx, Ny, Nz), dx
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER*4 :: i, j, k, unit
        REAL*8 :: z
        
        ! Open file for writing
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        
        ! Write header comment for clarity
        WRITE(unit, '(A)') '# Potential Cross-Section at center (i=Nx/2, j=Ny/2)'
        WRITE(unit, '(A)') '# z potential'
        
        ! Write cross-section at center
        i = Nx / 2
        j = Ny / 2
        DO k = 1, Nz
            z = (k-1) * dx / fnm2au
            WRITE(unit, '(200e20.12)') z, -potential(i,j,k) / feV2au
        END DO
        
        CLOSE(unit)
        
    END SUBROUTINE WRITE_POTENTIAL_CROSS_SECTION

END MODULE WRITERS