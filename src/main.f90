PROGRAM MAIN
    USE Poisson_Solver_Mod
    USE DIELECTRIC 
    IMPLICIT NONE
    REAL*8 :: n0_trapped, L_trapped, eps_0, dz, m1, m2, thickness
    INTEGER*4 :: nz
    REAL*8, ALLOCATABLE :: charge_trapped(:)
    REAL*8, ALLOCATABLE :: electric_field(:)
    REAL*8, ALLOCATABLE :: electric_field_new(:)

    n0_trapped = 5.0*1.0e13*fne2D2au
    L_trapped=15*fnm2au
    eps_0=100
    thickness=100.0*fnm2au
    dz = 0.1*fnm2au
    nz = ceiling(thickness/dz)
    m1=0.2
    m2=3.5

    ALLOCATE(charge_trapped(nz))
    ALLOCATE(electric_field(nz))
    ALLOCATE(electric_field_new(nz))

    CALL POISSON_ZDIRECTION_INIT(n0_trapped, L_trapped, eps_0, nz, dz, charge_trapped, electric_field)
    CALL POISSON_ZDIRECTION(electric_field_new, electric_field, charge_trapped, eps_0,  nz, dz)

    DEALLOCATE(charge_trapped)
    DEALLOCATE(electric_field)
    DEALLOCATE(electric_field_new)


END PROGRAM MAIN
