MODULE INDATA
    USE CONSTANTS
    IMPLICIT NONE
    
    ! calculation parameters
    INTEGER*4 :: Nx
    INTEGER*4 :: Ny
    INTEGER*4 :: Nz
    REAL*8 :: dx
    REAL*8 :: dz
    INTEGER*4 :: z0_indx
    INTEGER*4 :: MAX_ITER
    INTEGER*4 :: MAX_ITER_SCF
    REAL*8 :: tol
    REAL*8 :: tol_scf
    REAL*8 :: alfa ! relaxation parameter for poisson solver
    REAL*8 :: dt
    INTEGER*4 :: MAX_TIME

    ! physical parameters
    REAL*8 :: n0_trapped
    REAL*8 :: L_trapped
    REAL*8 :: m1
    REAL*8 :: m2
    INTEGER*4 :: norbital
    REAL*8 :: sigma ! size fo the initial gaussian

    NAMELIST /calculation_parameters/             &
       &  Nx,                                     &
       &  Ny,                                     &
       &  Nz,                                     &
       &  dx,                                     &
       &  dz,                                     &
       &  z0_indx,                                     &
       &  MAX_ITER,                               &
       &  MAX_ITER_SCF,                           &
       &  tol,                                    &
       &  tol_scf,                                &
       &  alfa,                                   &
       & dt,                                      &
       & MAX_TIME

    NAMELIST /physical_parameters/               &
    &  n0_trapped,                               &
    &  L_trapped,                               &
    &  m1,                                       &
    &  m2,                                       &
    &  norbital,                                 &
    &  sigma


    CONTAINS


    SUBROUTINE GET_INDATA(nmlfile)
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: nmlfile

        OPEN (33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
                &   STATUS="OLD")

        Nx = 0
        Ny = 0
        Nz = 0
        dx = 0.0
        dz = 0.0
        z0_indx = 10
        MAX_ITER = 0
        MAX_ITER_SCF = 0
        tol = 0.0
        tol_scf = 0.0
        alfa = 0.0
        dt = 0.0
        MAX_TIME =0

        n0_trapped = 0.0
        L_trapped = 0.0
        m1 = 0.0
        m2 = 0.0
        norbital = 0
        sigma = 0.0

        READ (33, NML=calculation_parameters)
        dx = dx * fnm2au
        dz = dz * fnm2au
        dt = dt * fns2au

        READ (33, NML=physical_parameters)
        n0_trapped = n0_trapped * 1.0e13*fne2D2au
        L_trapped = L_trapped * fnm2au
        sigma = sigma * fnm2au

    END SUBROUTINE GET_INDATA

END MODULE INDATA