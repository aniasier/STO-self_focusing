PROGRAM MAIN
USE Poisson_Solver_Mod
IMPLICIT NONE


double precision :: initial_psi
INTEGER :: Nx, Ny, Nz, MAX_ITER, NLAO
INTEGER :: i, j, k
double precision :: dx, x, y, z, x0, y0, z0, sigma, alfa, tol
double precision, ALLOCATABLE :: psi(:,:,:)
double precision, ALLOCATABLE :: density(:,:,:)
double precision, ALLOCATABLE :: potential(:,:,:)
double precision, ALLOCATABLE :: epsilon(:,:,:)
double precision, ALLOCATABLE :: charge_trapped(:,:,:)
double precision :: val
double precision :: L_trapped, n0_trapped

integer(8) :: start_count, end_count, count_rate
real(8)    :: elapsed_time

! has to be 2**n+1
Nx=256
Ny=256
Nz=256
NLAO=0
tol=1e-5
MAX_ITER=10000

alfa=1.5

dx=1.0*fnm2au
sigma=5*fnm2au

n0_trapped=1e20*fne2au
L_trapped=15*fnm2au

x0 = (Nx-1)*dx/2.0d0
y0 = (Ny-1)*dx/2.0d0
z0 = (Nz-1)*dx/2.0d0

allocate(psi(Nx,Ny,Nz))
allocate(density(Nx,Ny,Nz))
allocate(potential(Nx,Ny,Nz))
allocate(epsilon(Nx,Ny,Nz))
allocate(charge_trapped(Nx,Ny,Nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! initialization of state !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1+NLAO, Ny-NLAO
   do i = 1+NLAO, Nx-NLAO
      do k = 1, Nz
      z=(k-1)*dx
      charge_trapped(i,j,k)=(n0_trapped)*dexp(-z/L_trapped)
      enddo
   enddo
enddo


epsilon=300.0
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

!!!!!! normalizacja funkcji falowej
val=0.0
do k = 1, Nz
   do j = 1, Ny
      do i = 1, Nx         
         val = val+abs(psi(i,j,k))**2*dx*dx*dx
      enddo
   enddo
enddo
do k = 1, Nz
   do j = 1, Ny
      do i = 1, Nx         
         psi(i,j,k) = psi(i,j,k)/sqrt(val)
      enddo
   enddo
enddo

!!!!! inital electron concentration
do k = 1, Nz
   do j = 1, Ny
      do i = 1, Nx         
        density(i,j,k) = abs(psi(i,j,k))**2
      enddo
   enddo
enddo

!rzut na x--y i zapis do pliku
OPEN(1, FILE="density_inital.dat")
do i = 1, Nx
do j = 1, Ny
   val=0.0
   do k = 1, Nz
      val=val+density(i,j,k)*dx
   enddo
write(1, '(200e20.12)') (i-1)*dx/fnm2au, (j-1)*dx/fnm2au, val
enddo
write(1, '(200e20.12)') 
enddo
CLOSE(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! Poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call system_clock(start_count)
CALL Poisson_epsilon(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, tol, MAX_ITER, charge_trapped)
call system_clock(end_count, count_rate)
elapsed_time = real(end_count - start_count) / real(count_rate)
print *, "Executuon time:", elapsed_time

OPEN(1, FILE="potential_inital.dat")
do i = 1, Nx
do j = 1, Ny
   val=0.0
   do k = 1, Nz
      val=val-potential(i,j,k)*dx
   enddo
write(1, '(200e20.12)') (i-1)*dx/fnm2au, (j-1)*dx/fnm2au, val/feV2au
enddo
write(1, '(200e20.12)') 
enddo
CLOSE(1)

OPEN(1, FILE="potential_cros_section.dat")
i = Nx/2
j = Ny/2
val=0.0
do k = 1, Nz
   write(1, '(200e20.12)') (k-1)*dx/fnm2au, -potential(i,j,k)/feV2au
enddo
CLOSE(1)


! potential=0.0

! call system_clock(start_count)
! CALL Poisson_Multilevel(potential, density, epsilon, alfa, Nx, Ny, Nz, dx, 2, tol, MAX_ITER)
! call system_clock(end_count, count_rate)
! elapsed_time = real(end_count - start_count) / real(count_rate)
! print *, "Executuon time:", elapsed_time


! OPEN(1, FILE="potential_inital_multi.dat")
! do i = 1, Nx
! do j = 1, Ny
!    val=0.0
!    do k = 1, Nz
!       val=val+potential(i,j,k)*dx
!    enddo
! write(1, '(200e20.12)') (i-1)*dx/fnm2au, (j-1)*dx/fnm2au, val
! enddo
! enddo
! CLOSE(1)

END PROGRAM MAIN




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function initial_psi(x,y,z,x0,y0,z0,sigma) result(val)
    implicit none

    double precision, intent(in) :: x, y, z, x0, y0, z0, sigma
    double precision :: val
    double precision, parameter :: pi = 3.141592653589793d0
    double precision, parameter :: fnm2au = 18.897261d0

    val = 1.0d0 / ((2.0d0*pi)**(3.0d0/2.0d0) * sigma**3) * &
      exp(-((x-x0)**2 + (y-y0)**2 + (z-z0)**2) / (2.0d0*sigma**2))

end function initial_psi


function charge_Vg(V_gate)
    implicit double precision (a-h,o-z)
    fnm2au=18.897261
    fne2D2au=1.0/(1e14*fnm2au**2)
    feV2au=0.03674932587
    
    charge=0.0e13*fne2D2au
    charge_Vg=charge!*(1.0-exp(-V_gate/(200*feV2au)))
    return
end


function permitivity(eps_0,electric_field,T)
    implicit double precision (a-h,o-z)
    fm2au=18.897261*1e9
    feV2au=0.03674932587
    A=4.097*1e-5
    B=4.907*1e-10*(fm2au/feV2au)
    permitivity=eps_0+1.0/(A+B*dabs(electric_field))

    !Inny sposob
     !B=2.55*1e4
     !E0=8.22*1e4*feV2au/fm2au
     !permitivity=1+B/(1.0+(electric_field/E0)**2)**(1.0/3.0)
    
    return
end


function fermi_distribution_2D(E,EF,T,fmx,fmy)
    implicit double precision (a-h,o-z)
    stala_kb=3.1668152e-6
    pi=3.14159265
    energy=E-EF
    
    if(energy.lt.0.0) then 
        dlogarytm = dlog(1.0+exp(energy/stala_kb/T)) - energy/stala_kb/T
    else
        dlogarytm = dlog(1.0+exp(-energy/stala_kb/T))
    endif
    fermi_distribution_2D=dsqrt(fmx*fmy)*stala_kb*T/(2*pi)*dlogarytm
    return
end
	

function fermi_distribution(E,EF,T)
	implicit double precision (a-h,o-z)
	stala_kb=3.1668152e-6
	if(T.ne.0.0) then
	  beta=1.0/(stala_kb*T)
	  fermi_distribution=1.0/(dexp(beta*(E-EF))+1.0)
	else
	  if (E.GE.EF) then
		fermi_distribution=0.0
	  else if(E.eq.EF) then
		fermi_distribution=0.5
	  else
		fermi_distribution=1.0
	  endif
	endif
	return
end


!-----------------------------------------------------------------------
SUBROUTINE mix_broyden( ndim, deltaout, deltain, alphamix, iter, n_iter, conv )
    !-----------------------------------------------------------------------
    !!
    !! Modified Broyden's method for potential/charge density mixing
    !!             D.D.Johnson, PRB 38, 12807 (1988)
    !!
   
    IMPLICIT NONE
    !
    LOGICAL, INTENT (in) :: conv
    !! If true convergence reache
    !
    INTEGER, INTENT (in) :: ndim
    !! Dimension of arrays deltaout, deltain
    INTEGER, INTENT (in) :: iter
    !! Current iteration number
    INTEGER, INTENT (in) :: n_iter
    !! Number of iterations used in the mixing
    !
    REAL*8, INTENT (in) :: alphamix
    !! Mixing factor (0 < alphamix <= 1)
    REAL*8, INTENT (inout) :: deltaout(ndim)
    !! output Delta at current iteration
    REAL*8, INTENT (inout) :: deltain(ndim)
    !! Delta at previous iteration  
    !
    !   Here the local variables
    !
    ! max number of iterations used in mixing: n_iter must be .le. maxter
    INTEGER, PARAMETER :: maxter = 8
    !
    INTEGER ::  n, i, j, iwork(maxter), info, iter_used, ipos, inext 
    ! work space containing info from previous iterations:
    ! must be kept in memory and saved between calls
    REAL*8, ALLOCATABLE, SAVE :: df2(:,:), dv2(:,:)
    !
    REAL*8, ALLOCATABLE :: deltainsave(:)
    REAL*8 :: beta(maxter,maxter), gammamix, work(maxter), norm
    REAL*8, EXTERNAL :: DDOT, DNRM2
    ! adjustable PARAMETERs as suggested in the original paper
    REAL*8 wg(maxter), wg0
    DATA wg0 / 0.01d0 /, wg / maxter * 1.d0 /
    !
    !IF ( iter .lt. 1 ) print *, 'mix_broyden2','n_iter is smaller than 1'
    !IF ( n_iter .gt. maxter ) print *, 'mix_broyden2','n_iter is too big'
    !IF ( ndim .le. 0 ) print *, 'mix_broyden2','ndim .le. 0'
    !
    IF ( iter .eq. 1 ) THEN
       IF ( .not. ALLOCATED(df2) ) ALLOCATE( df2(ndim,n_iter) )    
       IF ( .not. ALLOCATED(dv2) ) ALLOCATE( dv2(ndim,n_iter) )    
    ENDIF
    IF ( conv ) THEN
       IF ( ALLOCATED(df2) ) DEALLOCATE(df2)
       IF ( ALLOCATED(dv2) ) DEALLOCATE(dv2)
       RETURN
    ENDIF
    IF ( .not. ALLOCATED(deltainsave) ) ALLOCATE( deltainsave(ndim) )    
    deltainsave(:) = deltain(:)
    !
    ! iter_used = iter-1  IF iter <= n_iter
    ! iter_used = n_iter  IF iter >  n_iter
    !
    iter_used = min(iter-1,n_iter)
    !
    ! ipos is the position in which results from the present iteraction
    ! are stored. ipos=iter-1 until ipos=n_iter, THEN back to 1,2,...
    !
    ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
    !
    DO n = 1, ndim
       deltaout(n) = deltaout(n) - deltain(n)
    ENDDO
    !
    IF ( iter .gt. 1 ) THEN
       DO n = 1, ndim
          df2(n,ipos) = deltaout(n) - df2(n,ipos)
          dv2(n,ipos) = deltain(n)  - dv2(n,ipos)
       ENDDO
       norm = ( DNRM2( ndim, df2(1,ipos), 1 ) )**2.d0
       norm = sqrt(norm)
       CALL DSCAL( ndim, 1.d0/norm, df2(1,ipos), 1 )
       CALL DSCAL( ndim, 1.d0/norm, dv2(1,ipos), 1 )
    ENDIF
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(i,j) = wg(i) * wg(j) * DDOT( ndim, df2(1,j), 1, df2(1,i), 1 )
       ENDDO
       beta(i,i) = wg0**2.d0 + wg(i)**2.d0
    ENDDO
    !
    ! DSYTRF computes the factorization of a real symmetric matrix 
    !
    CALL DSYTRF('U', iter_used, beta, maxter, iwork, work, maxter, info)
    !CALL errore('broyden', 'factorization', info)
    !
    ! DSYTRI computes the inverse of a real symmetric indefinite matrix
    !
    CALL DSYTRI('U', iter_used, beta, maxter, iwork, work, info)
    !CALL errore('broyden', 'DSYTRI', info)
    !
    DO i = 1, iter_used
       DO j = i + 1, iter_used
          beta(j, i) = beta(i, j)
       ENDDO
    ENDDO
    !
    DO i = 1, iter_used
       work(i) = DDOT( ndim, df2(1,i), 1, deltaout, 1 )
    ENDDO
    !
    DO n = 1, ndim
       deltain(n) = deltain(n) + alphamix * deltaout(n)
    ENDDO
    !
    DO i = 1, iter_used
       gammamix = 0.d0
       DO j = 1, iter_used
          gammamix = gammamix + beta(j,i) * wg(j) * work(j)
       ENDDO
       !
       DO n = 1, ndim
          deltain(n) = deltain(n) - wg(i) * gammamix * ( alphamix * df2(n,i) + dv2(n,i) )
       ENDDO
    ENDDO
    !
    inext = iter - ( ( iter - 1 ) / n_iter) * n_iter
    df2(:,inext) = deltaout(:)
    dv2(:,inext) = deltainsave(:)
    !
    IF ( ALLOCATED(deltainsave) ) DEALLOCATE(deltainsave)
    !
    RETURN
    !
    END SUBROUTINE mix_broyden
  
