PROGRAM MAIN
IMPLICIT NONE
!parametry
REAL*8, PARAMETER :: fnm2au=18.897261
REAL*8, PARAMETER :: feV2au=0.03674932587
REAL*8, PARAMETER :: pi=3.14159265359
REAL*8, PARAMETER :: fne2au=1.0/(1e21*fnm2au**3)
REAL*8, PARAMETER :: kb=3.1668152e-6
REAL*8, PARAMETER :: fne2D2au=1.0/(1e14*fnm2au**2)
REAL*8, PARAMETER :: epsilon0=1.0/(4.0*pi)

CALL SP()


END PROGRAM MAIN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Schrodinger Poisson method !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SP()

IMPLICIT NONE
!parametry
REAL*8, PARAMETER :: fnm2au=18.897261
REAL*8, PARAMETER :: feV2au=0.03674932587
REAL*8, PARAMETER :: pi=3.141592
REAL*8, PARAMETER :: fne2au=1.0/(1e21*fnm2au**3)
REAL*8, PARAMETER :: kb=3.1668152e-6	
REAL*8, PARAMETER :: fne2D2au=1.0/(1e14*fnm2au**2)
REAL*8, PARAMETER :: epsilon0=1.0/(4.0*pi)

REAL*8, ALLOCATABLE :: pos(:)
REAL*8, ALLOCATABLE :: psi(:)
REAL*8, ALLOCATABLE :: psi_old(:)
REAL*8, ALLOCATABLE :: potential(:)
REAL*8, ALLOCATABLE :: potential_old(:)

REAL*8 :: initial_psi
REAL*8 :: length, dx, xp, yp, sigma, d, m, dfi
INTEGER :: N, i, j, ip, jp, ii, nargc
REAL*8 :: epsilon, fi
REAL*8 :: tmp, abs_tol, var_tol, alfa, tol

!to solve eigenproblem
double precision, ALLOCATABLE :: A(:,:)
double precision, ALLOCATABLE :: VR(:,:)
double precision, ALLOCATABLE :: VL(:,:)
double precision, ALLOCATABLE :: WR(:)
double precision, ALLOCATABLE :: WI(:)
double precision, ALLOCATABLE :: WORK(:)
INTEGER  LWORK,INFO
INTEGER :: iterator, max_iteration, index
CHARACTER*50 text
LOGICAL :: conv

REAL*8 :: BC1, BC2

nargc=iargc()
if(nargc.lt.1) then
 write(*,*) "Brak argumentow wywolania"
 STOP
else
 call getarg(1,text)
 read(text,*) sigma
endif

sigma=sigma*fnm2au
m=0.28
!sigma=3*fnm2au
d=5*fnm2au
epsilon=1
length = 500*fnm2au
N=500 
dx=2*length/(2*N)
alfa=0.4
max_iteration=200
tol=1e-4
dfi=0.01

allocate(psi(-N:N))
allocate(psi_old(-N:N))
allocate(pos(-N:N))
allocate(potential(-N:N))
allocate(potential_old(-N:N))

!zakladamy warunki brzegowe Dirichleta stad nz-2
allocate(A(2*N+1,2*N+1))
allocate(VR(2*N+1,2*N+1))
allocate(VL(2*N+1,2*N+1))
allocate(WR(2*N+1))
allocate(WI(2*N+1))
LWORK=4*(2*N+1)
allocate(WORK(LWORK))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialization !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pos=0.0
potential=0.0
do i=-N,N
    pos(i)=i*dx
enddo

OPEN(1, FILE="psi_init.dat")
do i=-N,N
    write(1, *) pos(i)/fnm2au, initial_psi(pos(i),sigma)
enddo
CLOSE(1) 

tmp=0.0
do i=-N,N
   tmp=tmp+(initial_psi(pos(i),sigma))**2*dx
enddo
print *, 'Initial function normalization', tmp

do i=-N,N
   psi(i)=initial_psi(pos(i),sigma)
enddo

potential=0.0
do i=-N,N
   do ip=-N,N
      potential(i)=potential(i)+1/epsilon*abs(psi(ip))**2/sqrt((pos(i)-pos(ip))**2+4*d**2)*dx 
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! self-consistent procedure !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
iterator=1
conv=.FALSE.
do while (.not.conv .AND. iterator < max_iteration)

   do i=-N,N
      !psi_old(i)=psi(i)
      potential_old(i)=potential(i)
   enddo

!!!!!!!!!!!!  solve Schrodinger eqution !!!!!!!!!!!!
   WR=0.0
   WI=0.0	
   A=0.0
   VL=0.0
   VR=0.0
   WORK=0.0
         
   do i=-N,N
      index=i+(N+1)
      if (i.eq.(-N)) then
            A(index,index)=1/m/dx**2+(-potential(i))
            A(index,index+1)=-1.0/2.0/m/dx**2
      else if (i.eq.N) then
            A(index,index)=1/m/dx**2+(-potential(i))
            A(index,index-1)=-1.0/2.0/m/dx**2
      else
            A(index,index)=1/m/dx**2+(-potential(i))
            A(index,index+1)=-1.0/2.0/m/dx**2
            A(index,index-1)=-1.0/2.0/m/dx**2
      endif                    
   enddo

   CALL DGEEV( 'N', 'V', 2*N+1, A, 2*N+1, WR, WI, VL, N, VR, 2*N+1, WORK, LWORK, INFO )
   !write(*,*) "INFO=",INFO

   !sortowanie Energii
   do i=1,2*N+1
      do j=1,2*N+1-1
      if(WR(j).gt.WR(j+1)) then
      !energie
      tmp=WR(j)
      WR(j)=WR(j+1)
      WR(j+1)=tmp
      !wartosci wlasne
      do ii=1,2*N+1
            tmp=VR(ii,j)
            VR(ii,j)=VR(ii,j+1)
            VR(ii,j+1)=tmp
      enddo
      endif
      enddo
   enddo
         
   !funkcje falowe
   psi=0.0
   do i=-N,N
      index=i+(N+1)
      psi(i)=VR(index,1)
   enddo
   
   !normalizacja
   tmp=0.0
   do i=-N,N
      tmp=tmp+abs(psi(i))**2*dx
   enddo

   do i=-N,N
      psi(i)=psi(i)/sqrt(tmp)
   enddo


   !!!!!!!!!!!!! solve a poisson equation !!!!!!!!!!!!!!
   BC1=0.0
   BC2=0.0
   do ip=-N,N
         BC1=BC1+1/epsilon*abs(psi_old(ip))**2/sqrt((-N*dx-pos(ip))**2+4*d**2)*dx
         BC2=BC2+1/epsilon*abs(psi_old(ip))**2/sqrt((N*dx-pos(ip))**2+4*d**2)*dx
   enddo
   
   potential=0.0
   do i=-N,N
      do ip=-N,N
         potential(i)=potential(i)+1/epsilon*abs(psi(ip))**2/sqrt((pos(i)-pos(ip))**2+4*d**2)*dx 
      enddo
   enddo

   ! !BC
   ! do i=-N,N
   !    potential(i)=potential(i)-(N*dx-i*dx)/(2.0*length)*BC1-(N*dx+i*dx)/(2.0*length)*BC2
   ! enddo

   do i=-N,N
      potential(i)=alfa*potential(i)+(1-alfa)*potential_old(i)
   enddo

   
   WRITE(text,'(a3,i0.0,a4)') 'pot',iterator,'.dat'
   OPEN(1, FILE=text)
   do i=-N,N
      write(1, '(200e20.12)') pos(i)/fnm2au, -potential(i)/feV2au
   enddo
   CLOSE(1)

   WRITE(text,'(a3,i0.0,a4)') 'psi',iterator,'.dat'
   OPEN(1, FILE=text)
   do i=-N,N
      write(1, '(200e20.12)') pos(i)/fnm2au, psi(i)
   enddo
   CLOSE(1) 

   OPEN(1, FILE="energy.dat", position="append")
      write(1, '(200e20.12)') REAL(iterator), WR(1)/feV2au
   CLOSE(1)


   ! do i=-N,N
   !    psi(i)=alfa*psi(i)+(1-alfa)*psi_old(i)
   ! enddo

   var_tol=0.0
   abs_tol=0.0
   do i=-N,N
      ! var_tol=var_tol+abs(abs(psi(i))**2-abs(psi_old(i))**2)
      ! abs_tol=abs_tol+abs(psi(i))**2
      var_tol=var_tol+abs(potential(i)-potential_old(i))
      abs_tol=abs_tol+abs(potential(i))
   enddo

   if ((var_tol/abs_tol).lt.tol) then
      conv=.TRUE.
   endif

   print *, iterator, WR(1)/feV2au, var_tol/abs_tol
   iterator=iterator+1

enddo


OPEN(1, FILE="pot.dat")
do i=-N,N
    write(1, '(200e20.12)') pos(i)/fnm2au, -potential(i)/feV2au
enddo
CLOSE(1)

OPEN(1, FILE="psi.dat")
do i=-N,N
    write(1, '(200e20.12)') pos(i)/fnm2au, psi(i)
enddo
CLOSE(1) 
END SUBROUTINE SP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function initial_psi(x,sigma) result(val)
    implicit none

    real(8), intent(in) :: x, sigma
    real(8) :: val
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: fnm2au = 18.897261d0

    val = 1.0d0/(pi*sigma**2)**(1.0d0/4.0d0) * exp(-(x**2)/(2.0d0*sigma**2))

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
  
