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

!zmienne 
REAL*8 :: pot_chem, T, d_pot_chem, el_density, n_el_density
REAL*8 :: thickness, thickness_z
REAL*8 :: lever_arm
REAL*8 :: dz
INTEGER :: nz
CHARACTER*50 text
REAL*8 :: n0_trapped, L_trapped, eps_0, V_gate, L_tot, charge_bc
CHARACTER gate
double precision :: permitivity
REAL*8 :: m1, m2
INTEGER :: n_m1, n_m2

INTEGER :: n_base
REAL*8, ALLOCATABLE :: pot(:)
REAL*8, ALLOCATABLE :: ne(:,:)
REAL*8, ALLOCATABLE :: psi_base(:,:,:)
REAL*8, ALLOCATABLE :: energy_el(:,:)


REAL*8 :: tmp
INTEGER :: nargc
REAL*8 :: V_gate_1, n0_trapped_tmp
!!!!!!!!!!!!!!!!!!!obsuga argumentow wywolania!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nargc=iargc()
if(nargc.lt.1) then
  write(*,*) "Brak argumentow wywolania"
  STOP
else
  call getarg(1,text)
  read(text,*) n0_trapped_tmp
endif

m1=0.2
m2=3.5

thickness_z=100.0
L_tot=thickness_z*fnm2au
thickness=thickness_z*fnm2au
write(*,*) "Grubosc:", thickness_z

dz=0.1*fnm2au
n_base=50
nz=ceiling(thickness/dz)

T=0.1
pot_chem=0.0*feV2au
el_density=4.55*1.0e13*fne2D2au

n0_trapped=n0_trapped_tmp*1.0e13*fne2D2au
L_trapped=15*fnm2au
lever_arm=0.001
eps_0=100
V_gate=20.0*feV2au
gate = 'B'

allocate(pot(nz)) 
allocate(energy_el(n_base,2)) 
allocate(ne(nz,2))
allocate(psi_base(nz,n_base,2)) 

if (gate.eq.'T') then
    CALL SP(el_density, n0_trapped, L_trapped, L_tot, eps_0, V_gate, nz, dz, n_base, &
    & T, pot_chem, pot, ne, psi_base, energy_el, gate, lever_arm, m1, m2)
else 
    CALL SP(el_density, n0_trapped, L_trapped, L_tot, eps_0, V_gate, nz, dz, n_base, &
    & T, pot_chem, pot, ne, psi_base, energy_el, gate, lever_arm, m1, m2)
endif


END PROGRAM MAIN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Schrodinger Poisson method for nanofilms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SP(el_density, n0_trapped, L_trapped, L_tot, eps_0, V_gate, nz, dz, n_base, &
             T, pot_chem, pot, ne, psi, energy_el, gate, lever_arm, m1, m2)
IMPLICIT NONE
!parametry
REAL*8, PARAMETER :: fnm2au=18.897261
REAL*8, PARAMETER :: feV2au=0.03674932587
REAL*8, PARAMETER :: pi=3.141592
REAL*8, PARAMETER :: fne2au=1.0/(1e21*fnm2au**3)
REAL*8, PARAMETER :: kb=3.1668152e-6	
REAL*8, PARAMETER :: fne2D2au=1.0/(1e14*fnm2au**2)
REAL*8, PARAMETER :: epsilon0=1.0/(4.0*pi)

CHARACTER, intent(in) :: gate
REAL*8, intent(in) :: V_gate
REAL*8, intent(in) :: eps_0
REAL*8, intent(in) :: L_trapped
REAL*8, intent(in) :: L_tot
REAL*8, intent(in) :: n0_trapped
REAL*8, intent(in) :: el_density
REAL*8, intent(in) :: dz, T
INTEGER, intent(in) :: nz, n_base
REAL*8, intent(inout) :: pot_chem
REAL*8, dimension(nz), intent(inout) :: pot
REAL*8, dimension(nz,2), intent(out) :: ne
REAL*8, dimension(nz,n_base,2), intent(out) :: psi
REAL*8, dimension(n_base,2), intent(out) :: energy_el
REAL*8, intent (in) :: lever_arm
REAL*8, intent (in) :: m1, m2

!zmienne pomocnicze
REAL*8 :: suma, tmp, s_diff, diff_sp
double precision :: fermi_distribution, fermi_distribution_2D, permitivity, charge_Vg
REAL*8 :: al_sp
INTEGER :: i,j,iz, ii, jj, itmp
REAL*8 :: pot_chem_a, pot_chem_b, pot_chem_c, d_pot_chem
REAL*8 :: el_density_a,	el_density_b, el_density_c
REAL*8 :: mx_xy, mx_xz, mx_yz
REAL*8 :: my_xy, my_xz, my_yz
REAL*8 :: mz_xy, mz_xz, mz_yz
REAL*8 :: charge_bc

integer iterator_sp, iterator_poisson
REAL*8 :: tol_ne, tol_sp		
REAL*8 :: z

REAL*8, ALLOCATABLE :: charge_trapped(:)
REAL*8, ALLOCATABLE :: electric_field(:)
REAL*8, ALLOCATABLE :: pot_hartree_old(:)
REAL*8, ALLOCATABLE :: pot_hartree(:)
REAL*8, ALLOCATABLE :: pot_Vbg(:)

double precision, ALLOCATABLE :: A(:,:)
double precision, ALLOCATABLE :: VR(:,:)
double precision, ALLOCATABLE :: VL(:,:)
double precision, ALLOCATABLE :: WR(:)
double precision, ALLOCATABLE :: WI(:)
double precision, ALLOCATABLE :: WORK(:)
INTEGER  LWORK,INFO
CHARACTER*50 text

!PARDISO
INTEGER, ALLOCATABLE :: nmat(:)        ! row/column index from positions
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

!parametry do metody Broydena dla delty
INTEGER :: n_iter_broyden
INTEGER :: ndim_broyden
REAL*8, ALLOCATABLE :: broyden(:)
REAL*8, ALLOCATABLE :: broyden_p(:)
LOGICAL :: conv


!masy efektywne dla poszczegolnych pasm 
mx_xy=m1
my_xy=m1
mz_xy=m2

mx_xz=m1
my_xz=m2 
mz_xz=m1

mx_yz=m2
my_yz=m1
mz_yz=m1

d_pot_chem=0.2*feV2au
tol_ne=0.0001e11*fne2D2au

al_sp=0.5
iterator_sp=0
tol_sp=1e-5*feV2au

!Broyden method
ndim_broyden=nz
n_iter_broyden=4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    tablice     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(electric_field(nz))
allocate(charge_trapped(nz))
allocate(pot_hartree_old(nz))
allocate(pot_hartree(nz))
allocate(pot_Vbg(nz))

!zakladamy warunki brzegowe Dirichleta stad nz-2
allocate(A(nz-2,nz-2))
allocate(VR(nz-2,nz-2))
allocate(VL(nz-2,nz-2))
allocate(WR(nz-2))
allocate(WI(nz-2))
LWORK=4*(nz-2)
allocate(WORK(LWORK))

!broyden
ALLOCATE(broyden(ndim_broyden))
ALLOCATE(broyden_p(ndim_broyden))
  
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

nrhs_prd=1
n_prd=nz
maxnonzeroprd=(n_prd-2)*3+2+1

ALLOCATE( perm_prd(n_prd) )
ALLOCATE( ia_prd(n_prd+1) )
ALLOCATE( ja_prd(maxnonzeroprd) )
ALLOCATE( a_prd(maxnonzeroprd) )
ALLOCATE( b_prd(n_prd,nrhs_prd) )
ALLOCATE( x_prd(n_prd,nrhs_prd) )

pt_prd(:)= 0
maxfct_prd= 1
mnum_prd= 1
mtype_prd= 11   ! real nonsymmetric  
phase_prd= 13
msglvl_prd= 0
iparm_prd(:)= 0

perm_prd(:)= 0
ia_prd(:)= 0
ja_prd(:)= 0
a_prd(:)= 0.
b_prd(:,:)= 0.
x_prd(:,:)= 0.


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
OPEN(1, FILE="charge_trapped.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(1, '(200e20.12)') z/fnm2au, charge_trapped(iz)/fne2au   
enddo
CLOSE(1)

OPEN(1, FILE="potential_trapped.dat")
do iz=1,nz
   z=(iz-1)*dz
   write(1, '(200e20.12)') z/fnm2au, -pot_hartree(iz)/feV2au 
enddo
CLOSE(1)

OPEN(1, FILE="electric_field_trapped.dat")
do iz=1,nz
   z=(iz-1)*dz
   write(1, '(200e20.12)') z/fnm2au, electric_field(iz)*fnm2au/feV2au 
enddo
CLOSE(1)

OPEN(1, FILE="epsilon_trapped.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(1, '(200e20.12)') z/fnm2au, permitivity(eps_0,electric_field(iz),T)
enddo
CLOSE(1)


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!METODA SAMOUZGODNIONA SP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pot_hartree_old=0.0	
iterator_sp=1
diff_sp=1e6
conv=.FALSE.

!normalizacja do zera dla z=0
tmp=pot_hartree(1)
do iz=1,nz
    pot_hartree(iz)=pot_hartree(iz)-tmp
enddo

do iz=1,nz
    broyden(iz)=pot_hartree(iz)
    broyden_p(iz)=pot_hartree(iz)
    pot(iz)=-pot_hartree(iz)
    pot_hartree_old(iz)=pot_hartree(iz)
enddo

!petla po S-P
do while (.not.conv)
        
        do iz=1,nz-1
            electric_field(iz)=(pot(iz+1)-pot(iz))/dz
        enddo
        electric_field(nz)=electric_field(nz-1)
             
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!! Obliczenia dla pasma xy !!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        do i=1,nz-2
            WR(i)=0.0
            WI(i)=0.0	
            do j=1,nz-2
                A(i,j)=0.0
                VL(i,j)=0.0
                VR(i,j)=0.0
            enddo
        enddo
        do i=1,LWORK
            WORK(i)=0.0
        enddo
        
        do i=1,nz-2
            if (i.eq.1) then
                A(i,i)=1/mz_xy/dz**2+pot(i+1)
                A(i,i+1)=-1.0/2.0/mz_xy/dz**2
            else if (i.eq.(nz-2)) then
                A(i,i)=1/mz_xy/dz**2+pot(i+1)
                A(i,i-1)=-1.0/2.0/mz_xy/dz**2
            else
                A(i,i)=1/mz_xy/dz**2+pot(i+1)         
                A(i,i+1)=-1.0/2.0/mz_xy/dz**2
                A(i,i-1)=-1.0/2.0/mz_xy/dz**2
            endif                    
        enddo

        CALL DGEEV( 'N', 'V', nz-2, A, nz-2, WR, WI, VL, nz-2, VR, nz-2, WORK, LWORK, INFO )
        !write(*,*) "INFO=",INFO

        !sortowanie Energii
        do i=1,nz-2
            do j=1,nz-2-1
            if(WR(j).gt.WR(j+1)) then
            !energie
            tmp=WR(j)
            WR(j)=WR(j+1)
            WR(j+1)=tmp

            !wartosci wlasne
            do ii=1,nz-2
                tmp=VR(ii,j)
                VR(ii,j)=VR(ii,j+1)
                VR(ii,j+1)=tmp
            enddo
            endif
            enddo
        enddo
        
        !funkcje falowe
        psi(:,:,1)=0.0
        
        do i=1,n_base
            do iz=1,nz-2
                    psi(iz+1,i,1)=VR(iz,i)
            enddo
        enddo
        !warunki brzegowe
        do i=1,n_base
        psi(1,i,1)=0.0
        psi(nz,i,1)=0.0
        enddo

        !energie
        do i=1,n_base
            energy_el(i,1)=WR(i)
        enddo
        
        !normalizacja
        do i=1,n_base
            tmp=0.0
            do iz=1,nz
                tmp=tmp+dabs(psi(iz,i,1))**2*dz
            enddo
            do iz=1,nz
                psi(iz,i,1)=1.0/(dsqrt(tmp))*psi(iz,i,1)
            enddo
        enddo

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!! Obliczenia dla pasma xz yz !!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        do i=1,nz-2
            WR(i)=0.0
            WI(i)=0.0	
            do j=1,nz-2
                A(i,j)=0.0
                VL(i,j)=0.0
                VR(i,j)=0.0
            enddo
        enddo
        do i=1,LWORK
            WORK(i)=0.0
        enddo
        
        do i=1,nz-2
            if (i.eq.1) then
                A(i,i)=1/mz_xz/dz**2+pot(i+1)
                A(i,i+1)=-1.0/2.0/mz_xz/dz**2
            else if (i.eq.(nz-2)) then
                A(i,i)=1/mz_xz/dz**2+pot(i+1)
                A(i,i-1)=-1.0/2.0/mz_xz/dz**2
            else
                A(i,i)=1/mz_xz/dz**2+pot(i+1)         
                A(i,i+1)=-1.0/2.0/mz_xz/dz**2
                A(i,i-1)=-1.0/2.0/mz_xz/dz**2
            endif                    
        enddo

        CALL DGEEV( 'N', 'V', nz-2, A, nz-2, WR, WI, VL, nz-2, VR, nz-2, WORK, LWORK, INFO )
        !write(*,*) "INFO=",INFO

        !sortowanie Energii
        do i=1,nz-2
            do j=1,nz-2-1
            if(WR(j).gt.WR(j+1)) then
            !energie
            tmp=WR(j)
            WR(j)=WR(j+1)
            WR(j+1)=tmp

            !wartosci wlasne
            do ii=1,nz-2
                tmp=VR(ii,j)
                VR(ii,j)=VR(ii,j+1)
                VR(ii,j+1)=tmp
            enddo
            endif
            enddo
        enddo
        
        !funkcje falowe
        psi(:,:,2)=0.0
        do i=1,n_base
            do iz=1,nz-2
                    psi(iz+1,i,2)=VR(iz,i)
            enddo
        enddo
        ! !warunki brzegowe
        do i=1,n_base
        psi(1,i,2)=0.0
        psi(nz,i,2)=0.0
        enddo

        !energie
        do i=1,n_base
            energy_el(i,2)=WR(i)
        enddo
        
        !normalizacja
        do i=1,n_base
            tmp=0.0
            do iz=1,nz
                tmp=tmp+dabs(psi(iz,i,2))**2*dz
            enddo
            do iz=1,nz
                psi(iz,i,2)=1.0/(dsqrt(tmp))*psi(iz,i,2)
            enddo
        enddo            
        

            ! WRITE(text,'(a7,I3.3,a4)') 'energia',iterator_sp,'.dat'
            ! OPEN(9, FILE=text)
            ! do i=1,n_base
            !     write(9, '(I3.3 )', advance="no") i
            !     write(9, '(200e20.12)') energy_el(i)/feV2au
            ! enddo
            ! CLOSE(9)
            
            ! WRITE(text,'(a3,I3.3,a4)') 'psi',iterator_sp,'.dat'
            ! OPEN(10, FILE=text)
            ! do iz=0,nz
            !     z=iz*dz
            !     write(10, '(200e20.12)', advance="no") z/fnm2au
            !     do i=1,n_base
            !         write(10, '(200e20.12)', advance="no") psi(iz,i)
            !     enddo
            !     write(10,*)
            ! enddo
            ! CLOSE(10)
    
            
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!! Wyznaczanie potencjalu chem. metoda bisekcji !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !write(*,*) "Obliczanie potencjalu chemicznego..."
        pot_chem_a=pot_chem-d_pot_chem
        pot_chem_b=pot_chem+d_pot_chem
        
        el_density_a=0.0
        el_density_b=0.0
        el_density_c=0.0
        

        do while (dabs(dabs(el_density_a)-el_density).gt.tol_ne.or.dabs(dabs(el_density_b)-el_density).gt.tol_ne)
            
            ne=0.0
            do iz=1,nz
                !pasm xy
                do i=1,n_base
                    ne(iz,1)=ne(iz,1)+2.0*psi(iz,i,1)**2*fermi_distribution_2D(energy_el(i,1),pot_chem_a,T,mx_xy,my_xy)
                enddo
                !pasm xz, yz
                do i=1,n_base
                    ne(iz,2)=ne(iz,2)+2.0*2.0*psi(iz,i,2)**2*fermi_distribution_2D(energy_el(i,2),pot_chem_a,T,mx_xz,my_xz)
                enddo
            enddo
            
            el_density_a=0.0
            do iz=1,nz
                el_density_a=el_density_a+(ne(iz,1)+ne(iz,2))*dz
            enddo
        
            
            ne=0.0
            do iz=1,nz
                !pasm xy
                do i=1,n_base
                    ne(iz,1)=ne(iz,1)+2.0*psi(iz,i,1)**2*fermi_distribution_2D(energy_el(i,1),pot_chem_b,T,mx_xy,my_xy)
                enddo
                !pasm xz, yz
                do i=1,n_base
                    ne(iz,2)=ne(iz,2)+2.0*2.0*psi(iz,i,2)**2*fermi_distribution_2D(energy_el(i,2),pot_chem_b,T,mx_xz,my_xz)
                enddo
            enddo
            
            el_density_b=0.0
            do iz=1,nz
                el_density_b=el_density_b+(ne(iz,1)+ne(iz,2))*dz
            enddo
            
                    
            ne=0.0
            pot_chem_c=(pot_chem_b+pot_chem_a)/2.0
            
            ne=0.0
            do iz=1,nz
                !pasm xy
                do i=1,n_base
                    ne(iz,1)=ne(iz,1)+2.0*psi(iz,i,1)**2*fermi_distribution_2D(energy_el(i,1),pot_chem_c,T,mx_xy,my_xy)
                enddo
                !pasm xz, yz 
                do i=1,n_base
                    ne(iz,2)=ne(iz,2)+2.0*2.0*psi(iz,i,2)**2*fermi_distribution_2D(energy_el(i,2),pot_chem_c,T,mx_xz,my_xz)
                enddo
            enddo
            
            el_density_c=0.0
            do iz=1,nz
                el_density_c=el_density_c+(ne(iz,1)+ne(iz,2))*dz
            enddo
                
            
             !write(*,*) pot_chem_a/feV2au, pot_chem_b/feV2au, pot_chem_c/feV2au
             !write(*,*) el_density_a/fne2D2au, el_density_b/fne2D2au, el_density_c/fne2D2au
            
            if((dabs(el_density_a)-el_density).lt.0.0.and.(dabs(el_density_b)-el_density).lt.0.0) then
                pot_chem_b=pot_chem_b+d_pot_chem
                el_density_a=0.0
                el_density_b=0.0
                el_density_c=0.0
                cycle
            else if((dabs(el_density_a)-el_density).gt.0.0.and.(dabs(el_density_b)-el_density).gt.0.0) then
                pot_chem_a=pot_chem_a-d_pot_chem
                el_density_a=0.0
                el_density_b=0.0
                el_density_c=0.0
                cycle
            endif
            
            if((dabs(el_density_a)-el_density)*(dabs(el_density_c)-el_density).lt.0.0) then
                pot_chem_b=pot_chem_c
            else if ((dabs(el_density_b)-el_density)*(dabs(el_density_c)-el_density).lt.0.0) then
                pot_chem_a=pot_chem_c
            endif
                
            !stop
        enddo
                    
        pot_chem=pot_chem_a
                    
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!! rozwiazanie rownania Poisson - PARDISO ze zmiennym E !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        charge_bc=0.0
        do iz=1,nz
            charge_bc=charge_bc+(charge_trapped(iz)+ne(iz,1)+ne(iz,2))*dz
        enddo
	
        perm_prd(:)= 0
        ia_prd(:)= 0
        ja_prd(:)= 0
        a_prd(:)= 0.
        b_prd(:,:)= 0.
        x_prd(:,:)= 0.

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
                a_prd(nelem)= 0.5*(permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i-1),T)) 
                nelem= nelem + 1

                ja_prd( nelem )=  i
                a_prd( nelem )= -0.5*(permitivity(eps_0,electric_field(i+1),T)+ &
                            & 2.0*permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i-1),T)) 
                nelem= nelem + 1

                ja_prd( nelem )=i+1
                a_prd( nelem )= 0.5*(permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i+1),T)) 
                nelem= nelem + 1
            ENDIF
            enddo


        ia_prd( n_prd+1 )= nelem
        IF (nelem-1 /= maxnonzeroprd) STOP "PARDISO:  nelem /= maxnonzeroprd"

        !wektor wyrazow wolnych
        DO i=2,nz-1
               b_prd(i, 1)= -(-charge_trapped(i)-ne(i,1)-ne(i,2))*dz*dz/epsilon0
        ENDDO
        
        !warunki brzegowe
        b_prd( 1, 1 )= -charge_bc*dz/permitivity(eps_0,electric_field(1),T)/epsilon0
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

        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!! rozwiazanie rownania Poisson tylko dla napiecia bramki !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        perm_prd(:)= 0
        ia_prd(:)= 0
        ja_prd(:)= 0
        a_prd(:)= 0.
        b_prd(:,:)= 0.
        x_prd(:,:)= 0.

        !uzupelnianie macierzy
            nelem=1
            DO i=1,nz
            IF (i.eq.1) THEN
                ia_prd(i)= nelem

                ja_prd( nelem )=  i
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
                a_prd(nelem)= 0.5*(permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i-1),T)) 
                nelem= nelem + 1

                ja_prd( nelem )=  i
                a_prd( nelem )= -0.5*(permitivity(eps_0,electric_field(i+1),T)+ &
                            & 2.0*permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i-1),T)) 
                nelem= nelem + 1

                ja_prd( nelem )=i+1
                a_prd( nelem )= 0.5*(permitivity(eps_0,electric_field(i),T)+permitivity(eps_0,electric_field(i+1),T)) 
                nelem= nelem + 1
            ENDIF
            enddo


        ia_prd( n_prd+1 )= nelem
        IF (nelem-1 /= maxnonzeroprd-1) STOP "PARDISO:  nelem /= maxnonzeroprd"

        !wektor wyrazow wolnych
        DO i=2,nz-1
               b_prd(i, 1)= 0
        ENDDO
        
        !warunki brzegowe
        b_prd( 1, 1 )= 0.0
        b_prd( nz, 1 )= V_gate*lever_arm

        CALL pardiso (pt_prd, maxfct_prd, mnum_prd, mtype_prd, phase_prd,      &
                &        n_prd, a_prd, ia_prd, ja_prd, perm_prd, nrhs_prd,     &
                &        iparm_prd, msglvl_prd, b_prd, x_prd, error_prd)
        IF (error_prd /= 0)  THEN
            PRINT*, "pardiso: error_prd =", error_prd
            STOP
        END IF

        do iz=1,nz
            pot_Vbg(iz)=x_prd(iz,1)
        enddo


        !normalizacja potnecjalu hartree
        tmp=pot_hartree(1)
        do iz=1,nz    
            pot_hartree(iz)=pot_hartree(iz)-tmp
        enddo

        do iz=1,nz    
            pot_hartree(iz)=pot_hartree(iz)+pot_Vbg(iz)
        enddo

        diff_sp=0.0
        do iz=1,nz
            diff_sp=diff_sp+dabs(pot_hartree(iz)-pot_hartree_old(iz))/nz
        enddo

        if (diff_sp.lt.tol_sp) then
            conv=.TRUE.
        endif        
        write(*,*) iterator_sp, diff_sp/feV2au, conv, pot_chem/feV2au, el_density_a/fne2D2au
    
        !metoda Broydena dla delty
        do iz=1,nz
            broyden(iz)=pot_hartree(iz)
        enddo

        CALL mix_broyden( ndim_broyden, broyden, broyden_p, al_sp, iterator_sp, n_iter_broyden, conv )

        do iz=1,nz
            pot_hartree(iz)=broyden_p(iz)
        enddo

        do iz=1,nz
            pot(iz)=-pot_hartree(iz)
            pot_hartree_old(iz)=pot_hartree(iz) 
        enddo

        iterator_sp=iterator_sp+1

    !zamkniecie glowniej petli S-P
    enddo




!zapis do pliku
OPEN(554, FILE="Fermi_energy_d.dat")
write(554,'(200e20.12)') el_density/fne2D2au, (pot_chem-pot(1))/feV2au
CLOSE(554)

OPEN(1, FILE="pot.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(1, '(200e20.12)') z/fnm2au, (pot(iz)-pot(1))/feV2au, pot(iz)/feV2au
enddo
CLOSE(1)

OPEN(1, FILE="ne.dat") !kolejno pasma xy, xz(yz)
do iz=1,nz
    z=(iz-1)*dz
    write(1, '(200e20.12)') z/fnm2au, ne(iz,1)/fne2au, ne(iz,2)/fne2au, (ne(iz,1)+ne(iz,2))/fne2au
enddo
CLOSE(1)

OPEN(9, FILE="e.dat")
do i=1,n_base
    write(9, '(i10, 200e20.12, 200e20.12)') i, (energy_el(i,1)-pot_chem)/feV2au, (energy_el(i,2)-pot_chem)/feV2au
enddo
CLOSE(9)

OPEN(10, FILE="psi_xy.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(10, '(200e20.12)', advance="no") z/fnm2au
    do i=1,n_base
        write(10, '(200e20.12)', advance="no") psi(iz,i,1)
    enddo
    write(10,*)
enddo
CLOSE(10)

OPEN(10, FILE="psi_xz.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(10, '(200e20.12)', advance="no") z/fnm2au
    do i=1,n_base
        write(10, '(200e20.12)', advance="no") psi(iz,i,2)
    enddo
    write(10,*)
enddo
CLOSE(10)

OPEN(10, FILE="electric_field.dat")
do iz=1,nz
    z=(iz-1)*dz
    write(10, '(200e20.12)') z/fnm2au, electric_field(iz)
enddo
CLOSE(10)

deallocate(pot_hartree_old)
deallocate(pot_hartree)
deallocate(charge_trapped)
deallocate(electric_field)

deallocate(A)
deallocate(VR)
deallocate(VL)
deallocate(WR)
deallocate(WI)
deallocate(WORK)
		
END SUBROUTINE SP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
