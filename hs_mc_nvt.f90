!----------------------------------------------------------------------
! PROGRAM MONTE CARLO: HARD SPHERE
!----------------------------------------------------------------------
program mc_nvt
use variables
implicit none 

double precision :: en, vir

 call letrero
 call entrada
 call posiciones
 call cantidades
 call informa
 call ceros
 call etot(en,vir)
 call letrero1
 call ajusta
 call mc_cycle(en,vir)
 call promedios
 call perfil
 call guarda

stop
end
!----------------------------------------------------------------------
!  letrero
!----------------------------------------------------------------------
subroutine letrero
use variables

write(*,0001)
write(*,0002)
write(*,0003)

0001  format(/,t13,45('!'))
0002  format(t17,('Monte Carlo NVT:  Hard Sphere'))
0003  format( t13,45('!'),/)

end
!----------------------------------------------------------------------
!  letrero1
!----------------------------------------------------------------------
subroutine letrero1
use variables

write(*,0004)
write(*,0005)
write(*,0006)

0004  format(/,t10,56('-'))
0005  format(t14,'ncic',t22,'epo',t34,'vir',t46,'dr',t57,'acc')
0006  format(  t10,56('-'))

end
!----------------------------------------------------------------------
!  entrada
!----------------------------------------------------------------------
subroutine entrada
use variables
implicit none

open(1,file='run.dat',status='unknown')
 read(1,*)
 read(1,*)
 read(1,*)
 read(1,*)lattice
 read(1,*)structure
 read(1,*)ncyclos
 read(1,*)nmoves
 read(1,*)nprome
 read(1,*)nfile
 read(1,*)nprint
 read(1,*)najusta
 read(1,*)nghost
 read(1,*)
 read(1,*)temp
 read(1,*)rcut
 read(1,*)dr
 read(1,*)succ
 read(1,*)
 read(1,*)sig
 read(1,*)eps
 read(1,*)mass
 read(1,*)iseed
 read(1,*)
close(1)

end
!----------------------------------------------------------------------
!  POSICIONES
!----------------------------------------------------------------------
subroutine posiciones
use variables
implicit none

integer :: i

if(lattice) then
 if(structure) then
  write(6,001)'¿ fcc: nat=4*ncell**3 (ncell: # de caldas) ?'
  read(5,*)ncelda
  nat = 4*ncelda**3
  write(6,002)' las particulas quiere simular son:  ', nat
  write(6,*)
  write(6,003)'¿ cual es la densidad del sistema ?'
  read(5,*)dens
  call fcc
 else
  write(6,002)'¿cc: cuantas particulas quiere simular?'
  read(5,*)nat
  write(6,*)
  write(6,003)'¿ cual es la densidad del sistema ?'
  read(5,*)dens
  call cc
 endif
else
 open (1,file='fort.1',status='unknown')
  read(1,*) dr
  read(1,*) nat,boxx,boxy,boxz
  if (nat > maxnat ) stop 'nat es mayor que maxnat'
  do i = 1, nat
   read(1,*) rx(i), ry(i), rz(i)
  enddo
 close (1)
endif

0001  format(t15,a)
0002  format(t15,a,i7)
0003  format(t15,a,f10.5)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine fcc
use variables
implicit none

integer :: iz, iy, ix, iref, m
double precision :: cell, cell2

!---- fcc 
boxx = ( dble(nat)/dens)**(1.0d0/3.0d0)
boxy = boxx
boxz = boxx
write(6,001)'se construyo celda cubica de lado: ', boxx

!---calcula el lado de la celda unitaria
cell  = boxx/dble(ncelda)
cell2 = 0.50d0*cell

!---construye la celda unitaria
!---subcelda 'a'
rx(1) =  0.0
ry(1) =  0.0
rz(1) =  0.0
!---subcelda 'b'
rx(2) =  cell2
ry(2) =  cell2
rz(2) =  0.0
!---subcelda 'c'
rx(3) =  0.0
ry(3) =  cell2
rz(3) =  cell2
!---subcelda 'd'
rx(4) =  cell2
ry(4) =  0.0
rz(4) =  cell2
!---construye la celda apartir de la celda unitaria
m = 0
do iz = 1, ncelda
 do iy = 1, ncelda
  do ix = 1, ncelda
   do iref = 1, 4
    rx(iref+m) = rx(iref) + cell*dble(ix-1)
    ry(iref+m) = ry(iref) + cell*dble(iy-1)
    rz(iref+m) = rz(iref) + cell*dble(iz-1)
   enddo
   m = m + 4
  enddo
 enddo
enddo

0001  format(t15,a,f10.5,/)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine cc
use variables
implicit none

integer :: i, j, k, number, nplace
double precision :: place, size

!---pone las particulas  en una celda cubica de lado boxx
boxx = (dble(nat)/dens)**(1.0d0/3.0d0)
boxy = boxx
boxz = boxx
write(6,001)'se construyo celda cubica de lado: ', boxx

number = idint(dble(nat)**(1.0/3.0) + 1.0d0)
size   = boxx/dble(number)

nplace = 1

do i = 1,number
 do j = 1,number
  do k = 1,number
   if (nplace <= nat) then
    rx(nplace) = dble(i)*size
    ry(nplace) = dble(j)*size
    rz(nplace) = dble(k)*size
    nplace     = nplace + 1
   endif
  enddo
 enddo
enddo

0001  format(t15,a,f10.5,/)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
 subroutine cantidades
 use variables
 implicit none

 double precision :: r2i, r6i, boxt

vol = boxx*boxy*boxz

boxxinv = 1.0/boxx
boxyinv = 1.0/boxy
boxzinv = 1.0/boxz

hboxx = 0.50*boxx
hboxy = 0.50*boxy
hboxz = 0.50*boxz

rcut = min(rcut,hboxx)
rcut = min(rcut,hboxy)
rcut = min(rcut,hboxz)

beta = 1.0/temp
pi   = 4.0*atan(1.0d0)

rcut2 = rcut*rcut
eps4  = 4.0*eps
eps48 = 48.0*eps
sig2  = sig*sig

volx    = boxy*boxz*delta_r
voly    = boxx*boxz*delta_r
volz    = boxx*boxy*delta_r

boxt    = max(boxx,boxy)
max_box = max(boxt,boxz)

n_bin   = int(max_box/delta_r) + 1

if(n_bin > max_bin)stop 'cambiar max_bin'

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine informa
use variables
implicit none

if (lattice) then
 write(*,001)'inicie de una lattice  '
else
 write(*,001)'inicie de un archivo   '
end if
 write(*,002)'numero de atomos       ',nat
 write(*,003)'densidad               ',dble(nat)/vol
 write(*,003)'temperatura            ',temp
 write(*,003)'radio de corte         ',rcut
 write(*,002)'total de ciclos        ',ncyclos
 write(*,002)'imprime en pantalla    ',nprint
 write(*,002)'ajusta el dr           ',najusta
 write(*,002)'calcula promedios      ',nprome
 write(*,002)'guarda posiciones      ',nfile
 write(*,002)'intentos del pot. quim.',nghost

001   format(t15,a)
002   format(t15,a,i7)
003   format(t15,a,f10.6)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine ceros
use variables
implicit none

integer :: i

ac_en  = 0.0
ac_en2 = 0.0

ac_vr  = 0.0
ac_vr2 = 0.0

ac_pr  = 0.0
ac_pr2 = 0.0

intentos = 0
nacc     = 0

conta_sample = 0
ave_enp      = 0.0
ave_press    = 0.0
ave2_enp     = 0.0
ave2_press   = 0.0

con_pro = 0

rhoxt = 0.0
rhoyt = 0.0
rhozt = 0.0

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine etot(ener,vir)
use variables
implicit none
    
integer :: i, jb
double precision :: ener, vir, eni, viri
double precision :: rxi, ryi, rzi

ener = 0.0
vir  = 0.0

do i = 1,nat-1
 rxi = rx(i)
 ryi = ry(i)
 rzi = rz(i)
 jb  = i + 1
 call eneri(rxi, ryi, rzi, i, jb, eni, viri)
 ener = ener + eni
 vir  = vir + viri
enddo

write(*,004)' energia inicial:      ',ener/dble(nat)
write(*,005)' presion inicial:      ',vir/dble(nat)

004   format(/,t14,a,f10.6)
005   format(  t14,a,f10.6)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
 subroutine eneri(rxi, ryi, rzi, i, jb, en, vir)
 use variables
 implicit none

 integer :: i, j, jb
 double precision :: rxi, ryi, rzi, dx, dy, dz, r2
 double precision :: en, vir, virij, enij

en  = 0.0
vir = 0.0

do j = jb, nat
 if(j/=i) then
  dx = rxi - rx(j)
  dy = ryi - ry(j)
  dz = rzi - rz(j)
  if (dx>hboxx) then
   dx = dx - boxx
  elseif (dx<-hboxx) then
   dx = dx + boxx
  endif
  if (dy>hboxy) then
   dy = dy - boxy
  elseif (dy<-hboxy) then
   dy = dy + boxy
  endif
  if (dz>hboxz) then
   dz = dz - boxz
  elseif (dz<-hboxz) then 
   dz = dz + boxz
  endif
  r2 = dx*dx + dy*dy + dz*dz
  call ener(enij, virij, r2)
  en  = en  + enij
  vir = vir + virij
 endif
enddo

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine ener(en,vir,r2)
use variables
implicit none

double precision :: r2, r2i, r6i, en, vir

if(r2<rcut2) then
en = 1.0E4
vir = eps48*(r6i*r6i-0.5d0*r6i)
else
 en  = 0.0
 vir = 0.0
endif

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine mc_cycle(en,vir)
use variables
implicit none
integer :: icyc, imov
double precision :: en, vir
double precision :: enp, press

do icyc =1, ncyclos
 do imov=1,nmoves
  call mcmove(en,vir)
 enddo
 if(mod(icyc,najusta)==0) call ajusta
 if(mod(icyc,nprint)==0 ) call screen(icyc,en,vir,enp,press)
 if(mod(icyc,nfile )==0 ) call guarda
 if(mod(icyc,nprome)==0 ) then 
  call sample(enp,press)
  call profile
 endif
enddo

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine mcmove(en,vir)
use variables
implicit none
integer :: jb, o
double precision :: en, vir, enn, eno, viro, virn
double precision :: rxn, ryn, rzn, ranf, arg
logical :: aceptado

intentos = intentos + 1
jb = 1
o = int(nat*ranf(iseed)) + 1

call eneri(rx(o), ry(o), rz(o), o, jb, eno, viro)

rxn = rx(o) + (2*ranf(iseed)-1.0d0)*dr
ryn = ry(o) + (2*ranf(iseed)-1.0d0)*dr
rzn = rz(o) + (2*ranf(iseed)-1.0d0)*dr

call eneri(rxn, ryn, rzn, o, jb, enn, virn)

arg = enn - eno

if(arg < 0.0)then
 aceptado = .true.
elseif(ranf(iseed) < exp(-beta*arg)) then
 aceptado = .true.
else
 aceptado = .false.
endif

if(aceptado) then
 nacc = nacc + 1
 en   = en  + arg
 vir  = vir + (virn-viro)
 if (rxn<0   ) rxn = rxn + boxx
 if (rxn>boxx) rxn = rxn - boxx
 if (ryn<0   ) ryn = ryn + boxy
 if (ryn>boxy) ryn = ryn - boxy
 if (rzn<0   ) rzn = rzn + boxz
 if (rzn>boxz) rzn = rzn - boxz
 rx(o) = rxn
 ry(o) = ryn
 rz(o) = rzn
endif

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine screen(i,en,vir,enp,press)
use variables
implicit none
integer :: i
double precision :: en,enp, vir, press, corp, fracc

if (nat/=0) then
 enp   = en/dble(nat)
 press = dens/beta + vir/(3.0*vol)
else
 enp   = 0.0
 press = 0.0
endif

fracc = 100.*dble(nacc)/dble(intentos)

write (6, 0008) i, enp, press, dr, fracc
write (99,0008)i, enp, press, dr, fracc

0008  format(t10,i7,4(f12.6))
end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine sample(enp,press)
use variables
implicit none
double precision :: enp, press

conta_sample = conta_sample + 1

ave_enp   = ave_enp   + enp
ave_press = ave_press + press

ave2_enp   = ave2_enp   + enp*enp
ave2_press = ave2_press + press*press

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine promedios
use variables
implicit none
double precision :: muestras, dev_e, var_e, dev_p, var_p

muestras  = 1.0d0/dble(conta_sample)

ave_enp   = ave_enp*muestras
ave_press = ave_press*muestras

var_e     = ave2_enp*muestras   - ave_enp**2
var_p     = ave2_press*muestras - ave_press**2

if(var_e /= 0) dev_e = sqrt(var_e)
if(var_p /= 0) dev_p = sqrt(var_p)

write(*,0009)
write(*,0010)
write(*,0011)

0009   format(/,t10,45('!'))
0010   format(t13,('--------- r e s u l t a d o s ----------'))
0011   format(  t10,45('!'),/)

write(*,0012)' ave. energy  : ', ave_enp,   ' +/- ', dev_e
write(*,0012)' ave. pressure: ', ave_press, ' +/- ', dev_p
write(*,*)

0012   format(a,f10.6,a,f10.6)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine ajusta
use variables
implicit none
integer :: int_tmp, nacc_tmp
double precision :: dro, frac
save             nacc_tmp, int_tmp

if( (intentos==0) .or. (int_tmp>=intentos)) then
 nacc_tmp = nacc
 int_tmp  = intentos
else
 frac = dble(nacc-nacc_tmp)/dble(intentos-int_tmp)
 dro  = dr
 dr   = dr*abs(frac/(succ/100.0))
 if (dr/dro > 1.5d0)       dr = dro*1.5d0
 if (dr/dro < 0.5d0)       dr = dro*0.5d0
 if (dr     > hboxx/2.d0)  dr = hboxx/2.d0
 nacc_tmp = nacc
 int_tmp  = intentos
endif

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine profile
use variables
implicit none
integer :: i, rx_bin, ry_bin, rz_bin

con_pro = con_pro + 1

nx_bin = 0
ny_bin = 0
nz_bin = 0

do i =1, nat
 rx_bin = rx(i)/delta_r + 1
 ry_bin = ry(i)/delta_r + 1
 rz_bin = rz(i)/delta_r + 1
 nx_bin(rx_bin) = nx_bin(rx_bin) + 1
 ny_bin(ry_bin) = ny_bin(ry_bin) + 1
 nz_bin(rz_bin) = nz_bin(rz_bin) + 1
enddo

do i =1, n_bin
 rhoxt(i) = rhoxt(i) + dble(nx_bin(i))
 rhoyt(i) = rhoyt(i) + dble(ny_bin(i))
 rhozt(i) = rhozt(i) + dble(nz_bin(i))
enddo

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine perfil
use variables
implicit none
integer :: i
double precision :: ir, veces

veces = 1.0/dble(con_pro)

do i =1, n_bin-1
 rhoxt(i) = rhoxt(i)*veces/volx
 rhoyt(i) = rhoyt(i)*veces/voly
 rhozt(i) = rhozt(i)*veces/volz
 ir = dble(i)*delta_r
 write(10,0013)ir, rhoxt(i), rhoyt(i), rhozt(i)
enddo

0013  format(4(1X,f10.6))

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine guarda
use variables
implicit none
integer :: i

write(6,0001)' <<< guardando configuracion >>> '

open (1,file='fort.1',status='unknown')
 write(1,*) dr
 write(1,*) nat,boxx,boxy,boxz
 do i = 1, nat
  write(1,*) rx(i), ry(i), rz(i)
 enddo
close (1)

open (2,file='fort.xyz',status='unknown')
 write(2,*) nat
 write(2,*) 
 do i=1,nat
  write(2,*) 'Li', rx(i), ry(i), rz(i)
 enddo
close (2)

0001  format(t25,a)

end
!----------------------------------------------------------------------
!     random number generator
!----------------------------------------------------------------------
function ranf(idum)
integer :: idum
double precision :: ranf, ran2

ranf = ran2(idum)

end
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
function ran2(idum)
integer :: idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
double precision :: ran2, am, eps, rnmx
parameter (im1 = 2147483563,im2 = 2147483399,am  = 1.0d0/im1,imm1 = im1-1,ia1 = 40014)
parameter (ia2 = 40692,iq1 = 53668,iq2 = 52774,ir1 = 12211,ir2 = 3791,ntab = 32)
parameter (ndiv = 1+imm1/ntab,eps = 1.2d-7, rnmx = 1.0d0-eps)

integer :: idum2, j, k, iv(ntab), iy
save  iv, iy, idum2
data idum2 /123456789/, iv /ntab*0/, iy /0/

if (idum <= 0) then
 idum  = max(-idum,1)
 idum2 = idum
 do j = ntab+8,1,-1
  k   = idum/iq1
  idum= ia1*(idum-k*iq1)-k*ir1
  if (idum < 0) idum = idum+im1
  if (j<=ntab) iv(j) = idum
 enddo
 iy = iv(1)
endif
k   = idum/iq1
idum= ia1*(idum-k*iq1)-k*ir1

if (idum<0) idum =idum+im1

k     = idum2/iq2
idum2 = ia2*(idum2-k*iq2)-k*ir2

if (idum2<0) idum2 = idum2+im2

j    = 1+iy/ndiv
iy   = iv(j)-idum2
iv(j)= idum

if(iy<1) iy = iy+imm1
ran2 = min(am*iy,rnmx)

end
