!--------------------------------------
!  PROGRAM MONTE CARLO: HARD SPHERE 
!---------------------------------------
program mc_nvt
implicit none
include 'lj_nvt.inc'
double precision en, vir

 call letrero
 call entrada
 call posiciones
 call cantidades
 call informa
 call ceros
 call etot(en,vir)
 call letrero1
 call ajusta
 call widom(0)
 call mc_cycle(en,vir)
 call promedios
 call perfil
 call widom(2)
 call guarda
 
stop
end
!-----------------------
!  LETRERO
!-----------------------

!-----------------------
!  ENTRADA
!-----------------------

!-----------------------
!  POSICIONES
!-----------------------

!-----------------------
!  CANTIDADES
!-----------------------

!-----------------------
!  INFORMA
!-----------------------

!-----------------------
!  CEROS
!-----------------------

!-----------------------
!  ETOT
!-----------------------

!-----------------------
!  LETRERO 1
!-----------------------

!-----------------------
!  AJUSTA
!-----------------------

!-----------------------
!  WIDOM
!-----------------------

!-----------------------
!  MC_CYCLE
!-----------------------

!-----------------------
!  PROMEDIOS
!-----------------------

!-----------------------
!  PERFIL
!-----------------------

!-----------------------
!  GUARDA
!-----------------------
