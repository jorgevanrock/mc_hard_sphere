module variables
 implicit none
 integer   :: maxnat
 parameter (maxnat=5000)
 double precision :: rx(maxnat), ry(maxnat), rz(maxnat)

 common /bck00/ rx, ry, rz

 logical   :: lattice, structure

 common /bck01/  lattice, structure

 double precision :: boxx, boxy, boxz, vol,boxxinv, boxyinv, boxzinv
 double precision :: hboxx, hboxy, hboxz,eps, eps4, eps48, sig, sig2, mass, pi
 double precision :: rcut, rcut2, temp, dens, dr, succ, beta, wtest
 integer :: nat, ncyclos, ncelda, nmoves, nprint, nfile,nprome, najusta, nghost, nwidom

 common /bck02/  boxx,boxy,boxz,vol,boxxinv,boxyinv,boxzinv,hboxx,hboxy,hboxz
 common /bck02/  eps,eps4,eps48,sig,sig2,mass,pi,rcut,rcut2,temp,dens,dr,succ,beta,wtest 
 common /bck02/  nat,ncyclos,ncelda,nmoves,nprint,nfile,nprome,najusta,nghost,nwidom

 double precision :: ac_en, ac_en2, ac_vr, ac_vr2, ac_pr, ac_pr2
 double precision :: ave_enp, ave2_enp, ave_press, ave2_press
 integer  :: nacc, iseed, intentos, conta_sample

 common /bck03/  ac_en,ac_en2,ac_vr,ac_vr2,ac_pr,ac_pr2,ave_enp,ave2_enp,ave_press,ave2_press
 common /bck03/  nacc,iseed,intentos,conta_sample

 integer :: max_bin
 double precision  :: delta_r
 parameter (delta_r =0.050, max_bin =5000)
 integer :: nx_bin(max_bin), ny_bin(max_bin),nz_bin(max_bin), n_bin, con_pro
 double precision :: rhoxt(max_bin), rhoyt(max_bin),rhozt(max_bin), max_box, volx, voly, volz 

 common /bck04/  rhoxt,rhoyt,rhozt,max_box,nx_bin,ny_bin,nz_bin,volx,voly,volz,n_bin,con_pro
endmodule