************************************************************************

      parameter ( nn=10000,ntp=20000,nx=410,ny=410,nz=101 )

************************************************************************

      implicit real*8(a-h,o-z)
      
      common/m/      grav,amax,epsh

      common/reals/  prate, tsat_min

      common/reals/  flux, tp, start, finish

      common/reals/  tmax,dt,t,tc,t_rain,beta, 
     &               xk, dz, dtinfl, 
     &               r8fluxin, r8fluxout, r8error, r8newmass

      common/reals/  dxdum, r8mV,r8etaV, r8alphaV,
     &               r8mB, r8etaB, r8alphaB

      common/array/  xc(0:nx,0:ny), yc(0:nx,0:ny), zc(0:nx,0:ny), 
     &               dzz(0:nx,0:ny,2), xnc(0:nx,0:ny), 
     &               vegc(0:nx,0:ny)

      common/array/  tt(0:ntp)
      
      common/array/  f(0:nx,0:ny,3,2), q(0:nx,0:ny,3), fdum(3), qs(3),
     &               allf(0:nx,0:ny), xflux0(0:nx,0:ny),
     &               yflux0(0:nx,0:ny),xflux1(0:nx,0:ny),
     &               yflux1(0:nx,0:ny)      
    
      common/array/  dh(0:nx,0:ny,2), du(0:nx,0:ny,2), dv(0:nx,0:ny,2)
      
      common/array/  itype(nx,ny,4), ipos(nx,ny,4), inum(nx,ny),
     &               nop(nx,ny,4), kbeg(nx), kend(nx)

      common/array/  sn(0:nx,0:ny,2), cn(0:nx,0:ny,2), 
     &               ds(0:nx,0:ny,2), dxi(0:nx,0:ny,2),
     &               deta(0:nx,0:ny,2), area(0:nx,0:ny), 
     &               zinflmap2(0:nx, 0:ny),
     &               Tfluxin(0:nx, 0:ny),Tfluxout(0:nx, 0:ny),
     &               Tnewmass(0:nx, 0:ny),Terror(0:nx, 0:ny),
     &               Sfluxin(0:nx, 0:ny),Sfluxout(0:nx, 0:ny),
     &               Snewmass(0:nx, 0:ny),Serror(0:nx, 0:ny)

      common/srctrm/ fix(nx,ny,3),
     &               sx(0:nx,0:ny), sy(0:nx,0:ny), 
     &               t0(0:nx,0:ny)

      common/ints/   ncol, np, ne, nbcell, ilim,
     &               nprt, it, itp, nt, iprt,
     &               nrow, isurf, iBC,
     &               isweed_theta, isweed_flux

      common/ints/   inv,  isetflux, niter, iscouple, iscale,
     &               ifixh, infL, isveg,
     &               jveg,kveg,jbare,kbare,ipass

      common/reals/  htop, stop_tol, stop_tol0

      common/array/  h(0:nx,0:ny), u(0:nx,0:ny), v(0:nx,0:ny),
     &               hn(nn), un(nn), vn(nn),
     &               hp(0:nx,0:ny), up(0:nx,0:ny), vp(0:nx,0:ny),
     &               eta(0:nx,0:ny),  etap(0:nx,0:ny)
    
      common/array/  alpha(nz,2), theta_S(nz,2), theta_R(nz,2),
     &               r8lambda(nz,2), r8n(nz,2), r8m(nz,2),
     &               r8Ksat(nz,2)

      common/array/  DeltaPlus(nz,nz), DeltaMinus(nz,nz),
     &               APlus(nz,nz), AMinus(nz,nz)
     
      common/array/  r8H(0:nx,0:ny,nz), r8THETA(0:nx,0:ny,nz), 
     &               r8K(0:nx,0:ny,nz),  oldTHETA(0:nx,0:ny,nz)
                
      common/array/  z(nz), iterations(ntp)   
      