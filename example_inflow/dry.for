      program swe2d
C======================================================================
C  SVE (Saint Venant) Compontent
C.   time variables:
C         t :  time, ranges from 0 to tmax, updated at start of time step
C         nt :   number of time steps   (ntp>nt)        
C         it   : iteration number  
C         dt   : time increment
C         itp  : print iteration
C         tmax : simulation final time
C         tp  :  time of ponding (printed in console)
C         nprt : print frequency : SVE output is saved every nprt timesteps, 
C               where nprt = dt_p/dt
C    Richards Parameters: 
C       nz:   number of soil layers
C
C     Common variables for putting infiltration on a new time grid:
C         iscale :  scale factor for shallow water to infiltration time steps  (probably 10)
C         dtinfl :  infiltration time step  (dt in python notebook)
C
C       Precipitation variables: 
C         prate(nt)   :  m/s
C         tr   : time when rain stops
C       
C
C     Common SW variables:
C        grav  : acceleration due to gravity.
C        epsh : depth tolerence for dry bed problems.  in m
C        xni, xnv  :  Manning bed roughness coefficients (interspace and veg)
C
C        inum  : number of boundary faces in a boundary cell.
C        ipos  : orientation of boundary face:=1 bottom:=2 right:=3 top:=4 left.
C        nbcell : number of boundary cells.
C   
C    Common variables for infiltration
C        flux :  defined in source subroutine, used by timestep subroutine
C                 if isetflux=1, set flux (in cm/s)
C
C    Non-common infiltration variables:
C        PI  m/s
C        prate(nt)  cm/s
C     
C   Notes:  Richards solver is called in the corrector step.
C           source subroutine sets i (infiltration) for the current and following 
C           iscale timesteps, via qs.
C           Once h --> 0, winflt is set to zero.
C
C     Flux variables:
C         f(j,k,,) : m2/s
C         xflux0, xflux1, yflux0, yflux1:  m3
************************************************************************
      include 'dry.inc'    
      
C		  SWE input files
      open(2,file = 'input/coords.dat')
      open(3,file = 'input/boundary.dat')     
      open(4,file = 'input/params.dat') 
      open(5,file = 'input/veg.dat')
      open(11,file = 'input/nodes.dat')        
C     Richards input files
      open(304,file='input/vanG.dat')
C     Output files: saved every nprt time steps (dt_p seconds)
      open(101,file= 'output/time.out')                           
      open(100,file= 'output/h.out')                
      open(102,file= 'summary.txt')                  
      open(103,file= 'output/cfl.out')   
      open(104,file= 'output/hydro.out')  
      open(107,file= 'output/ptsTheta.out')
      open(110,file= 'output/fluxes1234.out')   ! boundary fluxes      
      open(111,file= 'output/dvol.out')  ! SVE volume tracking              

C 	  Richards output files
      call cpu_time(start)       
C      Read input, set up grid, initial conditions.
      call init
      call initvg


C     write headings for fortran output files
      write(101,*)  "     time (s) |   CFL  "
      write(104,*)  "     time (s) |  hydrograph (m3/s)" 
C     heading for volume tracking file (dvol.out)
      write(111,* )  "time (s)      vol               zflux         ",
     &                "     zinfl               zrain "
C     heading for horizontal flux tracking file (fluxes1234.out)
      write(110,* )  "      flux1               flux2        ",
     &                "    flux3               flux4 "

      write(100, *) "   j   k    h              u               v ",
     &    "           zinflmap2        xflux0         yflux0",
     &    "         xflux1          yflux1"

C     Initialize variables
      flux1 = 0.d0
      flux2 = 0.d0
      flux3 = 0.d0
      flux4 = 0.d0
      
      xflux0 = 0.0
      yflux0 = 0.0
      xflux1 = 0.0
      yflux1 = 0.0

C     tracks infiltration volume
      zinfl = 0.d0  
      zinflmap2 = 0.d0

      zrain = 0.d0
      dvol = 0.d0
      zflux = 0.d0

      vol = 0.D0
      do j=1,ncol 
        do k=kbeg(j),kend(j)
          vol = vol + h(j,k)*area(j,k)
        enddo
      enddo 
      
C     Write ICs to output files
      write(110,202) flux1, flux2, flux3, flux4  
      write(111,200) t, dvol, zflux, zinfl, zrain
      call myoutput

C    Begin time loop.
      do it = 0,nt-1
        t = t + dt
        if (t .gt. tr) prate = 0.		
        amax = 0.d0
    
C       Compute predictor.
        do j=1,ncol
          do k=kbeg(j),kend(j)
            call bconds(j,k,h,u,v)                              
            call predict(j,k)
          enddo
        enddo
        
C Loop over cells to compute fluxes.
        do j=1,ncol
        do k=kbeg(j),kend(j)
          call bconds(j,k,hp,up,vp)
          call fluxes(j-1,j,k,k,1)          ! vertical faces.   
          call fluxes(j,j,k-1,k,2)          ! horizontal faces.
          
          do l=1,inum(j,k)       
            if(ipos(j,k,l) .eq. 3) then
              call fluxes(j,j,k,k+1,2)       ! top boundaries. 
              flux3 =  flux3 - f(j,k+1,1,2)*ds(j,k+1,2)*dt
            
            elseif(ipos(j,k,l) .eq. 2) then
              call fluxes(j,j+1,k,k,1)              ! right boundaries.
              flux2 = flux2  - f(j+1,k,1,1)*ds(j+1,k,1)*dt
            
            elseif(ipos(j,k,l) .eq. 1) then         ! lower boundary
              flux1 = flux1 + f(j,k,1,2)*ds(j,k,2)*dt
            
            elseif(ipos(j,k,l) .eq. 4) then         ! left boundary
              flux4 = flux4 + f(j,k,1,1)*ds(j,k,1)*dt
            endif
          enddo

C         Increment the fluxes at each grid cell  
          xflux0(j,k) = xflux0(j,k) + f(j,k,1,1)*ds(j,k,1)*dt        
          xflux1(j,k) = xflux1(j,k) + f(j+1,k,1,1)*ds(j+1,k,1)*dt   
          yflux0(j,k) = yflux0(j,k) + f(j,k,1,2)*ds(j,k,2)*dt   
          yflux1(j,k) = yflux1(j,k) + f(j,k+1,1,2)*ds(j,k+1,2)*dt       
        enddo
        enddo
      

C Compute corrector solution.
        do j=1,ncol
        do k=kbeg(j),kend(j)
          call source(j,k,hp(j,k),up(j,k),vp(j,k), 1)
C         increment source (p-i) volume   
          zinfl = zinfl + dt*qs(1)*area(j,k)  
C         increment volume of rain inputs (as an added mass balance check)
          zrain = zrain + dt*prate*area(j,k)   
C         compute infiltration for this grid cell and timestep. 
          dzinfl =  dt*qs(1)*area(j,k)- dt*prate*area(j,k) 
          zinflmap2(j,k) = zinflmap2(j,k) + dzinfl  !  units : m^3 

         do l=1,3
            q(j,k,l) = q(j,k,l) + dt/area(j,k)*(
     &        f(j,k,l,2)*ds(j,k,2) + f(j,k,l,1)*ds(j,k,1) -  
     &        f(j+1,k,l,1)*ds(j+1,k,1) - f(j,k+1,l,2)*ds(j,k+1,2)) 
     &        + dt*qs(l)
          enddo
       
        enddo
        enddo
                
C Store solution. 
C Solve continuity equation even in all cells, but 
C solve momentum equations in wet cells only.
        do j=1,ncol
          do k=kbeg(j),kend(j)
C Check for negative depth.
            if(q(j,k,1) .ge. 0.D0) then
              h(j,k) = q(j,k,1)
            else
              zinfl = zinfl - q(j,k,1)*area(j,k)
              zinflmap2(j,k) = zinflmap2(j,k) - q(j,k,1)*area(j,k)   							                      
              q(j,k,1) = 0.D0
              h(j,k) = 0.D0
            endif
C Neglect momentum in nearly dry cells.
           if(h(j,k) .lt. epsh) then  
           	  u(j,k) = 0.d0
              v(j,k) = 0.d0
C               t0(j,k) = t
              do l=2,3
                q(j,k,l) = 0.d0
              enddo
            elseif (h(j,k) .ge. epsh) then  
              u(j,k) = q(j,k,2)/h(j,k)
              v(j,k) = q(j,k,3)/h(j,k)
            endif
           enddo
       enddo 

C      Below here only executed every nprt time steps        
        if (mod(it, nprt) .eq. nprt-1)  then         
          vol0 = vol
          vol = 0.D0
          do j=1,ncol 
            do k=kbeg(j),kend(j)
              vol = vol + h(j,k)*area(j,k)
            enddo
          enddo 
        
          dvol = vol - vol0
          
          zinfl = zinfl - zrain   
C         sum of all lateral fluxes at boundaries
          zflux =  flux1 + flux2 + flux3 + flux4   
C        fluxes are positive  out of the domain
          write(110,202) flux1, flux2, flux3, flux4  
          
          write(111,200) t, dvol, zflux, zinfl, zrain
                      
C         update print increment  
          itp = itp + 1     
          call myoutput
C           write(*,203) t, amax*dt, prate*3600*100, ! cm/hr
C      &              maxval(h)*100  ! cm

          flux1 = 0.d0
          flux2 = 0.d0
          flux3 = 0.d0
          flux4 = 0.d0    
          
C         allf = 0.0
          xflux0 = 0.0     
          yflux0 = 0.0  
          xflux1 = 0.0     
          yflux1 = 0.0               
          zinfl = 0.d0
          zrain = 0.d0
     
     
        endif
        
        nprthydro = 1/dt 
        if (mod(it, nprthydro) .eq. nprthydro-1) then
          call myhydro
        endif 
         
C        minimum volume to
         r8minvol = 0.001d0*epsh*dxdum**2*nrow*ncol

         if ((vol .le. r8minvol) .and. (t .gt. tr)) then

          call gracefulExit  
          write(102,*) 'No more water, t=', t/60
          return
        endif

        if(amax*dt.gt. 10.0d0) then
          call gracefulExit
          write(*,*) 'AMAX too big'
          write(102,*) 'AMAX too big'
     
          return   ! leave early          
        endif
        if(t .gt. tmax) then
          call gracefulExit
          write(102,*) 'out of time'
     
          return   ! leave early          
        endif

        
      enddo
      call gracefulExit 
      
      stop
      
 200  format(' ', f8.1, 5e19.8)  ! for writing mass components to dvol.out
 211  format(' ', 5A12 ) 
 201  format(' ', A14, f12.2) 
 202  format(' ', 5e19.8)  ! for writing fluxes to allfluxes.out
 203  format(' ','time =', f8.1,'s; max CFL =',f7.3,
     &       ' ; p =',f5.1, 'cm/hr', '; depth =', f8.2)

      end

************************************************************************
      subroutine gracefulExit
      include 'dry.inc'
      
      call cpu_time(finish)
!       print '("Time = ",f10.3," seconds.")',finish-start
      write(102,*) '  ' 
      write(102,*) 'Results summary:'
      write(102,211) 'runtime = ', (finish-start)
      write(102,211) '1st ponding  =  ', tp
      write(102,211) 'final time = ', t
      

      return 
      
 211  format(' ', A14, f12.2) 
 212  format(' ', A14, f12.3, A2)  
      end

************************************************************************
      subroutine myhydro
      include 'dry.inc'
      
      hydro = 0.d0
      do j=1,ncol
            hydro = hydro + f(j,1,1, 2)*ds(j,1,2)  ! m3/s
      enddo
      write(104, 204)  t, hydro      
      
      return
 202  format(' ', i8, f9.2)
!  203  format(' ', f13.2 , 4i10, 4i10 )
 204  format(' ', f15.2 , E20.10) 
      end
************************************************************************
      subroutine myoutput
      include 'dry.inc'
  
C     file 101 is 'output/time.out' - to keep track of the time stpes
      write(101, *) t, amax*dt 
      
C     file 100 is 'output/h.out' 
      do j=1,ncol
        do k=kbeg(j),kend(j)
          write(100, 201) j, k, h(j,k), u(j,k), v(j,k),
     &                  zinflmap2(j,k), xflux0(j,k), yflux0(j,k),
     &                  xflux1(j,k), yflux1(j,k)
          zinflmap2(j,k) = 0.d0       
        enddo
      enddo
	    write(100, 202) itp, t
      
C    Save H and Theta for two points - one vegetated and one bare
      do l = 1,nz
          write(107, 203) l,z(l),r8THETA(jveg,kveg,l), 
     &                  r8H(jveg,kveg,l), 
     &                  r8THETA(jbare,kbare,l), 
     &                  r8H(jbare,kbare,l)
      enddo  
      write(107, 202) itp, t

      return
 202  format(' ', i8, f9.2)
 201  format(' ', 2i4, 8e15.6) 
 203  format(' ', i4, f8.3, 4e15.6) 
      end
************************************************************************
      subroutine source(j,k,depth,udum,vdum,lstep)
      !   
      !   Input:
      !         j, k : grid cell
      !         depth  (h(j,k), hp(j,k))
      !         udum  (u(j,k), up(j,k))
      !         vdum  (v(j,k), vp(j,k)) 
      !   Modifies common variables: 
      !         qs(3):
      !              qs(1) = winflt = p-i  
      !         isetflux, htop, flux  
      !         r8THETA, r8H, r8K  (soil solution matrices)
      !
      !  Subroutines called from this subroutine:
      !         potential : 
      !             returns PI (potential infiltration)
      !         timestep:  
      !             returns hp1dum,thetap1dum, r8kp1dum

      ! isetflux :  1 if flux = 0 (no rain) or flux = prate 
      !             0 if ponding  (flux determined by potential infiltration)
      
      include 'dry.inc'      
      
      real ( kind = 8 )   PI, depth, fm, feta, falpha, htop
      real ( kind = 8 )   hdum(nz), thetadum(nz), r8kdum(nz)
      real ( kind = 8 )   hp1dum(nz), thetap1dum(nz), r8kp1dum(nz)  
      
      hdum = r8H(j,k,1:nz)
      thetadum = r8THETA(j,k,1:nz)   
      r8kdum = r8K(j,k,1:nz) 

      zold = depth*100.d0   ! depth in cm
      znew = zold           ! depth in cm
      PI = 0.d0
      
      iskip = 0

      ! determine whether the cell is vegetated or not
			! define the roughness parameters (fm, falpha, feta)
      if (vegc(j,k) .gt. 0) then
        isveg = 1
        fm = r8mV
        falpha = r8alphaV
        feta = r8etaV
      else
        isveg = 2 
        fm = r8mB
        falpha = r8alphaB 
        feta = r8etaB   
      endif  
      
      if (lstep .eq. 1)  then  ! only update depth in corrector step.
        
C       start Richards cases
C 	    Case 1: use infiltration rate from previous timestep if there's ponding, 
C               otherwise assume infiltration rate equals rainfall rate.

        if  (r8kdum(nz) .eq. 0) then   ! if Ks = 0, skip richards
           iskip = 1
           winflt = prate
        elseif (mod(it, iscale) .ne. 0) then  ! not solving Richards
          
          iskip = 1  
C         keep old infiltration rate if depth > 0            
          if (zold .gt. 0.d0) then    
                r8kt =  r8kdum(nz-1)
                znew = zold + prate*dt*100.d0 - r8kt*( 
     &                    (hdum(nz) - hdum(nz-1))/dz + 1.d0)*dt
                winflt = (znew - zold)/dt/100.d0
          else 

                winflt =  0.d0
          endif

C       Case 2: rain and no ponding.  
        elseif ((zold .le. 1e-7) .and. (prate .gt. 0.d0)) then  ! rain but no ponding

C         check whether rain (prate) exceeds the potential infiltration rate 
          call potential(hdum, thetadum, PI)  ! potential subroutine returns PI
C         If PI < rainfall rate, ponding starts
          if (PI .lt. prate*100.d0) then   
                        
            if (tp .gt. 1e5) then
                tp = t
             endif

            isetflux = 0
            htop = zold
            call timestep(hdum,thetadum,hp1dum,thetap1dum,r8kp1dum) 
            znew = zold + prate*dt*100.d0 - PI*dt

          else    ! no ponding,  flux = prate
                                  
            isetflux = 1  
            flux = - prate*100.d0 
            call timestep(hdum,thetadum,hp1dum,thetap1dum,r8kp1dum)

          endif  

          winflt = (znew - zold)/dt/100.d0
C       Case 3: no ponding and no rain
        elseif ((zold.le.1e-7).and.(prate .le. 1e-10)) then

            isetflux = 1
            flux = 0.d0
            htop = zold   !  not needed
            call timestep(hdum,thetadum,hp1dum,thetap1dum, r8kp1dum) 
            winflt = 0.d0      

C       Case 4:  ponding  (rain or no rain)      
        elseif (zold .gt. 0.d0) then       
            
            isetflux = 0
            htop = 0.0
            call timestep(hdum,thetadum,hp1dum,thetap1dum,r8kp1dum)   
            
            r8kt = r8kp1dum(nz-1)  ! update infiltation rate
            znew = zold + prate*dt*100.d0 - r8kt*( 
     &               (hp1dum(nz) - hp1dum(nz-1))/dz + 1.d0)*dt
            winflt =  (znew - zold)/dt/100.d0
        endif 
    

C       End of Richards solver cases 
C       If Richards equation was called  (cases 2-4, iskip = 0), update the solution matrices:
        if (iskip .eq. 0) then  
          
          r8THETA(j,k,1:nz) = thetap1dum  
          r8H(j,k,1:nz) = hp1dum
          r8K(j,k,1:nz) = r8kp1dum

C         Compute the surface flux, drainage and change in mass 
C           (this is not used for anything)
          r8kt = r8kp1dum(nz-1)
          r8fluxin = r8kt*
     &                ((hp1dum(nz) - hp1dum(nz-1))/dz+1)*dt*iscale
     
          if (ipass .eq. 0) then
            r8fluxout = - 0.5*(r8kp1dum(1)+r8kp1dum(2))*
     &                    ((hp1dum(2) - hp1dum(1))/dz+1)*dt*iscale
          else 
            r8fluxout = - r8fluxin         
          endif
     
          r8newmass = sum(r8THETA(j,k,1:nz) - oldTHETA(j,k,1:nz))*dz   
          
          oldTHETA(j,k,1:nz) = r8THETA(j,k,1:nz)    

        endif
        
        else   ! in the predictor step, there is no infiltration or precipitation
          winflt = 0.d0 
        
      endif

      if (depth .gt. epsh) then 
        
        vmag = dsqrt(udum*udum + vdum*vdum) 

        if ((feta .eq. 0.5) .and. (fm .eq. 0.d0)) then  ! cylinder

          fricSx = (falpha)**(2.d0)*udum*vmag 
          fricSy = (falpha)**(2.d0)*vdum*vmag 

	      elseif (feta .eq. 0.5) then  ! most schemes
          
	         fricSx = (falpha/depth**fm)**(2.d0)*udum*vmag 
	         fricSy = (falpha/depth**fm)**(2.d0)*vdum*vmag 

C         elseif (fm .eq. 0.0) then  ! the cylinder array is weird
C
C           fricSx = falpha**(1.57d0)*vmag**0.57d0*udum
C           fricSy = falpha**(1.57d0)*vmag**0.57d0*vdum
C
        elseif (feta .eq. 1) then  ! special case for poisseuille 
        
          fricSx = falpha*udum/depth**fm 
          fricSy = falpha*vdum/depth**fm 

          Rel = vmag*depth/1.e-6   

          if (Rel .gt. 500) then 
              ffact = 0.5  ! okay for smooth surfaces
              fricSx = ffact*udum*vmag/8./grav/depth
              fricSy = ffact*vdum*vmag/8./grav/depth      
          endif

        endif

        qs(1) = winflt

        qs(2) = 0.5D0*udum*winflt - grav*depth*fricSx -  
     &          grav*depth*sx(j,k)
        qs(3) = 0.5D0*vdum*winflt - grav*depth*fricSy - 
     &          grav*depth*sy(j,k)
       else
          
        qs(1) = winflt
        qs(2) = 0.d0
        qs(3) = 0.d0
       endif
  
      return

!  208  format(' ', ' it =',i5,'; t =', f8.1
!      &    ' ; depth =', f10.6, ' ', A14 )
 209  format(' ', ' time of ponding is =',f8.1,'s ')           
!  207   format(' ','it =',i8,'; niter =',i7,
!      &    ' ; prate =',f7.1, 'cm/hr',
!      &    ' ; infl =',f7.1, 'cm/hr',
!      &    ' ; depth =', f10.6, 'cm' ,
!      &    ' ; winflt =', f10.6 )

      end
      
    
************************************************************************
      subroutine fluxes(jl,jr,kl,kr,i1)
      include 'dry.inc'
      !   
      !   Input:
      !             jl, jr, kl, kr : left/right grid cell indices
      !             i1 :  interface type (horizontal or vertical)
      
      !   Modifies common variables:
      !             f(0:nx,0:ny,3,2)
      !             fdum(3) 
                  
C MUSCL extrapolation at cell interface.
      hl = hp(jl,kl) + 0.5D0*dh(jl,kl,i1)
      ul = up(jl,kl) + 0.5D0*du(jl,kl,i1)
      vl = vp(jl,kl) + 0.5D0*dv(jl,kl,i1)
      hr = hp(jr,kr) - 0.5D0*dh(jr,kr,i1)
      ur = up(jr,kr) - 0.5D0*du(jr,kr,i1)
      vr = vp(jr,kr) - 0.5D0*dv(jr,kr,i1)
      snn = sn(jr,kr,i1)
      cnn = cn(jr,kr,i1)
      
	    if(i1 .eq. 1) then
        dx =  deta(jr,kr,2)*area(jr,kr) 
	      dy = - deta(jr,kr,1)*area(jr,kr)
      else
        dx = - dxi(jr,kr,2)*area(jr,kr) 
	      dy =   dxi(jr,kr,1)*area(jr,kr)
      endif
C Needed for dry bed problems.
	    if(hl .lt. 0.D0) hl = 0.D0
      if(hr .lt. 0.D0) hr = 0.D0
C Compute arithmatic averages for source terms.
      havg = 0.5D0*(hl + hr)
      uavg = 0.5D0*(ul + ur)
      vavg = 0.5D0*(vl + vr)
C Prevent leakage into cells with higher bed elevation.
      etal = hp(jl,kl) + zc(jl,kl)
      etar = hp(jr,kr) + zc(jr,kr)
C Fluxes and source terms.
      if(havg .le. 0.D0 ) then
	      do i=1,3  
          f(jr,kr,i,i1) = 0.D0        
        enddo
      else   
        call solver(hl,hr,ul,ur,vl,vr,fdum,snn,cnn,dx,dy)
	      do i=1,3
          f(jr,kr,i,i1) = fdum(i)
        enddo
      endif

	    return
      end 
      
************************************************************************
      subroutine solver(hl,hr,ul,ur,vl,vr,ff,sndum,cndum,dx,dy)        
      !
      !   Input:
      !             hl,hr,ul,ur,vl,vr : left/right h,u,v
      !             i1 :  face numbers
    
      !   Output: 
      !             f(3)  --> fdum in fluxes

      include 'dry.inc'
      !implicit real*8(a-h,o-z)
      dimension ws(3), e(3,3), a(3), astar(3), da(3), ff(3), dum(3)
      !common/m/ grav, amax, epsh

C Compute Roe averages at cell face.
      duml  = dsqrt(hl)
      dumr  = dsqrt(hr)
      hhat  = duml*dumr
      uhat  = (duml*ul + dumr*ur)/(duml + dumr)
      vhat  = (duml*vl + dumr*vr)/(duml + dumr)
      chat  = dsqrt(0.5D0*grav*(hl + hr))
      uperp = uhat*cndum + vhat*sndum
C Compute eigenvalues.  Lambdahat
      a(1) = uperp - chat
      a(2) = uperp
      a(3) = uperp + chat
C Compute approximate wave strengths.
      dhdum   = hr - hl
      dudum   = ur - ul
      dvdum   = vr - vl
      dupar = -dudum*sndum + dvdum*cndum
      duperp=  dudum*cndum + dvdum*sndum
      !  deltaVhat (eqn 31 in 1999)
      ws(1) = 0.5D0*(dhdum - hhat*duperp/chat)
      ws(2) = hhat*dupar     
      ws(3) = 0.5D0*(dhdum + hhat*duperp/chat)
C Compute right eigenvectors.  Rhat
      e(1,1) = 1.D0
      e(2,1) = uhat - chat*cndum   
      e(3,1) = vhat - chat*sndum
      e(1,2) = 0.D0
      e(2,2) = -sndum
      e(3,2) =  cndum
      e(1,3) = 1.D0
      e(2,3) = uhat + chat*cndum
      e(3,3) = vhat + chat*sndum
C Entropy fix.
	    dl = dsqrt(dx*dx + dy*dy)
      cl = dsqrt(grav*hl)
      cr = dsqrt(grav*hr)
      uperpl = ul*cndum + vl*sndum
      uperpr = ur*cndum + vr*sndum
      al1 = uperpl - cl
      al3 = uperpl + cl
      ar1 = uperpr - cr
      ar3 = uperpr + cr
      da(1) = dmax1(0.D0, 4.D0*(ar1 - al1))
      da(2) = 0.d0
      da(3) = dmax1(0.D0, 4.D0*(ar3 - al3))
      do i=1,3   
         if(dabs(a(i)) .lt. 0.5D0*da(i)) then
           astar(i) = a(i)*a(i)/da(i) + 0.25D0*da(i)  
          else
           astar(i) = dabs(a(i))
         endif
         if(astar(i)/dl .gt. amax) amax = astar(i)/dl
      enddo
C Compute flux increments.
      do i=1,3
        dum(i) = 0.D0
        do l=1,3
          dum(i) = dum(i) + (astar(l)*ws(l))*e(i,l) 
        enddo
      enddo
C Add flux to appropriate cell faces.
      ff(1) = 0.5D0*(f1(hl,uperpl) + f1(hr,uperpr) - dum(1))
      ff(2) = 0.5D0*(f2(hl,ul,uperpl,cndum) 
     &            + f2(hr,ur,uperpr,cndum) - dum(2))
      ff(3) = 0.5D0*(f3(hl,vl,uperpl,sndum)   
     &            + f3(hr,vr,uperpr,sndum) - dum(3))

      return
      end
      
************************************************************************
      subroutine predict(j,k)
      include 'dry.inc'
      
      !   Predictor step 
      !   Input:
      !             j, k : grid cell
      
      !   Modifies common variables:
      !             dh, du, dv  
      !             qs(3)  -  source term
      !             hp, up, vp  - predictor variables

      do kk=1,2                  ! loop over coord. directons.
        if(kk .eq. 1) then
          jr = j + 1
          jl = j - 1
          kr = k
          kl = k
        else
          jr = j
          jl = j
          kr = k + 1
          kl = k - 1
        endif
C Compute gradients only in wet cells.
	      if(h(j,k) .ge. epsh) then
C Limit afree surface elevation to reduce dissipation..
          dh1 = h(j,k) + zc(j,k) - h(jl,kl) - zc(jl,kl)
          dh2 = h(jr,kr) + zc(jr,kr) - h(j,k) - zc(j,k)
C Needed to minimize dissipation at wet/dry interfaces.
          if(h(jl,kl) .lt. epsh) dh1 = 2.d0*dh1 
          if(h(jr,kr) .lt. epsh) dh2 = 2.d0*dh2
          call limitr(ilim,beta,dh1,dh2,dhh)
          dh(j,k,kk) = dhh - dzz(j,k,kk)
C U velocity.
          du1 = u(j,k) - u(jl,kl)
          du2 = u(jr,kr) - u(j,k)
          call limitr(ilim,beta,du1,du2,duu)
          du(j,k,kk) = duu
C V velocity.
          dv1 = v(j,k) - v(jl,kl)
          dv2 = v(jr,kr) - v(j,k)
          call limitr(ilim,beta,dv1,dv2,dvv) 
          dv(j,k,kk) = dvv
	      else
          dh(j,k,kk) = 0.D0
          du(j,k,kk) = 0.D0
          dv(j,k,kk) = 0.D0
        endif
      enddo
C Generalized velocities.
	    uxi = u(j,k)*dxi(j,k,1) + v(j,k)*dxi(j,k,2)
	    ueta = u(j,k)*deta(j,k,1) + v(j,k)*deta(j,k,2)
C Predictor.
! 		Call source before updating the predictor valuess
      call source(j, k, h(j,k), u(j,k), v(j,k), 0)
      if(h(j,k) .ge. epsh**(0.75d0)) then           
        qs(2) = qs(2)/h(j,k)
        qs(3) = qs(3)/h(j,k)
      else
        qs(2) = 0.d0
        qs(3) = 0.d0
	    endif
      hp(j,k) = h(j,k) - 0.5D0*dt*(
     &    uxi*dh(j,k,1) + h(j,k)*(dxi(j,k,1)*du(j,k,1) + 
     &                            dxi(j,k,2)*dv(j,k,1)) +   
     &    ueta*dh(j,k,2) + h(j,k)*(deta(j,k,1)*du(j,k,2) +
     &                             deta(j,k,2)*dv(j,k,2)) + qs(1))   
      up(j,k) = u(j,k) - 0.5D0*dt*(      
     &    grav*dxi(j,k,1)*dh(j,k,1) + uxi*du(j,k,1) +
     &    grav*deta(j,k,1)*dh(j,k,2) + ueta*du(j,k,2) + qs(2))
      vp(j,k) = v(j,k) - 0.5D0*dt*(
     &    grav*dxi(j,k,2)*dh(j,k,1) + uxi*dv(j,k,1) +
     &    grav*deta(j,k,2)*dh(j,k,2) + ueta*dv(j,k,2) + qs(3))
C Correct any negative depths.
	    if(hp(j,k) .lt. 0.d0) then
        hp(j,k) = 0.d0
        dh(j,k,1) = 0.d0
        dh(j,k,2) = 0.d0
      endif
C Neglect momentum in nearly dry cells.
	    if(hp(j,k) .le. epsh) then
        up(j,k) = 0.D0
        vp(j,k) = 0.D0
        do i=1,2
          du(j,k,i) = 0.d0
          dv(j,k,i) = 0.d0
        enddo
      endif

      return
      end
      
************************************************************************
      subroutine bconds(j,k,hdum,udum,vdum)
      include 'dry.inc'

      dimension hdum(0:nx,0:ny), udum(0:nx,0:ny), vdum(0:nx,0:ny)
C     Loop over all boundary faces in the cell. 
      
      do i=1, inum(j,k)
        if(ipos(j,k,i) .eq. 1) then ! front face.         
            jj = j
            kk = k-1
            jl = j
            kl = k 
            j2 = j
            k2 = k+1
            io = 2 
        elseif(ipos(j,k,i) .eq. 2) then ! right face.
            jj = j+1
            kk = k 
            jl = j+1 
            kl = k  
            j2 = j-1
            k2 = k
            io = 1
        elseif(ipos(j,k,i) .eq. 3) then ! back face.
            jj = j
            kk = k+1
            jl = j
            kl = k+1
            j2 = j
            k2 = k-1 
            io = 2 
        elseif(ipos(j,k,i) .eq. 4) then ! left face.
            jj = j-1
            kk = k 
            jl = j
            kl = k 
	          j2 = j+1
            k2 = k  
            io = 1  
        endif
	      t0(jj,kk) = t0(j,k)
C   Open boundary.
        if(itype(j,k,i) .eq. 0)then
          
          dh(jj,kk,io) = dh(j,k,io)
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)
          
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
          udum(jj,kk) = 2.D0*udum(j,k) - udum(j2,k2)
          vdum(jj,kk) = 2.D0*vdum(j,k) - vdum(j2,k2)
                    
          if (hdum(j,k) .gt. hdum(j2,k2) + 1e-7) then
          
            if (iBC .eq. 0) then  ! normal!
              hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
              
            elseif (iBC .eq. 1) then  ! same depth (not same height)
              hdum(jj,kk) = hdum(j,k)
            
            elseif (iBC .eq. 2) then   ! mirror
              hdum(jj,kk) = hdum(j2,k2) 

            endif 

          endif
               
          
C Wall boundary.
       elseif(itype(j,k,i) .eq. 1)then
          dh(jj,kk,io) = dh(j,k,io)
          du(jj,kk,io) = du(j,k,io)
          dv(jj,kk,io) = dv(j,k,io)
          
          hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)          

          udum(jj,kk) = udum(j,k)*(sn(jl,kl,io)*sn(jl,kl,io) - 
     &                                 cn(jl,kl,io)*cn(jl,kl,io)) -
     &                   2.D0*vdum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
          vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,io)*cn(jl,kl,io) - 
     &                      sn(jl,kl,io)*sn(jl,kl,io)) -
     &                      2.D0*udum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
        
C Specified depth and velocity (supercritical).
       elseif(itype(j,k,i) .eq. 2)then
          
          dh(jj,kk,io) = 0.D0
          du(jj,kk,io) = 0.D0
          dv(jj,kk,io) = 0.D0                 
          hdum(jj,kk) = fix(j,k,1) 
          udum(jj,kk) = fix(j,k,2)
          vdum(jj,kk) = fix(j,k,3)
          if(isurf .eq. 1) hdum(jj,kk) = hdum(jj,kk) - zc(j,k)
      
C Specified flow rate (subcritical).
       elseif((itype(j,k,i) .eq. 4) .and. (t .lt. tr)) then 
          du(jj,kk,io) =  0.d0
          dv(jj,kk,io) =  0.d0 
          if(hdum(j,k) .ge. epsh) then
            hdum(jj,kk) = hdum(j,k) 
            dh(jj,kk,io) =  0.d0
           elseif(hdum(j,k)*hdum(j2,k2) .ge. epsh) then
            hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)
            dh(jj,kk,io) = dh(j,k,io)
          endif
	        if(hdum(jj,kk) .ge. epsh) then
            udum(jj,kk) = fix(j,k,2)/hdum(jj,kk)
            vdum(jj,kk) = fix(j,k,3)/hdum(jj,kk)
          else
            udum(jj,kk) = 0.D0
            vdum(jj,kk) = 0.D0
          endif   
C Prescribe wall boundary after tr     
       elseif ((itype(j,k,i) .eq. 4) .and. (t .ge. tr)) then 
         dh(jj,kk,io) = dh(j,k,io)
         du(jj,kk,io) = du(j,k,io)
         dv(jj,kk,io) = dv(j,k,io)
         
         hdum(jj,kk) = 2.D0*hdum(j,k) - hdum(j2,k2)          

         udum(jj,kk) = udum(j,k)*(sn(jl,kl,io)*sn(jl,kl,io) - 
     &                                 cn(jl,kl,io)*cn(jl,kl,io)) -
     &                   2.D0*vdum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
         vdum(jj,kk) = vdum(j,k)*(cn(jl,kl,io)*cn(jl,kl,io) - 
     &                      sn(jl,kl,io)*sn(jl,kl,io)) -
     &                      2.D0*udum(j,k)*sn(jl,kl,io)*cn(jl,kl,io)
       
       endif


        if (hdum(jj,kk)  .lt. 0.D0) hdum(jj,kk) = 0
        
        if(hdum(jj,kk) .lt. epsh) then
	        udum(jj,kk) = 0.D0
          vdum(jj,kk) = 0.D0
	      endif
        
      enddo

      return
      end
      
***********************************************************************
      subroutine grid
      include 'dry.inc'
      !     Read grid data from  'coords.dat', 'veg.dat' and 'nodes.dat'
      !       file 2 = 'input/coords.dat'
      !       file 5 = 'input/veg.dat'
      !       file 11 = 'input/nodes.dat'

      real ( kind = 8 ) x(nn), y(nn), zz(nn)
      real ( kind = 8 ) veg(nn)

			area = 0.d0
      read(2,*) np, ne
    	if(np .gt. nn) then
    	  write(*,*) np
        write(*,*) 'ERROR:  parameter nn in file dry.inc is too small'
    	  stop
    	endif
      do i=1,np
        read(2,*) x(i), y(i), zz(i)
      enddo
      do i=1,np
        read(5,*) veg(i)
      enddo 
      do j=1,ncol
        do k=kbeg(j),kend(j)
          read(11,*) (nop(j,k,i), i=1,4)
        enddo
      enddo
C Compute grid metrics.
      do j=1,ncol
      do k=kbeg(j),kend(j)
        n1 = nop(j,k,1)
        n2 = nop(j,k,2)
        n3 = nop(j,k,3)
        n4 = nop(j,k,4)
        xc(j,k) = 0.25D0*(x(n1) + x(n2) + x(n3) + x(n4))
        yc(j,k) = 0.25D0*(y(n1) + y(n2) + y(n3) + y(n4))
        vegc(j,k) = veg(n1) 
        if (vegc(j,k) .gt. 0) then
          xnc(j,k) = r8alphaV
        else
          xnc(j,k) = r8alphaB
        endif        
        dxdxi = 0.5D0*(-x(n1) + x(n2) + x(n3) - x(n4))
        dxdeta = 0.5D0*(-x(n1) - x(n2) + x(n3) + x(n4))
        dydxi = 0.5D0*(-y(n1) + y(n2) + y(n3) - y(n4))
        dydeta = 0.5D0*(-y(n1) - y(n2) + y(n3) + y(n4))
        area(j,k) = dxdxi*dydeta - dxdeta*dydxi
        if(area(j,k) .le. 0.D0)then
           write(*,*) 'area error in cell ',j,k
           stop
        endif
        dxi(j,k,1) = dydeta/area(j,k)
        deta(j,k,1) = -dydxi/area(j,k)
        dxi(j,k,2) = -dxdeta/area(j,k)
        deta(j,k,2) = dxdxi/area(j,k)
        sx(j,k)=((zz(n2)-zz(n4))*(y(n3)-y(n1))-(zz(n3)-zz(n1))
     &     *(y(n2)-y(n4)))/(2.D0*area(j,k))
        sy(j,k)=((zz(n3)-zz(n1))*(x(n2)-x(n4))-(zz(n2)-zz(n4))
     &     *(x(n3)-x(n1)))/(2.D0*area(j,k))
        zc(j,k) = 0.25D0*(zz(n1) + zz(n2) + zz(n3) + zz(n4))
    	  dzz(j,k,1) = sx(j,k)*dxdxi + sy(j,k)*dydxi
    	  dzz(j,k,2) = sx(j,k)*dxdeta + sy(j,k)*dydeta
      enddo
      enddo
      dxdum = dxdxi    
C Compute cell face angles.
      do j=1,ncol
      do k=kbeg(j),kend(j)
       ddx = x(nop(j,k,2)) - x(nop(j,k,1))
       ddy = y(nop(j,k,2)) - y(nop(j,k,1))
       ds(j,k,2) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,2) = ddx/ds(j,k,2)            ! Horizontal face.
       cn(j,k,2) = -ddy/ds(j,k,2)          ! Horizontal face.
       ddx = x(nop(j,k,4)) - x(nop(j,k,1))
       ddy = y(nop(j,k,4)) - y(nop(j,k,1))
       ds(j,k,1) = dsqrt(ddx*ddx + ddy*ddy)
       sn(j,k,1) = -ddx/ds(j,k,1)            ! Vertical face.
       cn(j,k,1) = ddy/ds(j,k,1)
       do i=1,inum(j,k)
        if(ipos(j,k,i) .eq. 3) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,4))
           ddy = y(nop(j,k,3)) - y(nop(j,k,4))
           ds(j,k+1,2) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j,k+1,2) = ddx/ds(j,k+1,2)     ! Top (boundary) faces.
           cn(j,k+1,2) = -ddy/ds(j,k+1,2)
         elseif(ipos(j,k,i) .eq. 2) then
           ddx = x(nop(j,k,3)) - x(nop(j,k,2))
           ddy = y(nop(j,k,3)) - y(nop(j,k,2))
           ds(j+1,k,1) = dsqrt(ddx*ddx + ddy*ddy)
           sn(j+1,k,1) = -ddx/ds(j+1,k,1)     ! Right (boundary) faces.
           cn(j+1,k,1) = ddy/ds(j+1,k,1)
        endif
       enddo
      enddo
      enddo
C Set some things in ghost cells.
      do j=1,ncol
      do k=kbeg(j),kend(j) 
        do i=1, inum(j,k)
          call findbc(i,j,k,jj,kk,j2,k2)
          area(jj,kk) = area(j,k)
          sx(jj,kk) = sx(j,k)
          sy(jj,kk) = sy(j,k)
          dxi(jj,kk,1) =dxi(j,k,1)
          deta(jj,kk,1) = deta(j,k,1)
          dxi(jj,kk,2) = dxi(j,k,2)
          deta(jj,kk,2) = deta(j,k,2)
          xc(jj,kk) = 2.D0*xc(j,k) - xc(j2,k2)
          yc(jj,kk) = 2.D0*yc(j,k) - yc(j2,k2)
          zc(jj,kk) = 2.D0*zc(j,k) - zc(j2,k2)
          xnc(jj,kk) = 2.D0*xnc(j,k) - xnc(j2,k2)
          vegc(jj,kk) = vegc(j,k)
	      enddo
      enddo
      enddo


      return
      
 210  format(' ', A14, f12.2)
      end 
      
************************************************************************
      subroutine findbc(i,j,k,jj,kk,j2,k2)
      include 'dry.inc'
      if(ipos(j,k,i) .eq. 1) then  
         jj = j
         kk = k-1
	       j2 = j
         k2 = k+1
      elseif(ipos(j,k,i) .eq. 2) then
         jj = j+1
         kk = k 
	       j2 = j-1
         k2 = k
       elseif(ipos(j,k,i) .eq. 3) then 
         jj = j
         kk = k+1
	       j2 = j
         k2 = k-1
       elseif(ipos(j,k,i) .eq. 4) then 
         jj = j-1
         kk = k
	       j2 = j+1
	       k2 = k 
      endif
      return
      end   
************************************************************************
      subroutine limitr(i,beta,dq1,dq2,dq)
C     					
      implicit real*8(a-h,o-z)
C Lax-Wendroff.
      if(i .eq. 1)then
        dq = dq2 
C Beam-Warming.
      elseif(i .eq. 2)then
        dq = dq1 
C Fromm
      elseif(i .eq. 3)then
        dq = 0.5D0*(dq1 + dq2) 
C Double minmod.
      elseif(i .eq. 4)then
        a = 0.5D0*(dq1 + dq2)
        b = 2.D0*dq1
        c = 2.D0*dq2
        if(a*b .gt. 0.D0 .and. b*c .gt. 0.D0) then
           dq = fmin1(a, b, c)  
         else
           dq = 0.D0
        endif
C Beta Family.
      else
        if(dq1*dq2 .le. 0.D0) then
           dq = 0.D0
         else
           dq = fmin2(fmax2(dq1, dq2), beta*fmin2(dq1, dq2))
        endif
      endif

      return     
      end      
************************************************************************
      subroutine init
             
      include 'dry.inc' 
      ! Read  file 'input/params.dat'.
      read(4,'(a72)') dum
      read(4,*) grav, dt
      read(4,'(a72)') dum
      read(4,*) tmax, tr
      read(4,'(a72)') dum
      read(4,*) prate, nt
      read(4,'(a72)') dum      
      read(4,*) epsh, beta
      read(4,'(a72)') dum
      read(4,*) nprt 
      read(4,'(a72)') dum
      read(4,*) iscale
      read(4,'(a72)') dum
      read(4,*) stop_tol0   
      read(4,'(a72)') dum
      read(4,*) h0, u0, v0
      read(4,'(a72)') dum
      read(4,*) r8mV,r8etaV, r8alphaV  
      read(4,'(a72)') dum
      read(4,*) r8mB, r8etaB, r8alphaB
      read(4,'(a72)') dum
      read(4,*) jveg,  kveg               
      read(4,'(a72)') dum      
      read(4,*) jbare,  kbare
      read(4,'(a72)') dum 
      read(4,*) iBC, tsat_min


      ! Read  file 'input/boundary.dat'.
      read(3,'(a72)') dum
      read(3,*) nbcell
      read(3,'(a72)') dum
      do ii=1,nbcell
         read(3,*)j,k,inum(j,k),
     &   (itype(j,k,i),i=1,inum(j,k)),(ipos(j,k,i),i=1,inum(j,k))
      enddo
      read(3,'(a72)') dum
      read(3,*) ncol
      read(3,'(a72)') dum
      read(3,*) nrow
      read(3,'(a72)') dum
      do j=1,ncol       
         read(3,*) idum, kbeg(j), kend(j)
      enddo
      kbeg(ncol+1) = kbeg(ncol)
      kend(ncol+1) = kend(ncol)		
      read(3,'(a72)') dum
      read(3,*) ndir
      read(3,'(a72)') dum
      do i=1,ndir
       read(3,*) j,k,fix(j,k,1),fix(j,k,2),fix(j,k,3)
      enddo

      isurf = 2    ! 2 to set flow depth
      ilim =  5			
C Set up computational grid.
      call grid

C Set initial conditions.    
      t = 0.D0
      do j=1,ncol
        do k=kbeg(j),kend(j)
          
          n1 = nop(j,k,1)
          n2 = nop(j,k,2)
          n3 = nop(j,k,3)
          n4 = nop(j,k,4)
          
          if(isurf .eq. 1) then
            h(j,k) = h0 - zc(j,k)
           else
            h(j,k) = h0
          endif
          
          u(j,k) = u0
          v(j,k) = v0
        enddo
      enddo     

      do j=1,ncol
      do k=kbeg(j),kend(j)
C       if (h(j,k) .lt. 0.D0) h(j,k) = 0.D0
        q(j,k,1) = h(j,k)
        q(j,k,2) = h(j,k)*u(j,k)
        q(j,k,3) = h(j,k)*v(j,k)
        
        do i=1,inum(j,k)
C For fixed flux BCs, must initialize depth.
          if(itype(j,k,i) .eq. 4) then
            
            call findbc(i,j,k,jj,kk,j2,k2)
            if(ipos(j,k,i) .eq. 1 .or. ipos(j,k,i) .eq. 3) then
              qflux = fix(j,k,2)*cn(j,k,2) + fix(j,k,3)*sn(j,k,2)
              dx = -dxi(j,k,2)*area(j,k)
              dy = dxi(j,k,1)*area(j,k)  
              dss = dsqrt(dx*dx + dy*dy)
							
              if(h(j,k) .lt. epsh) then
                 write(*,*) ' *** bed adjacent to the boundary ',
     &                 'is dry ***'
                 if(qflux*dzz(j,k,2).lt.0.D0) then
                   qflux = dabs(qflux)
                   hnorm=(qflux*xnc(j,k)/dsqrt(dabs(dzz(j,k,2)/dss)))
     &           	          **(3.D0/5.D0)     
		 							
                   write(*,*) ' *** normal depth = ', hnorm*100.,
     &                              'cm is specified'
                   h(jj,kk) = hnorm 
                   hp(jj,kk) = hnorm     
	              else
                  write(*,*)  qflux, dzz(j,k,2)
                  write(*,*) 'adverse slope or zero Manning n in cell'
                  write(*,*) 'enter initial flow depth at specified'
                  write(*,*) 'flux boundary ', j,k,ipos(j,k,i)
!                 read(*,*) h(jj,kk)
                  hnorm  = 0.001
                  h(jj,kk) = hnorm                         
                  hp(jj,kk) = hnorm
	              endif
              endif
              qflux = dabs(qflux)
              
            endif
          endif
	      enddo
        
      enddo
      enddo
      
!     Initialize common time variables 
!     Make this every second instead
!       read(303,*) nt
!         do i=0,nt
!           read(303,*) tt(i), prate(i)
!         enddo
!
      
!     Define scaled time variables  
      dtinfl = dt*iscale  
      
      if (ntp .lt. nt) then 
        write(*,*) 'nt too small'
        write(*,*) 'ntp = ', ntp
        write(*,*) 'nt = ', nt                
        stop
      endif
                   
      
      return
 210  format(' ', A14, f12.3)
 211  format(' ', A14, f12.3, A4) 
 212  format(' ', A14, f12.3, A2)  
      end
      

      
************************************************************************
      subroutine initVG
    !   Specify van Genuchten soil parameters
    !   Creat time variables
    !   Input:      
    !   Output:           
    !   
    !   Comments:  
    !             
      include 'dry.inc' 

      real ( kind = 8 ) hinit(nz,2), thetainit(nz,2), 
     &                  r8Kinit(nz,2), gCinit(nz,2)

      inv = 0          ! inv = 0 for pseudo -inverse 
      tp = 1e10        !  time of ponding  (placeholder value) 
             
      htop =  -0.0d0

      !  option to read soil parameters from a file: 
      read(304,'(a72)') dum 
      read(304,*) dz, zmax
      read(304,'(a72)') dum 
      read(304,*) nzp
      read(304,'(a72)') dum 
      read(304,*) r8L	
					
      read(304,'(a72)') dum 			
      do i=1,nz
        read(304,*) alpha(i,1), theta_S(i,1), theta_R(i,1),
     &           r8lambda(i,1), r8Ksat(i,1), hinit(i,1)
      enddo

      read(304,'(a72)') dum 
 
      do i=1,nz
        read(304,*) alpha(i,2), theta_S(i,2), theta_R(i,2),
     &           r8lambda(i,2),r8Ksat(i,2), hinit(i,2)
      enddo

      r8n(1:nz,1) = r8lambda(1:nz,1) + 1.d0   !  n --> r8n
      r8m(1:nz,1) = r8lambda(1:nz,1)/r8n(1:nz,1)       ! m -->  r8m 
      r8n(1:nz,2) = r8lambda(1:nz,2) + 1.d0   !  n --> r8n
      r8m(1:nz,2) = r8lambda(1:nz,2)/r8n(1:nz,2)       ! m -->  r8m 

      !  set up grid
!       dz = 0.1d0
      z(1) = 0
      do j=1,nz-1
        z(j+1) = z(j) + dz
      enddo


C     Define matrices that we'll need in solution      
      DeltaPlus = 0.d0
      DeltaMinus = 0.d0
      do i = 2,nz-1
        DeltaPlus(i,i) = -1.d0
        DeltaPlus(i-1,i) = 1.d0  
        DeltaMinus(i,i) = 1.d0
        DeltaMinus(i,i-1) = -1.d0
      enddo
      DeltaPlus(1,2) = 0.d0
      DeltaPlus(nz-1, nz) = 1.d0
      
			APlus = 0.d0
			AMinus = 0.d0
      do i = 2,nz-1
        APlus(i,i) = 1.d0
        APlus(i-1,i) = 1.d0  
      enddo
      APlus(1,1) = 2.d0  ! MPlus --> APlus
      APlus(nz,nz) = 2.d0  
      APlus(nz-1,nz) = 1.d0            
      APlus(1,2) = 0.d0
            
      APlus(nz-1, nz-1) = 2.d0
      APlus(nz-1, nz) = 0.d0
      
      do i = 2,nz-1
        AMinus(i,i) = 1.d0
        AMinus(i,i-1) = 1.d0        
      enddo      
      AMinus(1,1) = 2.d0 ! MMinus --> AMinus
      AMinus(nz,nz) = 2.d0    
             
      do k=1,nz
        do i=1,2
          call vanGenuchten(k, hinit(k,i), 
     &         gCinit(k,i), r8Kinit(k,i), thetainit(k,i),i)
      enddo
      enddo

			r8H = 0.d0
			r8THETA = 0.d0
			r8K = 0.d0
			oldTHETA = 0.d0
      do j=1,ncol
        do k=kbeg(j),kend(j)
          if (vegc(j,k) .gt. 0) then
            r8H(j,k,1:nz) = hinit(1:nz,1)
            r8THETA(j,k,1:nz) = thetainit(1:nz,1)  
            r8K(j,k,1:nz) = r8Kinit(1:nz,1)
            oldTHETA(j,k,1:nz) = thetainit(1:nz,1)
          else
            r8H(j,k,1:nz) = hinit(1:nz,2)
            r8THETA(j,k,1:nz) = thetainit(1:nz,2)  
            r8K(j,k,1:nz) = r8Kinit(1:nz,2)
            oldTHETA(j,k,1:nz) = thetainit(1:nz,2)
          endif
        enddo
      enddo
      

      return
      end

************************************************************************
      subroutine vanGenuchten(k,gh,gC,gK,gtheta,iveg)
      include 'dry.inc'      

      !   vanGenuchten parameters
      !
      !   Input:
      !             k  (integer)  -  soil layer                          
      !             gh (real, kind = 8)    -  dimension
      !   Output:
      !             gtheta  (real, kind= 8) - volumetric moisture content
      !             gK  (real, kind= 8) - hydraulic conductivity
      !             gC  (real, kind= 8) - specific moisture storage
      !   
      !   Comments:   uses common variables theta_S, theta_R
      !             
			
      if (gh .lt. 0.d0) then
    
        gtheta = theta_R(k,iveg) +(theta_S(k,iveg)-theta_R(k,iveg))/
     &      (1.d0+(alpha(k,iveg)*abs(gh))**r8n(k,iveg))**r8m(k,iveg)
C     Compute the effective saturation      
        Se =  (gtheta - theta_R(k,iveg))/(theta_S(k,iveg) -
     &            theta_R(k,iveg))   
        
        gK =  r8Ksat(k,iveg)*Se**(r8L)*(1.d0-(1.d0-
     &            Se**(1.d0/r8m(k,iveg)))**r8m(k,iveg))**2.d0        
        
        gC =  alpha(k,iveg)*r8n(k,iveg)*(1.d0/r8n(k,iveg)-1.d0)*
     &        (alpha(k,iveg)*abs(gh))**(r8n(k,iveg)-1.d0)*
     &        (theta_R(k,iveg)-theta_S(k,iveg))*((alpha(k,iveg)
     &         *abs(gh))**r8n(k,iveg)+1)**(1/r8n(k,iveg)-2)
      
      elseif (gh .ge. 0.d0) then
        
        gtheta =  theta_S(k,iveg)
        gK = r8Ksat(k,iveg)
        gC = 0.d0
        
      endif


      return
      end   

************************************************************************
      subroutine potential( hnp1m_dum,thetan_dum, potI)
    !   Calculate potential infiltration with surface h = 0
    !   Input:
    !             hnp1m_dum,thetan_dum (real, kind = 8)    -  ICs
    !   Output:
    !             potential infiltration
            
      include 'dry.inc'
      integer ( kind = 4 )  info
      real ( kind = 8 )  r8gkt
      real ( kind = 8 )  hnp1m_dum(nz), thetan_dum(nz)      ! input   
      real ( kind = 8 )  hnp1m(nz), thetan(nz),     
     &          r8cnp1m(nz), r8knp1m(nz),  r8thetanp1m(nz),  
     &          hnp1mp1(nz), r8knp1mp1(nz),     ! time step m+1
     &          C(nz, nz), barKplus(nz), barKminus(nz),  !  matrices used by solver
     &          barplusK(nz, nz), barminusK(nz, nz), 
     &          dumPlus(nz, nz), dumMinus(nz, nz),     
     &          dumPlusH(nz), dumMinusH(nz), A(nz, nz),       
     &          R_MPFD(nz), Ainv(nz, nz), deltam(nz),
     &          uu(nz, nz), vv(nz, nz), s(nz, nz)
    
        istop_flag = 0
        niter = 0 
        stop_tol = stop_tol0
        hnp1m = hnp1m_dum 
        hnp1m(nz) = 0.0
        thetan = thetan_dum

        do while (istop_flag.eq.0)

          do k=1,nz             
            call vanGenuchten(k,hnp1m(k),
     &               r8cnp1m(k),r8knp1m(k),r8thetanp1m(k),isveg)
          enddo
          call fdiag(r8cnp1m, C)   
          !if (r8cnp1m(10) .ne. C(10, 10)) write (*,*) 'fdiag fail!'
          
          barKplus = 0.5d0*matmul(APlus, r8knp1m)  ! kbarplus -> barKplus
          call fdiag(barKplus, barplusK)  ! Kbarplus -> barplusK
          barKminus =  0.5d0*matmul(AMinus, r8knp1m)  ! kbarminus --> barKminus
          call fdiag(barKminus, barminusK)  ! Kbarminus --> barminusK
    
          dumPlus = matmul(barplusK, DeltaPlus)  ! Kbarplus.dot(DeltaPlus)
          dumPlusH = matmul(dumPlus, hnp1m)  ! (Kbarplus.dot(DeltaPlus).dot(hnp1m)
          dumMinus = matmul(barminusK, DeltaMinus) 
          dumMinusH = matmul(dumMinus, hnp1m)
    
          A = (1.d0/dtinfl)*C-1.d0/(dz**2.d0)*(dumPlus - dumMinus) 
          
C         Compute the residual of MPFD (RHS)
          R_MPFD = (1.d0/(dz**2.d0))*(dumPlusH - dumMinusH) +
     &         (1.d0/dz)*(barKplus - barKminus) - 
     &          1.d0/dtinfl*(r8thetanp1m - thetan)
C         Compute deltam for iteration level m+1

          call r8mat_svd_lapack ( nz, nz, A, uu, s, vv, info)

          call pseudo_inverse ( nz, nz, uu, s, vv, Ainv )

          deltam = matmul(Ainv, R_MPFD)
          niter = niter + 1           
          
          t2b_theta = (thetan(nz) - thetan(1))
          r8mean_theta = sum(thetan)/nz
          r8var_theta = sum((thetan - r8mean_theta)**2)/nz
          r8std_theta = sqrt(r8var_theta)


          r8mean_H = sum(hnp1m)/nz
          r8var_H = sum((hnp1m - r8mean_H)**2)/nz
          r8std_H  = sqrt(r8var_H)
          r8std_tol = 1.0

C         Saturation catch!  before allowing the soil to saturate: 
C             wait for at least tsat_min minutes, and check that depth > 0
          if ((t .gt. tsat_min) .and. (r8std_H  .lt. r8std_tol))  then  
C          if soil moisture difference between top and bottom is < 0.01
              potI =  - r8knp1m(nz) 
              istop_flag = 1  ! exit Richards solver - don't wait for convergence
           	  return
          endif

          if (niter .gt. 100) then
            write(*,*) 'PI:niter > 100,stop_tol = ',stop_tol						 
            if (stop_tol .gt. 100.0) then
              write(*,*) 'stop tol too big'
              return
            endif
            stop_tol = stop_tol*10
            niter = 0
          endif

          if (maxval(abs(deltam(2:(nz-1)))).lt. stop_tol) then
            istop_flag = 1 
            hnp1mp1 = hnp1m + deltam              
            hnp1mp1(nz) = 0.d0 
            hnp1mp1(1) = hnp1mp1(2)
            do k=1,nz         
              call vanGenuchten(k, hnp1mp1(k),  
     &              r8cnp1m(k), r8knp1m(k), r8thetanp1m(k), isveg)    !   knp1m --> r8knp1m
            enddo
            r8knp1mp1 = r8knp1m
            hnp1m = hnp1mp1
            r8gkt =  r8knp1mp1(nz-1)
            potI = r8gkt*((hnp1m(nz) - hnp1m(nz-1))/dz + 1.d0)

          else
              hnp1mp1 = hnp1m + deltam  
              hnp1m = hnp1mp1  
              hnp1m(nz) = 0.0
              hnp1m(1) = hnp1m(2)
          endif
          
        enddo
      return
      end   

************************************************************************
      subroutine timestep(hnp1m,thetan,hnp1mp1,r8thetanp1m,r8knp1mp1)
    !
    !   Input:
    !             hnp1m, thetan  (real, kind = 8)   - initial
    !   Output:
    !             hnp1mp1,r8thetanp1m,r8knp1mp1   -  h, theta and k at time m+1
    !             
    !  Comments:  uses common variables  nz, stop_tol, htop as surface h, 
    !             modifies the common variable ipass
      include 'dry.inc'
      
      real ( kind = 8 )   r8gkt, htop                   !  surface k
      real (kind = 8)   r8mean_H, r8var_H, r8std_H,r8std_tol
      real ( kind = 8 )   hnp1m(nz), thetan(nz),   ! input arrays
     &          r8cnp1m(nz), r8knp1m(nz),         !  output arrays
     &          hnp1mp1(nz), r8knp1mp1(nz), r8thetanp1m(nz),     
     &          C(nz, nz), barKplus(nz), barKminus(nz),
     &          barplusK(nz, nz), barminusK(nz, nz), 
     &          dumPlus(nz, nz), dumMinus(nz, nz),     
     &          dumPlusH(nz), dumMinusH(nz), A(nz, nz),       
     &          R_MPFD(nz), Ainv(nz, nz), deltam(nz),
     &          uu(nz, nz), vv(nz, nz), s(nz, nz)

      istop_flag = 0   ! indicator variable switches to 1 when convergence criteria is met
      niter = 0  			 ! number of iterations
      stop_tol = stop_tol0    ! 
      ipass = 0                   ! skip Richards solver for saturated flow

      do while (istop_flag .eq. 0)
        do k=1,nz    
            call vanGenuchten(k,hnp1m(k),
     &      			r8cnp1m(k),r8knp1m(k),r8thetanp1m(k), isveg)
          enddo
          
          call fdiag(r8cnp1m, C)   
          barKplus = 0.5d0*matmul(APlus, r8knp1m)  ! kbarplus -> barKplus
          call fdiag(barKplus, barplusK)  ! Kbarplus -> barplusK
          barKminus =  0.5d0*matmul(AMinus, r8knp1m)  ! kbarminus --> barKminus
          call fdiag(barKminus, barminusK)  ! Kbarminus --> barminusK
            
          dumPlus = matmul(barplusK, DeltaPlus)  ! Kbarplus.dot(DeltaPlus)
          dumPlusH = matmul(dumPlus, hnp1m)  ! (Kbarplus.dot(DeltaPlus).dot(hnp1m)
          dumMinus = matmul(barminusK, DeltaMinus) 
          dumMinusH = matmul(dumMinus, hnp1m)

          A = (1.d0/dtinfl)*C-1.d0/(dz**2.d0)*(dumPlus - dumMinus)
    
C         Compute the residual of MPFD (RHS)
          R_MPFD = (1.d0/(dz**2.d0))*(dumPlusH - dumMinusH) +
     &         (1.d0/dz)*(barKplus - barKminus) - 
     &          1.d0/dtinfl*(r8thetanp1m - thetan)

C         Compute deltam for iteration level m+1
          call r8mat_svd_lapack ( nz, nz, A, uu, s, vv, info )
          
          call pseudo_inverse ( nz, nz, uu, s, vv, Ainv )

          deltam = matmul(Ainv, R_MPFD)
          niter = niter + 1

C          t2b_theta = (thetan(nz) - thetan(1))
C           r8mean_theta = sum(thetan)/nz
C           r8var_theta = sum((thetan - r8mean_theta)**2)/nz
C           r8std_theta = sqrt(r8var_theta)

          r8mean_H = sum(hnp1m)/nz
          r8var_H = sum((hnp1m - r8mean_H)**2)/(nz-1)
          r8std_H  = sqrt(r8var_H)
          r8std_tol = 1.0

C         Saturation catch!  before allowing the soil to saturate: 
C             wait for at least tsat_min minutes, and check that depth > 0
          if ((t .gt. tsat_min) .and. (r8std_H  .lt. r8std_tol))  then  
C          if soil moisture difference between top and bottom is < 0.01
            if (htop .gt.0) then  

              hnp1m(nz) = 0
              hnp1mp1 = hnp1m     ! update h to n+1,m+1
              r8thetanp1m = thetan ! update theta to n,m+1

              do k = 1,nz
                call vanGenuchten(k, hnp1mp1(k),  r8cnp1m(k),
     &              r8knp1mp1(k), r8thetanp1m(k), isveg)
              enddo

C             set surface flux to K(surface) and adjust H to match 
              flux =   - r8knp1mp1(nz)
C               hnp1mp1(nz) = hnp1m(nz) !  - flux*dt ! ?
              hnp1mp1(nz-1) = hnp1mp1(nz) + dz + flux*dz/r8knp1m(nz-1)

              ipass = 1    
              istop_flag = 1  ! exit Richards solver - don't wait for convergence

            elseif ((r8knp1m(nz).ge.prate/100) .and.(prate.gt.0)) then

              hnp1m(nz) = 0
              hnp1mp1 = hnp1m     ! update h to n+1,m+1
              r8thetanp1m = thetan ! update theta to n,m+1
              do k = 1,nz
                call vanGenuchten(k, hnp1mp1(k),  r8cnp1m(k),
     &              r8knp1mp1(k), r8thetanp1m(k), isveg)
              enddo
C               write(*, *)  'H',hnp1m
C             set surface flux to K(surface) and adjust H to match 
              flux =   - r8knp1mp1(nz)
C               hnp1mp1(nz) = hnp1m(nz) !  - flux*dt ! ?
              hnp1mp1(nz-1) = hnp1mp1(nz)+ dz + flux*dz/r8knp1m(nz-1)


              ipass = 1    
              istop_flag = 1  ! exit Richards solver - don't wait for convergence

	           elseif (prate.le.1e-10) then

	              hnp1mp1 = hnp1m     ! update h to n+1,m+1
	              r8thetanp1m = thetan ! update theta to n,m+1

C	            write(*, *)  'H',hnp1m
C             set surface flux to K(surface) and adjust H to match 
              flux =   - r8knp1m(nz)
              hnp1mp1(nz) = hnp1m(nz) !  - flux*dt ! ?
              hnp1mp1(nz-1) = hnp1mp1(nz)+ dz + flux*dz/r8knp1m(nz-1)

              do k = 1,nz
                call vanGenuchten(k, hnp1mp1(k),  r8cnp1m(k),
     &              r8knp1mp1(k), r8thetanp1m(k), isveg)
              enddo
              ipass = 1    
              istop_flag = 1  ! exit Richards solver - don't wait for convergence
            endif

          endif

          if (niter .gt. 100) then 
            stop_tol = stop_tol*10
            niter = 0
            write(*,*) 'timestep: stop_tol=',stop_tol,  r8std_H,
     &               't=', t/60., hnp1m(1), hnp1m(nz), htop, isetflux	
            if (stop_tol .gt. 100.0) then
              write(*,*) 'stop tol too big'
C               write(102,*) 'stop_tol too big!'
C              call gracefulExit
              return
            endif
          endif
              
          if (ipass .eq. 0) then  ! if not saturated... execute richards solver

          if (maxval(abs(deltam(2:(nz-1)))).lt. stop_tol) then
            istop_flag = 1  
            hnp1mp1 = hnp1m + deltam              
            if (isetflux .eq. 0) then
              hnp1mp1(nz) = 0.0  !htop
            elseif (isetflux .eq. 1) then
              r8gkt = r8knp1m(nz-1)
              hnp1mp1(nz)  =  hnp1mp1(nz-1) - dz - flux*dz/r8gkt
            endif
            hnp1mp1(1) = hnp1mp1(2)
            
            do k=1,nz         
              call vanGenuchten(k, hnp1mp1(k),  
     &         r8cnp1m(k),r8knp1m(k),r8thetanp1m(k), isveg)    !   knp1m --> r8knp1m
            enddo
            r8knp1mp1 = r8knp1m
            ! hnp1m = hnp1mp1
            ! r8gkt =  r8knp1mp1(nz-1)                
            !flux = r8gkt*((hnp1m(nz) - hnp1m(nz-1))/dz + 1.d0)    
          
          else

            hnp1mp1 = hnp1m + deltam  
            hnp1m = hnp1mp1  
            if (isetflux .eq. 0) then
              
              hnp1m(nz) = 0.0  !htop
              
            elseif (isetflux .eq. 1) then
              
              r8gkt = r8knp1m(nz-1)
              hnp1m(nz) = hnp1m(nz-1) - dz-flux*dz/r8gkt
                
            endif
            hnp1m(1) = hnp1m(2)
          endif
          endif
C           if  (hnp1m(nz) .gt. 0.00)  then
C              write(*,*) t, hnp1m(nz),hnp1m(nz-1),htop,ipass,isetflux
C           endif
        enddo
        ! error estimate - not used for anythin
        r8error =  sum(matmul(A, deltam))*dt*dz*iscale

      return
 
      end 

************************************************************************
      real*8 function f1(h,up)
      implicit real*8(a-h,o-z)
      f1 = h*up
      return
      end

************************************************************************
      real*8 function f2(h,u,up,cn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh
      f2 = h*u*up + 0.5D0*grav*h*h*cn
      return
      end

************************************************************************
      real*8 function f3(h,v,up,sn)
      implicit real*8(a-h,o-z)
      common/m/ grav, amax, epsh
      f3 = h*v*up + 0.5D0*grav*h*h*sn
      return
      end
************************************************************************
      real*8 function fmin1(a,b,c)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
         fmin1 = -1.D0*dmin1(dabs(a),dabs(b),dabs(c)) 
        else
         fmin1 = dmin1(a,b,c) 
      endif 
      return
      end
************************************************************************
      real*8 function fmin2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmin2 = -dmin1(dabs(a),dabs(b))
         else
          fmin2 = dmin1(a,b)
      endif 
      return
      end
************************************************************************
      real*8 function fmax2(a,b)
      implicit real*8(a-h,o-z)
      if (a .lt. 0.D0) then
          fmax2 = -dmax1(dabs(a),dabs(b))
         else
          fmax2 = dmax1(a,b) 
      endif
      return
      end

************************************************************************
      subroutine fdiag(aa, bb)
    ! Make a diagonal 2D matrix (hh) from 1D array (dh)
    ! 
    ! Input:
    !             dh  (real,)
    !             gh (real, kind = 8)    -  dimension
    ! Output:
    !             gtheta  (real, kind= 8) - volumetric moisture content
    !             gK  (real, kind= 8) - hydraulic conductivity
      include 'dry.inc'
      
      real(kind=8), dimension(nz) :: aa
      real(kind=8), dimension(nz, nz) :: bb
      
      bb = 0.0
      do k=1,nz
        bb(k,k) = aa(k)
      enddo

      return
      end
      
************************************************************************
      subroutine r8mat_svd_lapack ( m, n, a, u, s, v, info )
    !
    !! R8MAT_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
    !
    !
    !    The singular value decomposition of a real MxN matrix A has the form:
    !
    !      A = U * S * V'
    !
    !    where
    !
    !      U is MxM orthogonal,
    !      S is MxN, and entirely zero except for the diagonal;
    !      V is NxN orthogonal.
    !
    !    Moreover, the nonzero entries of S are positive, and appear
    !    in order, from largest magnitude to smallest.
    !
    !    This routine calls the LAPACK routine DGESVD to compute the
    !    factorization.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Output, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
    !    that form the singular value decomposition of A.
    !
      implicit none

      integer ( kind = 4 ) m, n
      real ( kind = 8 ) a(m,n), a_copy(m,n)
      integer ( kind = 4 ) i, info, lda, ldu, ldv
      character jobu, jobv
      integer ( kind = 4 ) lwork
      real ( kind = 8 ) sdiag(min(m,n))
      real ( kind = 8 ) s(m,n)
      real ( kind = 8 ) u(m,m)
      real ( kind = 8 ) v(n,n)
      real ( kind = 8 ), allocatable, dimension ( : ) :: work

      lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

      allocate ( work(1:lwork) )
    !
    !  Compute the eigenvalues and eigenvectors.
    !
      jobu = 'A'
      jobv = 'A'
      lda = m
      ldu = m
      ldv = n
    !
    !  The input matrix is destroyed by the routine.  Since we need to keep
    !  it around, we only pass a copy to the routine.
    !
      a_copy(1:m,1:n) = a(1:m,1:n)

      call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v,ldv,
     &    work, lwork, info )
     
      if ( info /= 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info

        return
      end if
    !   
    !  Make the MxN matrix S from the diagonal values in SDIAG.
    !
      s(1:m,1:n) = 0.0D+00
      do i = 1, min ( m, n )
        s(i,i) = sdiag(i)
      end do
    !
    !  Transpose V.
    !
      v = transpose ( v )

      deallocate ( work )

      return
      end

************************************************************************
      subroutine pseudo_inverse ( m, n, u, s, v, a_pseudo )
    !
    !  PSEUDO_INVERSE computes the pseudoinverse.
    !
    !    Given the singular value decomposition of a real MxN matrix A:
    !
    !      A = U * S * V'
    !
    !    where
    !
    !      U is MxM orthogonal,
    !      S is MxN, and entirely zero except for the diagonal;
    !      V is NxN orthogonal.
    !
    !    the pseudo inverse is the NxM matrix A+ with the form
    !
    !      A+ = V * S+ * U'
    !
    !    where
    !
    !      S+ is the NxM matrix whose nonzero diagonal elements are
    !      the inverses of the corresponding diagonal elements of S.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
    !    that form the singular value decomposition of A.
    !
    !    Output, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
    !
      implicit none

      integer ( kind = 4 ) m
      integer ( kind = 4 ) n

      real ( kind = 8 )    a_pseudo(n,m)
      integer ( kind = 4 ) i
      real ( kind = 8 ) s(m,n)
      real ( kind = 8 ) sp(n,m)
      real ( kind = 8 ) u(m,m)
      real ( kind = 8 ) v(n,n)

      sp(1:n,1:m) = 0.0D+00
      do i = 1, min ( m, n )
        if ( s(i,i) /= 0.0D+00 ) then
          sp(i,i) = 1.0D+00 / s(i,i)
        end if
      end do

      a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n),
     &    matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )

      return
      end
  
