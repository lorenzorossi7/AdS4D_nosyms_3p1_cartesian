c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different from ex at grid points that are at large rho=sqrt(x**2+y**2+z**2) and next to excised points, so we're identifying the last not excised grid points near the boundary. We will use them to extrapolate the value of the quasi-local boundary stress-energy tensor at the boundary
c----------------------------------------------------------------------

        subroutine nexttobdypoints(
     &                  chrbdy,
     &                  numbdypoints,
     &                  x,y,z,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1
        real*8 maxxyzp1
        integer numbdypoints

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        numbdypoints=0

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if

            if (rhop1.lt.(1-5*dx/2)) then !points between rhobdy=1 and 1-dx/2 are excised, points between rhobdy=1 and 1-3*dx/2 use forward/backward stencils (so when we cannot expect convergence at these points because the stencils used are different for different resolutions), points between rhobdy=1 and 1-5*dx/2 use points that are set by forward/backward stencils, so we cannot expect convergence. If we want to check convergence at the grid points used for extrapolation at the boundary, we need to pick points that have rho<1-5*dx/2 
              chrbdy(i,j,k) =ex-1.0d0
              chrbdy2(i,j,k)=ex-1.0d0
            else
              chrbdy(i,j,k) =ex
              chrbdy2(i,j,k)=ex
            end if

! eliminate troublesome points
           if (chrbdy(i,j,k).ne.ex) then
            if (
     &          (abs(zp1).lt.10.0d0**(-10))  !removing z=0 implies that we remove in particular the troublesome points with chi=0,1 (which have y=z=0,x=1,-1) and points with xi=0,1 (which have z=0,y>0,any x). We will fill and impose regularity at these points in Mathematica.
     &         ) then
             chrbdy(i,j,k)=ex
            end if
           end if

           end do
          end do
         end do


         do i=is,ie
          do j=js,je
           do k=ks,ke

           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)

           !chrbdy(i,j,k) is not ex only for points near the boundary AND next to excised points   

           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

           if (chrbdy(i,j,k).ne.ex) then
            if (maxxyzp1.eq.abs(xp1)) then
!if we use derivatives to define near boundary quantities, we will only define them at points between is and ie (js and je, ks and ke).
!Therefore, for extrapolation, we can only select near boundary points whose neighbors in the direction of the bulk along the axes (i.e. the direction of extrapolation) are within that range
             if (xp1.gt.0) then 
              if (((i-1).lt.is).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &           ((chrbdy2(i+1,j,k).ne.ex)
     &           .and.(chrbdy2(i-1,j,k).ne.ex))) then
                   chrbdy(i,j,k)=ex
              end if
             else
              if (((i+1).gt.ie).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &           ((chrbdy2(i+1,j,k).ne.ex)
     &           .and.(chrbdy2(i-1,j,k).ne.ex))) then
                   chrbdy(i,j,k)=ex
              end if
             end if


            else if (maxxyzp1.eq.abs(yp1)) then
             if (yp1.gt.0) then
              if (((j-1).lt.js).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &            ((chrbdy2(i,j+1,k).ne.ex)
     &           .and.(chrbdy2(i,j-1,k).ne.ex))) then
                   chrbdy(i,j,k)=ex
              end if
             else
              if (((j+1).gt.je).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &           ((chrbdy2(i,j+1,k).ne.ex)
     &           .and.(chrbdy2(i,j-1,k).ne.ex))) then
                   chrbdy(i,j,k)=ex
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if (((k-1).lt.ks).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &           ((chrbdy2(i,j,k+1).ne.ex)
     &           .and.(chrbdy2(i,j,k-1).ne.ex))) then
               chrbdy(i,j,k)=ex
              end if
             else
              if (((k+1).gt.ke).or. !this condition is actually only necessary when using derivatives to define near bdy quantities, i.e. when no_derivatives is set to .false. below
     &           ((chrbdy2(i,j,k+1).ne.ex)
     &           .and.(chrbdy2(i,j,k-1).ne.ex))) then
                   chrbdy(i,j,k)=ex
              end if
             end if
            
            end if
           

          end if

           if (chrbdy(i,j,k).ne.ex) then
             numbdypoints=numbdypoints+1
           end if


           end do
         end do
        end do

        return
        end


c------------------------------------------------------------------------------------------------------


c------------------------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different from ex at grid points that are at large rho=sqrt(x**2+y**2+z**2) and next to excised points, so we're identifying the last not excised grid points near the boundary. We will use them to extrapolate the value of the quasi-local boundary stress-energy tensor at the boundary
c-------------------------------------------------------------------------------------------------------------------------

        subroutine xyzextrap(
     &                  xextrap,yextrap,zextrap,
     &                  chrbdy,
     &                  numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,ghost_width)

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz
        integer lind

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1
        real*8 maxxyzp1
        integer numbdypoints
        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1
              if (maxxyzp1.eq.abs(xp1)) then
               if (xp1.gt.0) then
                  xextrap(lind)=sqrt(1-yp1**2-zp1**2)
                  yextrap(lind)=yp1
                  zextrap(lind)=zp1
              else
                  xextrap(lind)=-sqrt(1-yp1**2-zp1**2)
                  yextrap(lind)=yp1
                  zextrap(lind)=zp1
              end if
             else if (maxxyzp1.eq.abs(yp1)) then
              if (yp1.gt.0) then
                  yextrap(lind)=sqrt(1-xp1**2-zp1**2)
                  xextrap(lind)=xp1
                  zextrap(lind)=zp1
              else
                  yextrap(lind)=-sqrt(1-xp1**2-zp1**2)
                  xextrap(lind)=xp1
                  zextrap(lind)=zp1
              end if
             else
                 if (zp1.gt.0) then
                  zextrap(lind)=sqrt(1-yp1**2-xp1**2)
                  yextrap(lind)=yp1
                  xextrap(lind)=xp1
              else
                  zextrap(lind)=-sqrt(1-yp1**2-xp1**2)
                  yextrap(lind)=yp1
                  xextrap(lind)=xp1
               end if
              end if

!               xextrap(lind)=xp1   !TEST
!               yextrap(lind)=yp1   !TEST
!               zextrap(lind)=zp1   !TEST

            end if
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the leading order coefficient of the scalar field phi1 in the near bdy expansion
c in powers of q=1-rho about  q=0.
c----------------------------------------------------------------------

        subroutine calc_leadordcoeff_phi1(leadordcoeff_phi1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical no_derivatives
        data no_derivatives/.false./

        real*8 df_drho

        real*8 test1(Nx,Ny,Nz)
        real*8 dtest1_drho
!----------------------------------------------------------------------

        if (no_derivatives) then
         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0

            if ((chr(i,j,k).ne.ex)
     &          .and.(
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then
               leadordcoeff_phi1(i,j,k)=phi1_n(i,j,k)/q
            else
               leadordcoeff_phi1(i,j,k)=0
            end if
           end do
          end do
         end do

        else !i.e. if .not.no_derivatives
         ! set index bounds for main loop
         is=2
         ie=Nx-1
         js=2
         je=Ny-1
         ks=2
         ke=Nz-1
 
         ! adjust index bounds to compensate for ghost_width
         if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
         if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
         if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
         if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
         if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
         if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1) 

         do i=is,ie
          do j=js,je
           do k=ks,ke
             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0

            if ((chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then
               leadordcoeff_phi1(i,j,k)=
     &                  -df_drho(phi1_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
            else
               leadordcoeff_phi1(i,j,k)=0
            end if

           end do
          end do
         end do

        end if

        return
        end
c-----------------------------------------------------------------------------------


c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the leading order coefficient of the scalar field phi1 at the boundary
c----------------------------------------------------------------------

        subroutine calc_bdyphi(bdyphi,
     &                  leadordcoeff_phi1,
     &                  xextrap,yextrap,zextrap,
     &                  chrbdy,numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        real*8 chrbdy(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        integer i,j,k,is,ie,js,je,ks,ke,lind
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1
        real*8 leadordcoeff_phi1_p2
        real*8 bdyphi(numbdypoints)

        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 xp1,yp1,zp1,xp2,yp2,zp2
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 maxxyzp1

        real*8 extrapalongx,extrapalongy,extrapalongz
!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)


        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)

           leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

                   xex=xextrap(lind)
                   yex=yextrap(lind)
                   zex=zextrap(lind)
                   rhoex=sqrt(xex**2+yex**2+zex**2)
                   chiex=(1/PI)*acos(xex/rhoex)
                if (zex.lt.0) then
                   xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                else
                   xiex=(1/(2*PI))*atan2(zex,yex)
                end if


              if (maxxyzp1.eq.abs(xp1)) then
                if (xp1.gt.0) then
                  xp2=x(i-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i-1,j,k)
                  bdyphi(lind)=extrapalongx(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,xp1,xp2,xex)
              else
                  xp2=x(i+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i+1,j,k)
                  bdyphi(lind)=extrapalongx(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,xp1,xp2,xex)
               end if
             else if (maxxyzp1.eq.abs(yp1)) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j-1,k)
                  bdyphi(lind)=extrapalongy(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,yp1,yp2,yex)
              else
                  yp2=y(j+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j+1,k)
                  bdyphi(lind)=extrapalongy(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,yp1,yp2,yex)
              end if
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k-1)
                  bdyphi(lind)=extrapalongz(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,zp1,zp2,zex)
                 else
                  zp2=z(k+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k+1)
                  bdyphi(lind)=extrapalongz(leadordcoeff_phi1_p1
     &                         ,leadordcoeff_phi1_p2,zp1,zp2,zex)
               end if
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if
          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------

        


c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D
c using a 1-rho expansion about rho=1 at points near the boundary. 
c The tensor components are given in spherical polar coordinates.
c----------------------------------------------------------------------

        subroutine quasiset_ll(
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_tracell,
     &                  quasiset_massdensityll,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)


        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 psi_np1(Nx,Ny,Nz),psi_n(Nx,Ny,Nz),psi_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)
        real*8 quasiset_tracell(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)



        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t, gb_tt_x, gb_tt_y
        real*8 gb_tt_z
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty
        real*8 gb_tt_tz
        real*8 gb_tt_xx,gb_tt_xy
        real*8 gb_tt_xz
        real*8 gb_tt_yy
        real*8 gb_tt_yz
        real*8 gb_tt_zz


        real*8 gb_tx_t, gb_tx_x, gb_tx_y
        real*8 gb_tx_z
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty
        real*8 gb_tx_tz
        real*8 gb_tx_xx,gb_tx_xy
        real*8 gb_tx_xz
        real*8 gb_tx_yy
        real*8 gb_tx_yz
        real*8 gb_tx_zz

        real*8 gb_ty_t, gb_ty_x, gb_ty_y
        real*8 gb_ty_z
        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty
        real*8 gb_ty_tz
        real*8 gb_ty_xx,gb_ty_xy
        real*8 gb_ty_xz
        real*8 gb_ty_yy
        real*8 gb_ty_yz
        real*8 gb_ty_zz

        real*8 gb_tz_t, gb_tz_x, gb_tz_y
        real*8 gb_tz_z
        real*8 gb_tz_tt,gb_tz_tx,gb_tz_ty
        real*8 gb_tz_tz
        real*8 gb_tz_xx,gb_tz_xy
        real*8 gb_tz_xz
        real*8 gb_tz_yy
        real*8 gb_tz_yz
        real*8 gb_tz_zz


        real*8 gb_xx_t, gb_xx_x, gb_xx_y
        real*8 gb_xx_z
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty
        real*8 gb_xx_tz
        real*8 gb_xx_xx,gb_xx_xy
        real*8 gb_xx_xz
        real*8 gb_xx_yy
        real*8 gb_xx_yz
        real*8 gb_xx_zz

        real*8 gb_xy_t, gb_xy_x, gb_xy_y
        real*8 gb_xy_z
        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty
        real*8 gb_xy_tz
        real*8 gb_xy_xx,gb_xy_xy
        real*8 gb_xy_xz
        real*8 gb_xy_yy
        real*8 gb_xy_yz
        real*8 gb_xy_zz

        real*8 gb_xz_t, gb_xz_x, gb_xz_y
        real*8 gb_xz_z
        real*8 gb_xz_tt,gb_xz_tx,gb_xz_ty
        real*8 gb_xz_tz
        real*8 gb_xz_xx,gb_xz_xy
        real*8 gb_xz_xz
        real*8 gb_xz_yy
        real*8 gb_xz_yz
        real*8 gb_xz_zz

        real*8 gb_yy_t, gb_yy_x, gb_yy_y
        real*8 gb_yy_z
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty
        real*8 gb_yy_tz
        real*8 gb_yy_xx,gb_yy_xy
        real*8 gb_yy_xz
        real*8 gb_yy_yy
        real*8 gb_yy_yz
        real*8 gb_yy_zz

        real*8 gb_yz_t, gb_yz_x, gb_yz_y
        real*8 gb_yz_z
        real*8 gb_yz_tt,gb_yz_tx,gb_yz_ty
        real*8 gb_yz_tz
        real*8 gb_yz_xx,gb_yz_xy
        real*8 gb_yz_xz
        real*8 gb_yz_yy
        real*8 gb_yz_yz
        real*8 gb_yz_zz

        real*8 psi_t, psi_x, psi_y
        real*8 psi_z
        real*8 psi_tt,psi_tx,psi_ty
        real*8 psi_tz
        real*8 psi_xx,psi_xy
        real*8 psi_xz
        real*8 psi_yy
        real*8 psi_yz
        real*8 psi_zz


        real*8 dtdt,dtdrho,dtdchi,dtdxi
        real*8 dxdt,dxdrho,dxdchi,dxdxi
        real*8 dydt,dydrho,dydchi,dydxi
        real*8 dzdt,dzdrho,dzdchi,dzdxi

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_tx_ads0,g0_ty_ads0,g0_tz_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_psi_ads0
        real*8 g0_xz_ads0,g0_yz_ads0

        real*8 dxcar_dxsph(4,4)
        real*8 hcar_n(4,4),g0car_n(4,4)
        real*8 hsph_n(4,4),gbsph_n(4,4),g0sph_n(4,4)
        real*8 gamma0sph_ll(4,4),gamma0sph_uu(4,4)
        real*8 gamma0sphbdy_ll(3,3),gamma0sphbdy_uu(3,3)
        real*8 detgamma3
        real*8 detgamma0sphbdy

        logical no_derivatives
        data no_derivatives/.false./
 
        real*8 df_drho
        real*8 gbsph_tt_n(Nx,Ny,Nz),gbsph_trho_n(Nx,Ny,Nz)
        real*8 gbsph_tchi_n(Nx,Ny,Nz),gbsph_txi_n(Nx,Ny,Nz)
        real*8 gbsph_rhorho_n(Nx,Ny,Nz),gbsph_rhochi_n(Nx,Ny,Nz)
        real*8 gbsph_rhoxi_n(Nx,Ny,Nz)
        real*8 gbsph_chichi_n(Nx,Ny,Nz),gbsph_chixi_n(Nx,Ny,Nz)
        real*8 gbsph_xixi_n(Nx,Ny,Nz)
        real*8 dgbsph_tt_drho_n,dgbsph_trho_drho_n
        real*8 dgbsph_tchi_drho_n,dgbsph_txi_drho_n
        real*8 dgbsph_rhorho_drho_n
        real*8 dgbsph_rhochi_drho_n
        real*8 dgbsph_rhoxi_drho_n
        real*8 dgbsph_chichi_drho_n
        real*8 dgbsph_chixi_drho_n
        real*8 dgbsph_xixi_drho_n
        real*8 gamma0sphbdy_uu_tt(Nx,Ny,Nz)
        real*8 gamma0sphbdy_uu_tchi(Nx,Ny,Nz)
        real*8 gamma0sphbdy_uu_txi(Nx,Ny,Nz)
        real*8 gamma0sphbdy_uu_chichi(Nx,Ny,Nz)
        real*8 gamma0sphbdy_uu_chixi(Nx,Ny,Nz)
        real*8 gamma0sphbdy_uu_xixi(Nx,Ny,Nz)

        real*8 gbsph_tt_n_x,gbsph_tt_n_y
        real*8 gb_tt_n_x,gb_tt_n_y
        real*8 gb_xx_n_x,gb_xx_n_y

        real*8 dergb_tt_x_n,dergbsph_tt_x_n
        real*8 dergb_tt_y_n,dergbsph_tt_y_n
        real*8 dergb_tt_z_n,dergbsph_tt_z_n
        real*8 dergb_xx_x_n
        real*8 dergb_xx_y_n
        real*8 dergb_xx_z_n
        real*8 dergb_tt_x_np1,dergb_tt_y_np1,dergb_tt_z_np1

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

! define quantities necessary to define the quasi local stress energy tensor over the whole grid
       do i=1,Nx
        do j=1,Ny
         do k=1,Nz

             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if
!calculate regularized metric components in spherical coordinates in terms of regularized metric components in Cartesian coordinates
! we use the following coordinate transformation (notice that the angles are rescaled w.r.t. the usual spherical coordinates): x=rho*cos(PI*chi),y=rho*sin(PI*chi)*cos(2*PI*xi),z=rho*sin(PI*chi)*sin(2*PI*xi)

          if ( (chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then

!transformation matrix

             dtdt=1
             dtdrho=0
             dtdchi=0
             dtdxi=0
             dxdt=0
             dxdrho=cos(PI*chi0)
             dxdchi=-PI*rho0*sin(PI*chi0)
             dxdxi=0
             dydt=0
             dydrho=sin(PI*chi0)*cos(2*PI*xi0)
             dydchi=PI*rho0*cos(PI*chi0)*cos(2*PI*xi0)
             dydxi=-2*PI*rho0*sin(PI*chi0)*sin(2*PI*xi0)
             dzdt=0
             dzdrho=sin(PI*chi0)*sin(2*PI*xi0)
             dzdchi=PI*rho0*cos(PI*chi0)*sin(2*PI*xi0)
             dzdxi=2*PI*rho0*sin(PI*chi0)*cos(2*PI*xi0)

             dxcar_dxsph(1,1)=dtdt
             dxcar_dxsph(1,2)=dtdrho
             dxcar_dxsph(1,3)=dtdchi
             dxcar_dxsph(1,4)=dtdxi
             dxcar_dxsph(2,1)=dxdt
             dxcar_dxsph(2,2)=dxdrho
             dxcar_dxsph(2,3)=dxdchi
             dxcar_dxsph(2,4)=dxdxi
             dxcar_dxsph(3,1)=dydt
             dxcar_dxsph(3,2)=dydrho
             dxcar_dxsph(3,3)=dydchi
             dxcar_dxsph(3,4)=dydxi
             dxcar_dxsph(4,1)=dzdt
             dxcar_dxsph(4,2)=dzdrho
             dxcar_dxsph(4,3)=dzdchi
             dxcar_dxsph(4,4)=dzdxi

        !metric components of pure AdS in Cartesian coordinates
        g0_tt_ads0 =-(4*rho0**2+L**2*(-1+rho0**2)**2)
     &               /L**2/(-1+rho0**2)**2
        g0_tx_ads0 =0
        g0_ty_ads0 =0
        g0_tz_ads0 =0
        g0_xx_ads0 =(8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &              +8*rho0**2+4*L**2*(1+rho0**4))
     &              /(-1+rho0**2)**2/(4*rho0**2+L**2*(-1+rho0**2)**2)
        g0_xy_ads0 =(16 *(-1 + L**2) *x0* y0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_xz_ads0 =(16 *(-1 + L**2) *x0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_yy_ads0 =(4*(4*(x0**2+z0**2)+L**2*(x0**4+(1+y0**2)**2
     &              +2*(-1+y0**2)*z0**2+z0**4
     &              +2*x0**2*(-1+y0**2+z0**2))))
     &              /(L**2*(-1+rho0**2)**4+4*(-1+rho0**2)**2*(rho0**2))
        g0_yz_ads0 =(16 *(-1 + L**2) *y0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_psi_ads0=(4*(4*(x0**2+y0**2)+L**2*((-1+x0**2+y0**2)**2
     &              +2*(1+x0**2+y0**2)*z0**2+z0**4)))
     &              /(L**2*(-1+rho0**2)**4
     &              +4*(-1+rho0**2)**2*(rho0**2))

      !deviation from pure AdS in Cartesian coordinates

        hcar_n(1,1)=gb_tt_n(i,j,k)
        hcar_n(1,2)=gb_tx_n(i,j,k)
        hcar_n(1,3)=gb_ty_n(i,j,k)
        hcar_n(1,4)=gb_tz_n(i,j,k)
        hcar_n(2,2)=gb_xx_n(i,j,k)
        hcar_n(2,3)=gb_xy_n(i,j,k)
        hcar_n(2,4)=gb_xz_n(i,j,k)
        hcar_n(3,3)=gb_yy_n(i,j,k)
        hcar_n(3,4)=gb_yz_n(i,j,k)
        hcar_n(4,4)=psi_n(i,j,k)

       do a=1,3
          do b=a+1,4
            hcar_n(b,a)=hcar_n(a,b)
          end do
        end do

       !full metric in Cartesian coordinates
        g0car_n(1,1)=g0_tt_ads0+hcar_n(1,1)
        g0car_n(1,2)=g0_tx_ads0+hcar_n(1,2)
        g0car_n(1,3)=g0_ty_ads0+hcar_n(1,3)
        g0car_n(1,4)=g0_tz_ads0+hcar_n(1,4)
        g0car_n(2,2)=g0_xx_ads0+hcar_n(2,2)
        g0car_n(2,3)=g0_xy_ads0+hcar_n(2,3)
        g0car_n(2,4)=g0_xz_ads0+hcar_n(2,4)
        g0car_n(3,3)=g0_yy_ads0+hcar_n(3,3)
        g0car_n(3,4)=g0_yz_ads0+hcar_n(3,4)
        g0car_n(4,4)=g0_psi_ads0+hcar_n(4,4)

       do a=1,3
          do b=a+1,4
            g0car_n(b,a)=g0car_n(a,b)
          end do
        end do


       !deviation from pure AdS in (rescaled) spherical coordinates
        do a=1,4
          do b=1,4
           hsph_n(a,b)=0.0d0
           do c=1,4
            do d=1,4
             hsph_n(a,b)=hsph_n(a,b)
     &                   +dxcar_dxsph(c,a)*dxcar_dxsph(d,b)*hcar_n(c,d)
            end do
           end do
          end do
        end do


        !regularised metric components in (rescaled) spherical coordinates
        gbsph_n(1,1)=hsph_n(1,1)
        gbsph_n(1,2)=hsph_n(1,2)/(1-rho0**2)
        gbsph_n(1,3)=hsph_n(1,3)
        gbsph_n(1,4)=hsph_n(1,4)
        gbsph_n(2,2)=hsph_n(2,2)
        gbsph_n(2,3)=hsph_n(2,3)/(1-rho0**2)
        gbsph_n(2,4)=hsph_n(2,4)/(1-rho0**2)
        gbsph_n(3,3)=hsph_n(3,3)
        gbsph_n(3,4)=hsph_n(3,4)
        gbsph_n(4,4)=hsph_n(4,4)

         gbsph_tt_n(i,j,k)    =gbsph_n(1,1)
         gbsph_trho_n(i,j,k)  =gbsph_n(1,2)
         gbsph_tchi_n(i,j,k)  =gbsph_n(1,3)
         gbsph_txi_n(i,j,k)   =gbsph_n(1,4)
         gbsph_rhorho_n(i,j,k)=gbsph_n(2,2)
         gbsph_rhochi_n(i,j,k)=gbsph_n(2,3)
         gbsph_rhoxi_n(i,j,k) =gbsph_n(2,4)
         gbsph_chichi_n(i,j,k)=gbsph_n(3,3)
         gbsph_chixi_n(i,j,k) =gbsph_n(3,4)
         gbsph_xixi_n(i,j,k)  =gbsph_n(4,4)

        do a=1,3
          do b=a+1,4
            gbsph_n(b,a)=gbsph_n(a,b)
          end do
        end do

        !full metric in (rescaled) spherical coordinates
        do a=1,4
          do b=1,4
           g0sph_n(a,b)=0.0d0
           do c=1,4
            do d=1,4
             g0sph_n(a,b)=g0sph_n(a,b)
     &                   +dxcar_dxsph(c,a)*dxcar_dxsph(d,b)*g0car_n(c,d)
            end do
           end do
          end do
        end do

        !induced metric on hypersurface at constant rho
        gamma0sph_ll(1,1)=g0sph_n(1,1)
        gamma0sph_ll(1,2)=0
        gamma0sph_ll(1,3)=g0sph_n(1,3)
        gamma0sph_ll(1,4)=g0sph_n(1,4)
        gamma0sph_ll(2,2)=0
        gamma0sph_ll(2,3)=0
        gamma0sph_ll(2,4)=0
        gamma0sph_ll(3,3)=g0sph_n(3,3)
        gamma0sph_ll(3,4)=g0sph_n(3,4)
        gamma0sph_ll(4,4)=g0sph_n(4,4)

        do a=1,3
          do b=a+1,4
            gamma0sph_ll(b,a)=gamma0sph_ll(a,b)
          end do
        end do

        !metric on conformal AdS boundary at rho=1
        gamma0sphbdy_ll(1,1)=g0sph_n(1,1)*q**2
        gamma0sphbdy_ll(1,2)=g0sph_n(1,3)*q**2
        gamma0sphbdy_ll(1,3)=g0sph_n(1,4)*q**2
        gamma0sphbdy_ll(2,2)=g0sph_n(3,3)*q**2
        gamma0sphbdy_ll(2,3)=g0sph_n(3,4)*q**2
        gamma0sphbdy_ll(3,3)=g0sph_n(4,4)*q**2

        do a=1,2
          do b=a+1,3
            gamma0sphbdy_ll(b,a)=gamma0sphbdy_ll(a,b)
          end do
        end do

!
!!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!!       write (*,*) 'gamma0sph_ll(1,1)=',gamma0sph_ll(1,1)
!!       write (*,*) 'gamma0sph_ll(1,2)=',gamma0sph_ll(1,2)
!!       write (*,*) 'gamma0sph_ll(1,3)=',gamma0sph_ll(1,3)
!!       write (*,*) 'gamma0sph_ll(1,4)=',gamma0sph_ll(1,4)
!!       write (*,*) 'gamma0sph_ll(2,2)=',gamma0sph_ll(2,2)
!!       write (*,*) 'gamma0sph_ll(2,3)=',gamma0sph_ll(2,3)
!!       write (*,*) 'gamma0sph_ll(2,4)=',gamma0sph_ll(2,4)
!!       write (*,*) 'gamma0sph_ll(3,3)=',gamma0sph_ll(3,3)
!!       write (*,*) 'gamma0sph_ll(3,4)=',gamma0sph_ll(3,4)
!!       write (*,*) 'gamma0sph_ll(4,4)=',gamma0sph_ll(4,4)
!
      !determinant of induced metric on hypersurface at constant rho in 3-D form
        detgamma3=-g0sph_n(1,4)**2*g0sph_n(3,3)
     &            +2*g0sph_n(1,3)*g0sph_n(1,4)*g0sph_n(3,4)
     &            -g0sph_n(1,3)**2*g0sph_n(4,4)
     &            +g0sph_n(1,1)*(-g0sph_n(3,4)**2
     &              +g0sph_n(3,3)*g0sph_n(4,4))

      !induced inverse metric on hypersurface at constant rho
        gamma0sph_uu(1,1)=-(g0sph_n(3,4)**2-g0sph_n(3,3)*g0sph_n(4,4))
     &                    /detgamma3
        gamma0sph_uu(1,2)=0
        gamma0sph_uu(1,3)=-(-g0sph_n(1,4)*g0sph_n(3,4)
     &                      +g0sph_n(1,3)*g0sph_n(4,4))
     &                      /detgamma3
        gamma0sph_uu(1,4)=-(g0sph_n(1,4)*g0sph_n(3,3)
     &                      -g0sph_n(1,3)*g0sph_n(3,4))
     &                      /detgamma3
        gamma0sph_uu(2,2)=0
        gamma0sph_uu(2,3)=0
        gamma0sph_uu(2,4)=0
        gamma0sph_uu(3,3)=-(g0sph_n(1,4)**2-g0sph_n(1,1)*g0sph_n(4,4))
     &                    /detgamma3
        gamma0sph_uu(3,4)=-(-g0sph_n(1,3)*g0sph_n(1,4)
     &                      +g0sph_n(1,1)*g0sph_n(3,4))
     &                      /detgamma3
        gamma0sph_uu(4,4)=-(g0sph_n(1,3)**2-g0sph_n(1,1)*g0sph_n(3,3))
     &                    /detgamma3

         do a=1,3
           do b=a+1,4
             gamma0sph_uu(b,a)=gamma0sph_uu(a,b)
           end do
         end do

!!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!!       write (*,*) 'gamma0sph_uu(1,1)=',gamma0sph_uu(1,1)
!!       write (*,*) 'gamma0sph_uu(1,2)=',gamma0sph_uu(1,2)
!!       write (*,*) 'gamma0sph_uu(1,3)=',gamma0sph_uu(1,3)
!!       write (*,*) 'gamma0sph_uu(1,4)=',gamma0sph_uu(1,4)
!!       write (*,*) 'gamma0sph_uu(2,2)=',gamma0sph_uu(2,2) 
!!       write (*,*) 'gamma0sph_uu(2,3)=',gamma0sph_uu(2,3) 
!!       write (*,*) 'gamma0sph_uu(2,4)=',gamma0sph_uu(2,4) 
!!       write (*,*) 'gamma0sph_uu(3,3)=',gamma0sph_uu(3,3) 
!!       write (*,*) 'gamma0sph_uu(3,4)=',gamma0sph_uu(3,4) 
!!       write (*,*) 'gamma0sph_uu(4,4)=',gamma0sph_uu(4,4)
!!

!     !determinant of metric on conformal AdS boundary at rho=1
       detgamma0sphbdy=
     &            -gamma0sphbdy_ll(1,3)**2*gamma0sphbdy_ll(2,2)
     &            +2*gamma0sphbdy_ll(1,2)*gamma0sphbdy_ll(1,3)
     &                *gamma0sphbdy_ll(2,3)
     &            -gamma0sphbdy_ll(1,2)**2*gamma0sphbdy_ll(3,3)
     &            +gamma0sphbdy_ll(1,1)*(-gamma0sphbdy_ll(2,3)**2
     &              +gamma0sphbdy_ll(2,2)*gamma0sphbdy_ll(3,3))

!     !inverse of metric on conformal AdS boundary at rho=1
       gamma0sphbdy_uu(1,1)=-(gamma0sphbdy_ll(2,3)**2
     &                       -gamma0sphbdy_ll(2,2)*gamma0sphbdy_ll(3,3))
     &                     /detgamma0sphbdy
       gamma0sphbdy_uu(1,2)=
     &                      -(-gamma0sphbdy_ll(1,3)*gamma0sphbdy_ll(2,3)
     &                      +gamma0sphbdy_ll(1,2)*gamma0sphbdy_ll(3,3))
     &                      /detgamma0sphbdy
       gamma0sphbdy_uu(1,3)=
     &                      -(gamma0sphbdy_ll(1,3)*gamma0sphbdy_ll(2,2)
     &                      -gamma0sphbdy_ll(1,2)*gamma0sphbdy_ll(2,3))
     &                      /detgamma0sphbdy
       gamma0sphbdy_uu(2,2)=-(gamma0sphbdy_ll(1,3)**2
     &                       -gamma0sphbdy_ll(1,1)*gamma0sphbdy_ll(3,3))
     &                    /detgamma0sphbdy
       gamma0sphbdy_uu(2,3)=
     &                      -(-gamma0sphbdy_ll(1,2)*gamma0sphbdy_ll(1,3)
     &                      +gamma0sphbdy_ll(1,1)*gamma0sphbdy_ll(2,3))
     &                      /detgamma0sphbdy
       gamma0sphbdy_uu(3,3)=
     &                       -(gamma0sphbdy_ll(1,2)**2
     &                       -gamma0sphbdy_ll(1,1)*gamma0sphbdy_ll(2,2))
     &                    /detgamma0sphbdy

         do a=1,2
           do b=a+1,3
             gamma0sphbdy_uu(b,a)=gamma0sphbdy_uu(a,b)
           end do
         end do

          gamma0sphbdy_uu_tt(i,j,k)=gamma0sphbdy_uu(1,1)
          gamma0sphbdy_uu_tchi(i,j,k)=gamma0sphbdy_uu(1,2)
          gamma0sphbdy_uu_txi(i,j,k)=gamma0sphbdy_uu(1,3)
          gamma0sphbdy_uu_chichi(i,j,k)=gamma0sphbdy_uu(2,2)
          gamma0sphbdy_uu_chixi(i,j,k)=gamma0sphbdy_uu(2,3)
          gamma0sphbdy_uu_xixi(i,j,k)=gamma0sphbdy_uu(3,3)

! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
! identified with the bdy CFT stress-energy tensor one-point function
!these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor

!initialise
              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0
              quasiset_tracell(i,j,k)=0

            if (no_derivatives) then
!             if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
!
               quasiset_tt_ll(i,j,k)=(12*(gbsph_n(3,3)/q)
     &                         + 8*PI**2*(gbsph_n(2,2)/q)
     &                         +  (3*(gbsph_n(4,4)/q)/(sin(PI*chi0))**2)
     &                         )/(64*PI**3)


               quasiset_tchi_ll(i,j,k)  = (3*(gbsph_n(1,3)/q))/(16*PI)


               quasiset_txi_ll(i,j,k)   = (3*(gbsph_n(1,4)/q))/(16*PI)

               quasiset_chichi_ll(i,j,k)=(3.0d0/16.0d0)*PI
     &                                   *(gbsph_n(1,1)/q)
     &                                  -(1.0d0/8.0d0)*PI
     &                                   *(gbsph_n(2,2)/q)
     &                           -(3*(gbsph_n(4,4)/q)/(sin(PI*chi0))**2)
     &                                   /(64*PI)


               quasiset_chixi_ll(i,j,k) =(3*(gbsph_n(3,4)/q))/(16*PI)


               quasiset_xixi_ll(i,j,k)  =((sin(PI*chi0))**2*(-3
     &                                   *(gbsph_n(3,3)/q)
     &                                  +PI**2*(3*(gbsph_n(1,1)/q)
     &                                  -2*(gbsph_n(2,2)/q)))
     &                                  )/(4*PI)

               quasiset_massdensityll(i,j,k)=(sin(PI*chi0))
     &                                    *(
     &                                     12*(gbsph_n(3,3)/q)
     &                                    +8*PI**2*(gbsph_n(2,2)/q)
     &                                    +3*(gbsph_n(4,4)/q)
     &                                     /((sin(PI*chi0))**2)
     &                                    )
     &                                    /(32*PI)


       !trace of quasi local stress-tensor from definition of trace
             quasiset_tracell(i,j,k)=(
     &           gamma0sphbdy_uu_tt(i,j,k)*quasiset_tt_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_tchi(i,j,k)*quasiset_tchi_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_txi(i,j,k)*quasiset_txi_ll(i,j,k)
     &          +gamma0sphbdy_uu_chichi(i,j,k)*quasiset_chichi_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_chixi(i,j,k)*quasiset_chixi_ll(i,j,k)
     &          +gamma0sphbdy_uu_xixi(i,j,k)*quasiset_xixi_ll(i,j,k)
     &                   )
!
!!        write(*,*) "i,j,k,x(i),y(j),z(k),rho0="
!!     &             ,i,j,k,x(i),y(j),z(k),rho0
!!        write(*,*) "TRACE: quasiset_tracell(i,j,k)="
!!     &             ,quasiset_tracell(i,j,k)
!!
!!!        write(*,*) "i,j,k,x(i),y(j),z(k),rho0="
!!!     &             ,i,j,k,x(i),y(j),z(k),rho0
!!!        write(*,*) "TRACE: quasiset_tracell(i,j,k)="
!!!     &             ,quasiset_tracell(i,j,k)
!!!        write(*,*) "TRACE: quasiset_tt_ll(i,j,k)="
!!!     &             ,quasiset_tt_ll(i,j,k) 
!!!
!!!
!!!       !trace of quasi local stress-tensor in terms of regularised metric components (this should be the same as the one above, within numerical error)
!!!              quasiset_tracell(i,j,k)=
!!!     &                               -(3/(8*PI))*(-(gbsph_n(1,1)/q)
!!!     &                               +(gbsph_n(2,2)/q)
!!!     &                               +(gbsph_n(3,3)/q)
!!!     &                               +(gbsph_n(4,4)/q)/(sin(PI*chi0)**2)
!!!     &                               )
!!!
!!!        write(*,*) "EXPANSION: quasiset_tracell(i,j,k)="
!!!     &             ,quasiset_tracell(i,j,k)
!
!!!!!!!!!!case where y0=z0=0
!!             else
!!              quasiset_tt_ll(i,j,k)= (x0**2*(2*(gb_xx_n(i,j,k)/q)
!!     &                                +3*rho0**2*(psi_n(i,j,k)/q))
!!     &                                +3*rho0**4*(gb_yy_n(i,j,k)/q) )
!!     &                               /(16*PI*rho0**2)
!!              quasiset_tchi_ll(i,j,k)  =(3/16)*x0*(gb_tz_n(i,j,k)/q)
!!              quasiset_txi_ll(i,j,k)   =0
!!           quasiset_chichi_ll(i,j,k)=-PI*(-3*rho0**2*(gb_tt_n(i,j,k)/q)
!!     &                                   +2*x0**2*(gb_xx_n(i,j,k)/q)
!!     &                                   +3*rho0**4*(gb_yy_n(i,j,k)/q))
!!     &                                   /(16*rho0**2)
!!              quasiset_chixi_ll(i,j,k) =0
!!              quasiset_xixi_ll(i,j,k)  =0
!!              quasiset_massdensityll(i,j,k)=0
!!              quasiset_tracell(i,j,k)=-((3*(-rho0**2*(gb_tt_n(i,j,k)/q)
!!     &                                +x0**2*((gb_xx_n(i,j,k)/q)
!!     &                                  +PI**2*rho0**2*(psi_n(i,j,k)/q))
!!     &                                  +4*(PI**2)*(rho0**4)
!!     &                                   *(gb_yy_n(i,j,k)/q)))
!!     &                                 /(8*PI*rho0**2))
!!             end if
            end if

          else  !excised points or points with y0=z0=0

              gbsph_tt_n(i,j,k)    =0
              gbsph_trho_n(i,j,k)  =0
              gbsph_tchi_n(i,j,k)  =0
              gbsph_txi_n(i,j,k)   =0
              gbsph_rhorho_n(i,j,k)=0
              gbsph_rhochi_n(i,j,k)=0
              gbsph_rhoxi_n(i,j,k) =0
              gbsph_chichi_n(i,j,k)  =0
              gbsph_chixi_n(i,j,k) =0
              gbsph_xixi_n(i,j,k)=0

              gamma0sphbdy_uu_tt(i,j,k)=0
              gamma0sphbdy_uu_tchi(i,j,k)=0
              gamma0sphbdy_uu_txi(i,j,k)=0
              gamma0sphbdy_uu_chichi(i,j,k)=0
              gamma0sphbdy_uu_chixi(i,j,k)=0
              gamma0sphbdy_uu_xixi(i,j,k)=0

              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0
              quasiset_tracell(i,j,k)=0

          end if


         end do
        end do
       end do


      if (.not.no_derivatives) then !we need to restric the range of definition to points owned only by the current process (not any other process) if we want to use derivatives
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1


!!!ghost_width not needed as long as we don't use derivatives
        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        do i=is,ie
         do j=js,je
          do k=ks,ke

             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if


!                ! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
!                ! identified with the bdy CFT stress-energy tensor one-point function
!                !these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor
            if ( (chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then
              dgbsph_tt_drho_n    =
     &             df_drho(gbsph_tt_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_trho_drho_n  =
     &             df_drho(gbsph_trho_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_tchi_drho_n  =
     &             df_drho(gbsph_tchi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_txi_drho_n   =
     &             df_drho(gbsph_txi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_rhorho_drho_n=
     &             df_drho(gbsph_rhorho_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_rhochi_drho_n=
     &             df_drho(gbsph_rhochi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_rhoxi_drho_n =
     &             df_drho(gbsph_rhoxi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_chichi_drho_n=
     &             df_drho(gbsph_chichi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_chixi_drho_n =
     &             df_drho(gbsph_chixi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              dgbsph_xixi_drho_n  =
     &             df_drho(gbsph_xixi_n,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

               quasiset_tt_ll(i,j,k)=
     &                         (12*(-dgbsph_chichi_drho_n)
     &                         + 8*PI**2*(-dgbsph_rhorho_drho_n)
     &                         +  (3*(-dgbsph_xixi_drho_n)
     &                               /(sin(PI*chi0))**2)
     &                         )/(64*PI**3)


               quasiset_tchi_ll(i,j,k)  =
     &                        (3*(-dgbsph_tchi_drho_n))
     &                                       /(16*PI)


               quasiset_txi_ll(i,j,k)   =
     &                        (3*(-dgbsph_txi_drho_n))
     &                                       /(16*PI)


               quasiset_chichi_ll(i,j,k)=
     &                                   (3.0d0/16.0d0)*PI
     &                                   *(-dgbsph_tt_drho_n)
     &                                  -(1.0d0/8.0d0)*PI
     &                                   *(-dgbsph_rhorho_drho_n)
     &                           -(3*(-dgbsph_xixi_drho_n)
     &                               /(sin(PI*chi0))**2)
     &                                   /(64*PI)


               quasiset_chixi_ll(i,j,k) =
     &                                   (3*(-dgbsph_chixi_drho_n))
     &                                      /(16*PI)



               quasiset_xixi_ll(i,j,k)=
     &                                  ((sin(PI*chi0))**2*(-3
     &                                   *(-dgbsph_chichi_drho_n)
     &                                  +PI**2*(3*(-dgbsph_tt_drho_n)
     &                                  -2*(-dgbsph_rhorho_drho_n)))
     &                                  )/(4*PI)

               quasiset_massdensityll(i,j,k)=
     &                                     (sin(PI*chi0))
     &                                    *(
     &                                     12*(-dgbsph_chichi_drho_n)
     &                                    +8*PI**2
     &                                        *(-dgbsph_rhorho_drho_n)
     &                                    +3*(-dgbsph_xixi_drho_n)
     &                                     /((sin(PI*chi0))**2)
     &                                    )
     &                                    /(32*PI)


       !trace of quasi local stress-tensor from definition of trace
             quasiset_tracell(i,j,k)=
     &                               (
     &           gamma0sphbdy_uu_tt(i,j,k)*quasiset_tt_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_tchi(i,j,k)*quasiset_tchi_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_txi(i,j,k)*quasiset_txi_ll(i,j,k)
     &          +gamma0sphbdy_uu_chichi(i,j,k)*quasiset_chichi_ll(i,j,k)
     &        +2*gamma0sphbdy_uu_chixi(i,j,k)*quasiset_chixi_ll(i,j,k)
     &          +gamma0sphbdy_uu_xixi(i,j,k)*quasiset_xixi_ll(i,j,k)
     &                   )


           else !excised points or points with y0=z0=0

              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0
              quasiset_tracell(i,j,k)=0
            end if

          end do
         end do
        end do

       end if

       return
       end
c--------------------------------------------------------------------------------------



c-------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 AT the boundary through extrapolation.
c The tensor components are given in spherical polar coordinates.
c-------------------------------------------------------------------------------------

        subroutine quasiset(
     &                  quasiset_tt,quasiset_tchi,quasiset_txi,
     &                  quasiset_chichi,quasiset_chixi,
     &                  quasiset_xixi,
     &                  quasiset_massdensity,
     &                  quasiset_trace,
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_massdensityll,
     &                  quasiset_tracell,
     &                  xextrap,yextrap,zextrap,
     &                  chrbdy,numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!!!!!!!!define variables!!!!!!!!!

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        real*8 chrbdy(Nx,Ny,Nz)
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 psi_np1(Nx,Ny,Nz),psi_n(Nx,Ny,Nz),psi_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex
        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)
        real*8 quasiset_tracell(Nx,Ny,Nz)

        real*8 quasiset_tt_p1,quasiset_tchi_p1
        real*8 quasiset_txi_p1,quasiset_chichi_p1
        real*8 quasiset_chixi_p1,quasiset_xixi_p1
        real*8 quasiset_massdensity_p1
        real*8 quasiset_trace_p1

        real*8 quasiset_tt_p2,quasiset_tchi_p2
        real*8 quasiset_txi_p2,quasiset_chichi_p2
        real*8 quasiset_chixi_p2,quasiset_xixi_p2
        real*8 quasiset_massdensity_p2
        real*8 quasiset_trace_p2

        real*8 quasiset_tt(numbdypoints),quasiset_tchi(numbdypoints)
        real*8 quasiset_txi(numbdypoints),quasiset_chichi(numbdypoints)
        real*8 quasiset_chixi(numbdypoints),quasiset_xixi(numbdypoints)
        real*8 quasiset_massdensity(numbdypoints)
        real*8 quasiset_trace(numbdypoints)

        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 rhoextrap
        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        integer bdy_Nchi,bdy_Nxi

!        real*8 chibdy(bdy_Nchi)
!        real*8 xibdy(bdy_Nxi)

        integer i,j,k,is,ie,js,je,ks,ke,lind,m,e,increase

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 xp1,yp1,zp1,rhop1,xp2,yp2,zp2
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 maxxyzp1

        real*8 extrapalongx,extrapalongy,extrapalongz

        real*8 dp1p2

        real*8 AdS_mass

!!!!!!!!!!!!!!!!

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)


        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)
       
 
        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)

           quasiset_tt_p1=quasiset_tt_ll(i,j,k)
           quasiset_tchi_p1=quasiset_tchi_ll(i,j,k)
           quasiset_txi_p1=quasiset_txi_ll(i,j,k)
           quasiset_chichi_p1=quasiset_chichi_ll(i,j,k)
           quasiset_chixi_p1=quasiset_chixi_ll(i,j,k)
           quasiset_xixi_p1=quasiset_xixi_ll(i,j,k)
           quasiset_massdensity_p1=quasiset_massdensityll(i,j,k)
           quasiset_trace_p1=quasiset_tracell(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then

              lind=lind+1

                  xex=xextrap(lind)
                  yex=yextrap(lind)
                  zex=zextrap(lind)
                  rhoex=sqrt(xex**2+yex**2+zex**2)
                  chiex=(1/PI)*acos(xex/rhoex)
                  if (zex.lt.0) then
                    xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                  else
                    xiex=(1/(2*PI))*atan2(zex,yex)
                  end if


              if (maxxyzp1.eq.abs(xp1)) then
                if (xp1.gt.0) then

                  xp2=x(i-1)
                  quasiset_tt_p2=quasiset_tt_ll(i-1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i-1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i-1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i-1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i-1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i-1,j,k)
                  quasiset_trace_p2=quasiset_tracell(i-1,j,k)
                  quasiset_massdensity_p2=
     &                           quasiset_massdensityll(i-1,j,k)


                  quasiset_tt(lind)=extrapalongx(quasiset_tt_p1
     &                         ,quasiset_tt_p2,xp1,xp2,xex)
                  quasiset_tchi(lind)=extrapalongx(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,xp1,xp2,xex)
                  quasiset_txi(lind)=extrapalongx(quasiset_txi_p1
     &                         ,quasiset_txi_p2,xp1,xp2,xex)
                  quasiset_chichi(lind)=extrapalongx(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,xp1,xp2,xex)
                  quasiset_chixi(lind)=extrapalongx(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,xp1,xp2,xex)
                  quasiset_xixi(lind)=extrapalongx(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,xp1,xp2,xex)
          quasiset_trace(lind)=
     &                         extrapalongx(quasiset_trace_p1
     &                         ,quasiset_trace_p2,xp1,xp2,xex)
          quasiset_massdensity(lind)=
     &                         extrapalongx(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,xp1,xp2,xex)
 
 
              else
                  xp2=x(i+1)
                  quasiset_tt_p2=quasiset_tt_ll(i+1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i+1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i+1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i+1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i+1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i+1,j,k)
                  quasiset_trace_p2=quasiset_tracell(i+1,j,k)
                  quasiset_massdensity_p2=
     &                                 quasiset_massdensityll(i+1,j,k)


                  quasiset_tt(lind)=extrapalongx(quasiset_tt_p1
     &                         ,quasiset_tt_p2,xp1,xp2,xex)
                  quasiset_tchi(lind)=extrapalongx(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,xp1,xp2,xex)
                  quasiset_txi(lind)=extrapalongx(quasiset_txi_p1
     &                         ,quasiset_txi_p2,xp1,xp2,xex)
                  quasiset_chichi(lind)=extrapalongx(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,xp1,xp2,xex)
                  quasiset_chixi(lind)=extrapalongx(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,xp1,xp2,xex)
                  quasiset_xixi(lind)=extrapalongx(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,xp1,xp2,xex)
          quasiset_trace(lind)=
     &                         extrapalongx(quasiset_trace_p1
     &                         ,quasiset_trace_p2,xp1,xp2,xex)
          quasiset_massdensity(lind)=
     &                         extrapalongx(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,xp1,xp2,xex)
 
 
               end if
             else if (maxxyzp1.eq.abs(yp1)) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j-1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j-1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j-1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j-1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j-1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j-1,k)
                  quasiset_trace_p2=quasiset_tracell(i,j-1,k)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j-1,k)


                 quasiset_tt(lind)=extrapalongy(quasiset_tt_p1
     &                         ,quasiset_tt_p2,yp1,yp2,yex)
                  quasiset_tchi(lind)=extrapalongy(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,yp1,yp2,yex)
                  quasiset_txi(lind)=extrapalongy(quasiset_txi_p1
     &                         ,quasiset_txi_p2,yp1,yp2,yex)
                  quasiset_chichi(lind)=extrapalongy(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,yp1,yp2,yex)
                  quasiset_chixi(lind)=extrapalongy(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,yp1,yp2,yex)
                  quasiset_xixi(lind)=extrapalongy(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,yp1,yp2,yex)
          quasiset_trace(lind)=
     &                         extrapalongy(quasiset_trace_p1
     &                         ,quasiset_trace_p2,yp1,yp2,yex)
          quasiset_massdensity(lind)=
     &                         extrapalongy(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,yp1,yp2,yex)
 
 
              else
                  yp2=y(j+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j+1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j+1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j+1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j+1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j+1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j+1,k)
                  quasiset_trace_p2=quasiset_tracell(i,j+1,k)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j+1,k)


                 quasiset_tt(lind)=extrapalongy(quasiset_tt_p1
     &                         ,quasiset_tt_p2,yp1,yp2,yex)
                  quasiset_tchi(lind)=extrapalongy(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,yp1,yp2,yex)
                  quasiset_txi(lind)=extrapalongy(quasiset_txi_p1
     &                         ,quasiset_txi_p2,yp1,yp2,yex)
                  quasiset_chichi(lind)=extrapalongy(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,yp1,yp2,yex)
                  quasiset_chixi(lind)=extrapalongy(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,yp1,yp2,yex)
                  quasiset_xixi(lind)=extrapalongy(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,yp1,yp2,yex)
          quasiset_trace(lind)=
     &                         extrapalongy(quasiset_trace_p1
     &                         ,quasiset_trace_p2,yp1,yp2,yex)
          quasiset_massdensity(lind)=
     &                         extrapalongy(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,yp1,yp2,yex)
 
 
              end if             
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k-1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k-1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k-1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k-1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k-1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k-1)
                  quasiset_trace_p2=quasiset_tracell(i,j,k-1)
                quasiset_massdensity_p2=quasiset_massdensityll(i,j,k-1)


                 quasiset_tt(lind)=extrapalongz(quasiset_tt_p1
     &                         ,quasiset_tt_p2,zp1,zp2,zex)
                  quasiset_tchi(lind)=extrapalongz(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,zp1,zp2,zex)
                  quasiset_txi(lind)=extrapalongz(quasiset_txi_p1
     &                         ,quasiset_txi_p2,zp1,zp2,zex)
                  quasiset_chichi(lind)=extrapalongz(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,zp1,zp2,zex)
                  quasiset_chixi(lind)=extrapalongz(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,zp1,zp2,zex)
                  quasiset_xixi(lind)=extrapalongz(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,zp1,zp2,zex)
          quasiset_trace(lind)=
     &                         extrapalongz(quasiset_trace_p1
     &                         ,quasiset_trace_p2,zp1,zp2,zex)
          quasiset_massdensity(lind)=
     &                         extrapalongz(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,zp1,zp2,zex)
 
              else
                  zp2=z(k+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k+1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k+1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k+1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k+1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k+1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k+1)
                  quasiset_trace_p2=quasiset_tracell(i,j,k+1)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j,k+1)


                  quasiset_tt(lind)=extrapalongz(quasiset_tt_p1
     &                         ,quasiset_tt_p2,zp1,zp2,zex)
                  quasiset_tchi(lind)=extrapalongz(quasiset_tchi_p1
     &                         ,quasiset_tchi_p2,zp1,zp2,zex)
                  quasiset_txi(lind)=extrapalongz(quasiset_txi_p1
     &                         ,quasiset_txi_p2,zp1,zp2,zex)
                  quasiset_chichi(lind)=extrapalongz(quasiset_chichi_p1
     &                         ,quasiset_chichi_p2,zp1,zp2,zex)
                  quasiset_chixi(lind)=extrapalongz(quasiset_chixi_p1
     &                         ,quasiset_chixi_p2,zp1,zp2,zex)
                  quasiset_xixi(lind)=extrapalongz(quasiset_xixi_p1
     &                         ,quasiset_xixi_p2,zp1,zp2,zex)
          quasiset_trace(lind)=
     &                         extrapalongz(quasiset_trace_p1
     &                         ,quasiset_trace_p2,zp1,zp2,zex)
          quasiset_massdensity(lind)=
     &                         extrapalongz(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,zp1,zp2,zex)
 
               end if
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if
          end do
         end do
        end do

!        do lind=1,numbdypoints
!         write(*,*) "lind,xextrap(lind),yextrap(lind),zextrap(lind)="
!     &              ,lind,xextrap(lind),yextrap(lind),zextrap(lind)
!         write(*,*) "quasiset_trace(lind)=",quasiset_trace(lind)
!        end do

        return
        end
c--------------------------------------------------------------------------------------

c--------------------------------------------------------------------------
c derivative w.r.t. rho=sqrt(x**2+y**2+z**2)
c--------------------------------------------------------------------------
        real*8 function df_drho(f,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        implicit none

        integer Nx,Ny,Nz,i,j,k
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,x(Nx),y(Ny),z(Nz)
        real*8 f_x,f_y,f_z
        real*8 dxdrho,dydrho,dzdrho
        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 PI
        parameter (PI=3.141592653589793d0)

!--------------------------------------------------------------------

        call df1_int_x(f,f_x,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_y(f,f_y,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_z(f,f_z,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if

        dxdrho=cos(PI*chi0)
        dydrho=sin(PI*chi0)*cos(2*PI*xi0)
        dzdrho=sin(PI*chi0)*sin(2*PI*xi0)

        df_drho=dxdrho*f_x+dydrho*f_y+dzdrho*f_z

        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c 2-point extrapolation along x axis function using values Txp1,Txp2 at xp1,xp2
c----------------------------------------------------------------------
        real*8 function extrapalongx(Txp1,Txp2,xp1,xp2,xex)
        implicit none
        real*8 Txp1,Txp2,xp1,xp2,xex

        !--------------------------------------------------------------

        extrapalongx=Txp2+((Txp1-Txp2)/(xp1-xp2))*(xex-xp2)

        return 
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c 2-point extrapolation along y axis function using values Typ1,Typ2 at yp1,yp2
c----------------------------------------------------------------------
        real*8 function extrapalongy(Typ1,Typ2,yp1,yp2,yex)
        implicit none
        real*8 Typ1,Typ2,yp1,yp2,yex

        !--------------------------------------------------------------

        extrapalongy=Typ2+((Typ1-Typ2)/(yp1-yp2))*(yex-yp2)

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c 2-point extrapolation along z axis function using values Tzp1,Tzp2 at zp1,zp2
c----------------------------------------------------------------------
        real*8 function extrapalongz(Tzp1,Tzp2,zp1,zp2,zex)
        implicit none
        real*8 Tzp1,Tzp2,zp1,zp2,zex

        !--------------------------------------------------------------

        extrapalongz=Tzp2+((Tzp1-Tzp2)/(zp1-zp2))*(zex-zp2)

        return
        end
c--------------------------------------------------------------------------------------

c------------------------------------------------------------------------------------------------------------------------------
c The following routine performs the double integral on the S^2 boundary. We approximate the integrals via the trapezoidal rule
c------------------------------------------------------------------------------------------------------------------------------

        subroutine doubleintegralonsphere(integral,density,
     &                  xextrap,yextrap,zextrap,numbdypoints,
     &                  rhobdy,chibdy,xibdy,
     &                  bdy_Nchi,bdy_Nxi)

        implicit none
        integer numbdypoints
        real*8 density(numbdypoints)
        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 integral

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 rhobdy
        real*8 chibdy(bdy_Nchi),xibdy(bdy_Nxi)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e,i,j,lind
        integer lind_chipxip,lind_chipxipp1
        integer lind_chipp1xip,lind_chipp1xipp1
        integer additions

        real*8 rhoextrap,rho2

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        real*8 x_chipxip,y_chipxip,z_chipxip
        real*8 x_chipxipp1,y_chipxipp1,z_chipxipp1
        real*8 x_chipp1xip,y_chipp1xip,z_chipp1xip
        real*8 x_chipp1xipp1,y_chipp1xipp1,z_chipp1xipp1

        integer trig_chipxip,trig_chipxipp1
        integer trig_chipp1xip,trig_chipp1xipp1

        real*8 dist_extr_bdypp_min,dist_extr_bdyppp1_min
        real*8 dist_extr_bdypp1p_min,dist_extr_bdypp1pp1_min

        real*8 dist_extr_bdypp_a,dist_extr_bdyppp1_a
        real*8 dist_extr_bdypp1p_a,dist_extr_bdypp1pp1_a


        !-----------------------------------------------------------------------

        integral=0.0d0
        do i=1,bdy_Nchi-1
         do j=1,bdy_Nxi-1
           x_chipxip=rhobdy*cos(PI*chibdy(i))
           y_chipxip=rhobdy*sin(PI*chibdy(i))*cos(2*PI*xibdy(j))
           z_chipxip=rhobdy*sin(PI*chibdy(i))*sin(2*PI*xibdy(j))

           x_chipxipp1=rhobdy*cos(PI*chibdy(i))
           y_chipxipp1=rhobdy*sin(PI*chibdy(i))*cos(2*PI*xibdy(j+1))
           z_chipxipp1=rhobdy*sin(PI*chibdy(i))*sin(2*PI*xibdy(j+1))

           x_chipp1xip=rhobdy*cos(PI*chibdy(i+1))
           y_chipp1xip=rhobdy*sin(PI*chibdy(i+1))*cos(2*PI*xibdy(j))
           z_chipp1xip=rhobdy*sin(PI*chibdy(i+1))*sin(2*PI*xibdy(j))

           x_chipp1xipp1=rhobdy*cos(PI*chibdy(i+1))
           y_chipp1xipp1=rhobdy*sin(PI*chibdy(i+1))*cos(2*PI*xibdy(j+1))
           z_chipp1xipp1=rhobdy*sin(PI*chibdy(i+1))*sin(2*PI*xibdy(j+1))

             dist_extr_bdypp_min=sqrt((xextrap(1)-x_chipxip)**2
     &                              +(yextrap(1)-y_chipxip)**2
     &                              +(zextrap(1)-z_chipxip)**2)
             lind_chipxip=1

             dist_extr_bdyppp1_min=sqrt((xextrap(1)-x_chipxipp1)**2
     &                              +(yextrap(1)-y_chipxipp1)**2
     &                              +(zextrap(1)-z_chipxipp1)**2)
             lind_chipxipp1=1

             dist_extr_bdypp1p_min=sqrt((xextrap(1)-x_chipp1xip)**2
     &                              +(yextrap(1)-y_chipp1xip)**2
     &                              +(zextrap(1)-z_chipp1xip)**2)
             lind_chipp1xip=1

             dist_extr_bdypp1pp1_min=sqrt((xextrap(1)-x_chipp1xipp1)**2
     &                              +(yextrap(1)-y_chipp1xipp1)**2
     &                              +(zextrap(1)-z_chipp1xipp1)**2)
             lind_chipp1xipp1=1

           do lind=2,numbdypoints

              dist_extr_bdypp_a=sqrt((xextrap(lind)-x_chipxip)**2
     &                              +(yextrap(lind)-y_chipxip)**2
     &                              +(zextrap(lind)-z_chipxip)**2)
              if (dist_extr_bdypp_a.lt.dist_extr_bdypp_min) then
                  dist_extr_bdypp_min=dist_extr_bdypp_a
                  lind_chipxip=lind
              end if

              dist_extr_bdyppp1_a=sqrt((xextrap(lind)-x_chipxipp1)**2
     &                              +(yextrap(lind)-y_chipxipp1)**2
     &                              +(zextrap(lind)-z_chipxipp1)**2)
              if (dist_extr_bdyppp1_a.lt.dist_extr_bdyppp1_min) then
                  dist_extr_bdyppp1_min=dist_extr_bdyppp1_a
                  lind_chipxipp1=lind
              end if

              dist_extr_bdypp1p_a=sqrt((xextrap(lind)-x_chipp1xip)**2
     &                              +(yextrap(lind)-y_chipp1xip)**2
     &                              +(zextrap(lind)-z_chipp1xip)**2)
              if (dist_extr_bdypp1p_a.lt.dist_extr_bdypp1p_min) then
                  dist_extr_bdypp1p_min=dist_extr_bdypp1p_a
                   lind_chipp1xip=lind
              end if

             dist_extr_bdypp1pp1_a=sqrt((xextrap(lind)-x_chipp1xipp1)**2
     &                              +(yextrap(lind)-y_chipp1xipp1)**2
     &                              +(zextrap(lind)-z_chipp1xipp1)**2)
              if (dist_extr_bdypp1pp1_a.lt.dist_extr_bdypp1pp1_min) then
                  dist_extr_bdypp1pp1_min=dist_extr_bdypp1pp1_a
                   lind_chipp1xipp1=lind
              end if
            end do

              integral=integral+
     &              (chibdy(i+1)-chibdy(i))/2 * (xibdy(j+1)-xibdy(j))/2
     &              *(density(lind_chipxip)+density(lind_chipxipp1)
     &              +density(lind_chipp1xip)+density(lind_chipp1xipp1))

         end do
        end do
 
        return
        end
!----------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c Routine to calculate the value of angular coordinates chi and xi at the boundary points where we extrapolate the quasi-local stress energy tensor
c--------------------------------------------------------------------------------------

        subroutine chixiextrap(
     &                         rhoextrap,chiextrap,xiextrap,
     &                         xextrap,yextrap,zextrap,
     &                         numbdypoints)

        integer numbdypoints

        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 rhoextrap
        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e

                !--------------------------------------------------------------

        rhoextrap=1.0d0
        do e=1,numbdypoints
         chiextrap(e)=(1/PI)*acos(xextrap(e)/rhoextrap)
         if (zextrap(e).lt.0) then
             xiextrap(e)=(1/(2*PI))*(atan2(zextrap(e),yextrap(e))+2*PI)
         else
             xiextrap(e)=(1/(2*PI))*atan2(zextrap(e),yextrap(e))
         end if
        end do

        return
        end
!----------------------------------------------------------------------


c---------------------------------------------------------------------------------------------------------------------------------------------------------
c The following routine calculate the number of different values taken by chi and xi at the boundary points where the stress-energy tensor is extrapolated
c---------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine bdyn(
     &                  bdy_Nchi,bdy_Nxi,
     &                  numbdypoints,
     &                  chiextrap,xiextrap)

        implicit none
        integer numbdypoints
        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e

        real*8 rhoextrap

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        !-----------------------------------------------------------------------

        chiextrap_min=minval(chiextrap)-1
        chiextrap_max=maxval(chiextrap)
        xiextrap_min=minval(xiextrap)-1
        xiextrap_max=maxval(xiextrap)

!calculate the number of different values of the chi and xi coordinates taken on extrapolated bdy points
        bdy_Nchi=0
        do while (chiextrap_min.lt.chiextrap_max)
            bdy_Nchi=bdy_Nchi+1
            chiextrap_min=
     &                minval(chiextrap,mask=chiextrap.gt.chiextrap_min)
        end do

        bdy_Nxi=0
        do while (xiextrap_min.lt.xiextrap_max)
            bdy_Nxi=bdy_Nxi+1
            xiextrap_min=minval(xiextrap,mask=xiextrap.gt.xiextrap_min)
        end do

        return
        end
!----------------------------------------------------------------------

c---------------------------------------------------------------------------------------------------------------------------------------------------------
c The following routine returns the arrays chibdy(bdy_Nchi) and xibdy(bdy_Nxi) containing the values of chi and xi at the extrapolated boundary points, in increasing order
c---------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine chibdy_xibdy(
     &                  chibdy,xibdy,
     &                  xextrap,yextrap,zextrap,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi)

        implicit none
        integer numbdypoints
        real*8 xextrap(numbdypoints)
        real*8 yextrap(numbdypoints)
        real*8 zextrap(numbdypoints)

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 chibdy(bdy_Nchi)
        real*8 xibdy(bdy_Nxi)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer i,j

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        !-----------------------------------------------------------------------

        chiextrap_min=minval(chiextrap)-1
        chiextrap_max=maxval(chiextrap)
        xiextrap_min=minval(xiextrap)-1
        xiextrap_max=maxval(xiextrap)

        do i=1,bdy_Nchi
            chiextrap_min=
     &                minval(chiextrap,mask=chiextrap.gt.chiextrap_min)
            chibdy(i)=chiextrap_min
        end do

        do j=1,bdy_Nxi 
            xiextrap_min=minval(xiextrap,mask=xiextrap.gt.xiextrap_min)
            xibdy(j)=xiextrap_min
        end do

        return
        end
!----------------------------------------------------------------------
