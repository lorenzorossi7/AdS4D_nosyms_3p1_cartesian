c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different from ex at grid points that are at large rho=sqrt(x**2+y**2+z**2) and next to excised points, so we're identifying the last not excised grid points near the boundary. We will use them to extrapolate the value of the quasi-local boundary stress-energy tensor at the boundary
c----------------------------------------------------------------------

        subroutine nexttobdypoints(
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

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1
        real*8 maxxyzp1
        integer numbdypoints

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

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

        numbdypoints=0

        do i=is,ie
         do j=js,je
           do k=ks,ke

           ! sets xp1,yp1,zp1 to values at i,j,k
           x0=x(i)
           y0=y(j)
           z0=z(k)
           rho0=sqrt(x0**2+y0**2+z0**2)

!           !chrbdy(i,j,k) is not ex for points near the boundary
           if ((chr(i,j,k).ne.ex).and.(rho0.ge.(1.0d0-3*dx/2))) then
              chrbdy(i,j,k)=ex-1
           else
              chrbdy(i,j,k)=ex
           end if

           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

!           !chrbdy(i,j,k) is not ex only for points near the boundary AND next to excised points   

           if (chrbdy(i,j,k).ne.ex) then
            if (maxxyzp1.eq.abs(xp1)) then
             if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
               chrbdy(i,j,k)=ex
             end if
            else if (maxxyzp1.eq.abs(yp1)) then
             if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
               chrbdy(i,j,k)=ex
             end if
            else
             if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
               chrbdy(i,j,k)=ex
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

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)

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
            end if
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------
