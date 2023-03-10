c----------------------------------------------------------------------
c miscellaneous numerical routines
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c specific x/y first derivative routines, including support
c for excision. (separate AdS boundaries no longer supported)
c
c only "method 4" copied over here from gh3d ... if that's 
c good enough simplifies routines a bit. Also,
c can easily add back y derivatives.
c-----------------------------------------------------------------------
        subroutine df1_int_x(f,f_x,x,y,i,j,chr,ex,Nx,Ny)
        implicit none
        integer Nx,Ny,i,j
        real*8 f(Nx,Ny),chr(Nx,Ny),ex,f_x,x(Nx),y(Ny)

        real*8 dx

        !--------------------------------------------------------------

        logical first
        save first
        data first/.true./

        dx=x(2)-x(1)

        if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
           if (i.le.(Nx-3)
     &         .and.((chr(i+1,j).ne.ex
     &         .and.chr(i+2,j).ne.ex
     &         .and.chr(i+3,j).ne.ex))) then
             f_x=(-4*f(i,j)+7*f(i+1,j)-4*f(i+2,j)+f(i+3,j))/2/dx
           else if (i.le.(Nx-2)
     &              .and.((chr(i+1,j).ne.ex
     &              .and.chr(i+2,j).ne.ex))) then
              f_x=(-3*f(i,j)+4*f(i+1,j)-f(i+2,j))/2/dx
           else if (i.le.(Nx-1).and.chr(i+1,j).ne.ex) then
              f_x=(-f(i,j)+f(i+1,j))/dx
!              write(*,*) 'df1_int_x: warning ... i=1 first order'
!              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (A)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
           if (i.ge.4
     &         .and.((chr(i-1,j).ne.ex
     &         .and.chr(i-2,j).ne.ex
     &         .and.chr(i-3,j).ne.ex))) then
             f_x=(4*f(i,j)-7*f(i-1,j)+4*f(i-2,j)-f(i-3,j))/2/dx
           else if (i.ge.3
     &              .and.((chr(i-1,j).ne.ex
     &              .and.chr(i-2,j).ne.ex))) then
              f_x=(3*f(i,j)-4*f(i-1,j)+f(i-2,j))/2/dx
           else if (i.ge.2.and.chr(i-1,j).ne.ex) then
              f_x=(f(i,j)-f(i-1,j))/dx
!              write(*,*) 'df1_int: warning ... i=Nx first order'
!              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (B)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else
           if ((chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex)) then
              f_x=(f(i+1,j)-f(i-1,j))/2/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (C)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        end if

        return
        end
!----------------------------------------------------------------------
        subroutine df1_int_y(f,f_y,x,y,i,j,chr,ex,Nx,Ny)
        implicit none
        integer Nx,Ny,i,j
        real*8 f(Nx,Ny),chr(Nx,Ny),ex,f_y,x(Nx),y(Ny)

        real*8 dy
        logical first
        save first
        data first/.true./

        !--------------------------------------------------------------

        dy=y(2)-y(1)

        if ((j.eq.1).or.(chr(i,j-1).eq.ex)) then
           if (j.le.(Ny-3)
     &         .and.((chr(i,j+1).ne.ex
     &         .and.chr(i,j+2).ne.ex
     &         .and.chr(i,j+3).ne.ex))) then
             f_y=(-4*f(i,j)+7*f(i,j+1)-4*f(i,j+2)+f(i,j+3))/2/dy
           else if (j.le.(Ny-2).and.((chr(i,j+1).ne.ex
     &              .and.chr(i,j+2).ne.ex))) then
              f_y=(-3*f(i,j)+4*f(i,j+1)-f(i,j+2))/2/dy
           else if (j.le.(Ny-1).and.chr(i,j+1).ne.ex) then
              f_y=(-f(i,j)+f(i,j+1))/dy
!              write(*,*) 'df1_int_y: warning ... j=1 first order'
!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (D)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        else if ((j.eq.Ny).or.(chr(i,j+1).eq.ex)) then
           if (j.ge.4
     &         .and.((chr(i,j-1).ne.ex
     &         .and.chr(i,j-2).ne.ex
     &         .and.chr(i,j-3).ne.ex))) then
             f_y=(4*f(i,j)-7*f(i,j-1)+4*f(i,j-2)-f(i,j-3))/2/dy
           else if (j.ge.3.and.((chr(i,j-1).ne.ex
     &              .and.chr(i,j-2).ne.ex))) then
              f_y=(3*f(i,j)-4*f(i,j-1)+f(i,j-2))/2/dy
           else if (j.ge.2.and.chr(i,j-1).ne.ex) then
              f_y=(f(i,j)-f(i,j-1))/dy
!              write(*,*) 'df1_int_y: warning ... j=Ny first order'
!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (E)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        else
           if ((chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex)) then
              f_y=(f(i,j+1)-f(i,j-1))/2/dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (F)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        end if

        return
        end

c----------------------------------------------------------------------
c the following computes all first derivatives of f,
c at a point i,j, at time level n.
c----------------------------------------------------------------------
        subroutine df1_int(f_np1,f_n,f_nm1,f_t,f_x,f_y,
     &                     x,y,dt,i,j,chr,ex,Nx,Ny,name)
        implicit none
        integer Nx,Ny,i,j
        real*8 f_np1(Nx,Ny),f_n(Nx,Ny),f_nm1(Nx,Ny)
        real*8 f_t,f_x,f_y,x(Nx),y(Ny),dt,ex,chr(Nx,Ny)
        character*(*) name

        logical first
        save first
        data first/.true./

        logical ltrace
        parameter (ltrace=.false.)

        !--------------------------------------------------------------

        if (chr(i,j).eq.ex) then
          write(*,*) 'df1_int: error ... point excised'
          stop
        end if

        f_t=(f_np1(i,j)-f_nm1(i,j))/2/dt
  
        call df1_int_x(f_n,f_x,x,y,i,j,chr,ex,Nx,Ny)
        call df1_int_y(f_n,f_y,x,y,i,j,chr,ex,Nx,Ny)

        if (ltrace) then
           write(*,*) 'df1_int for ',name
           write(*,*) ' f_t=',f_t
           write(*,*) ' f_x=',f_x
           write(*,*) ' f_y=',f_y
        end if

        return
        end

c----------------------------------------------------------------------
c same as df1_int above, but computes all second derivatives as well 
c
c if (extrap), then we use TE matched first derivative
c operators, plus 2nd order accurate 2nd derivatives ... this
c (with *NO* boundary dissipation) is equivalent to using 
c 4th order extrapolation (along the given direction where the 
c derivative is evaluated) 2-points into the excision surface, and then
c using the standard update scheme (now with interior dissipation)
c applied to the 'old' boundaries.
c
c For AdS: forget now whether extrap was essential for gh3d or not.
c          leaving on for now
c
c CALLING grid function f_n f here, to avoid conflict
c with f_n in the include stuff
c----------------------------------------------------------------------
        subroutine df2_int(f_np1,f,f_nm1,f_t,f_x,f_y,
     &                     f_tt,f_tx,f_ty,f_xx,f_xy,
     &                     f_yy,x,y,dt,i,j,chr,ex,Nx,Ny,name)
        implicit none
        integer Nx,Ny,i,j
        real*8 f_np1(Nx,Ny),f(Nx,Ny),f_nm1(Nx,Ny)
        real*8 f_t,f_x,f_y,f_tt,f_tx,f_ty,f_xx,f_xy
        real*8 f_yy,x(Nx),y(Ny),dt,ex,chr(Nx,Ny)
        character*(*) name

        real*8 dx,dy

        logical first,extrap
        save first
        data first/.true./
        parameter (extrap=.true.)

        logical ltrace
        parameter (ltrace=.false.)

        real*8 f_x_np1,f_x_nm1,f_y_np1,f_y_nm1
        real*8 f_y_ip1,f_y_ip2
        real*8 f_y_im1,f_y_im2

        ! initialize fixed-size variables
        data f_x_np1,f_x_nm1,f_y_np1,f_y_nm1/0.0,0.0,0.0,0.0/
        data f_y_ip1,f_y_ip2/0.0,0.0/
        data f_y_im1,f_y_im2/0.0,0.0/

        !--------------------------------------------------------------
       
        call df1_int(f_np1,f,f_nm1,f_t,f_x,f_y,
     &               x,y,dt,i,j,chr,ex,Nx,Ny,name)

        f_tt=(f_np1(i,j)-2*f(i,j)+f_nm1(i,j))/dt/dt 

        f_xx=0
        f_xy=0
        f_yy=0

        if (chr(i,j).eq.ex) then
         write(*,*) 'df2_int: error ... point excised'
         stop
        end if

        call df1_int_x(f_np1,f_x_np1,x,y,i,j,chr,ex,Nx,Ny)
        call df1_int_x(f_nm1,f_x_nm1,x,y,i,j,chr,ex,Nx,Ny)

        f_tx=(f_x_np1-f_x_nm1)/2/dt

        call df1_int_y(f_np1,f_y_np1,x,y,i,j,chr,ex,Nx,Ny)
        call df1_int_y(f_nm1,f_y_nm1,x,y,i,j,chr,ex,Nx,Ny)

        f_ty=(f_y_np1-f_y_nm1)/2/dt

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        !i

        if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
           if (i.ge.(Nx-1).or.chr(i+1,j).eq.ex.or.chr(i+2,j).eq.ex) then
              if (first) then
                 first=.false.
                 write(*,*) 'df2_int: error in chr (A)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              return
           end if

           if (i.eq.(Nx-2).or.chr(i+3,j).eq.ex) then
              ! write(*,*) 'df2_int: warning ... first order i=1'
              ! write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
              f_xx=(f(i+2,j)-2*f(i+1,j)+f(i,j))/dx/dx 
           else if ((.not.extrap).and.i.lt.(Nx-3)
     &               .and.chr(i+4,j).ne.ex) then
              f_xx=(3*f(i,j)-9*f(i+1,j)+
     &             10*f(i+2,j)-5*f(i+3,j)+
     &             f(i+4,j))/dx/dx
           else
              f_xx=(2*f(i,j)-5*f(i+1,j)+
     &              4*f(i+2,j)-f(i+3,j))/dx/dx
           end if

           call df1_int_y(f,f_y_ip1,x,y,i+1,j,chr,ex,Nx,Ny)
           call df1_int_y(f,f_y_ip2,x,y,i+2,j,chr,ex,Nx,Ny)
           f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
        else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
           if (i.le.2.or.
     &        chr(i-1,j).eq.ex.or.chr(i-2,j).eq.ex) then
              if (first) then
                 first=.false.
                 write(*,*) 'df2_int: error in chr (B)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              return
           end if

           if (i.eq.3.or.chr(i-3,j).eq.ex) then
              ! write(*,*) 'df2_int: warning ... first order i=Nx'
              ! write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
              f_xx=(f(i,j)-2*f(i-1,j)+f(i-2,j))/dx/dx 
           else if ((.not.extrap).and.i.gt.4.and.chr(i-4,j).ne.ex) then
              f_xx=(3*f(i,j)-9*f(i-1,j)+
     &              10*f(i-2,j)-5*f(i-3,j)+
     &              f(i-4,j))/dx/dx
           else
              f_xx=(2*f(i,j)-5*f(i-1,j)+
     &              4*f(i-2,j)-f(i-3,j))/dx/dx
           end if

           call df1_int_y(f,f_y_im1,x,y,i-1,j,chr,ex,Nx,Ny)
           call df1_int_y(f,f_y_im2,x,y,i-2,j,chr,ex,Nx,Ny)
           f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
        else if (chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex) then

           f_xx=(f(i+1,j)-2*f(i,j)+f(i-1,j))/dx/dx 

           call df1_int_y(f,f_y_im1,x,y,i-1,j,chr,ex,Nx,Ny)
           call df1_int_y(f,f_y_ip1,x,y,i+1,j,chr,ex,Nx,Ny)
           f_xy=(f_y_ip1-f_y_im1)/2/dx
        else
           if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (C)'
              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
              write(*,*) '    (first error only)'
           end if
           return
        end if

        !j

        if (j.eq.1.or.(chr(i,j-1).eq.ex)) then
           if (j.ge.(Ny-1).or.chr(i,j+1).eq.ex.or.chr(i,j+2).eq.ex) then
              if (first) then
                 first=.false.
                 write(*,*) 'df2_int: error in chr (D)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              return
           end if

           if (j.eq.(Ny-2).or.chr(i,j+3).eq.ex) then
              !write(*,*) 'df2_int: warning ... first order j=1'
              !write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
              f_yy=(f(i,j+2)-2*f(i,j+1)+f(i,j))/dy/dy 
           else if ((.not.extrap).and.j.lt.(Ny-3).and.
     &               chr(i,j+4).ne.ex) then
              f_yy=(3*f(i,j)-9*f(i,j+1)+
     &              10*f(i,j+2)-5*f(i,j+3)+
     &              f(i,j+4))/dy/dy
           else
              f_yy=(2*f(i,j)-5*f(i,j+1)+
     &              4*f(i,j+2)-f(i,j+3))/dy/dy
           end if

        else if (j.eq.Ny.or.(chr(i,j+1).eq.ex)) then
           if (j.le.2.or.chr(i,j-1).eq.ex.or.chr(i,j-2).eq.ex) then
              if (first) then
                 first=.false.
                 write(*,*) 'df2_int: error in chr (E)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    name=',name
                 write(*,*) '    (first error only)'
              end if
              return
           end if

           if (j.eq.3.or.chr(i,j-3).eq.ex) then
              !write(*,*) 'df2_int: warning ... first order j=Ny'
              !write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
              f_yy=(f(i,j)-2*f(i,j-1)+f(i,j-2))/dy/dy 
           else if ((.not.extrap).and.j.gt.4.and.chr(i,j-4).ne.ex) then
              f_yy=(3*f(i,j)-9*f(i,j-1)+
     &              10*f(i,j-2)-5*f(i,j-3)+
     &              f(i,j-4))/dy/dy
           else
              f_yy=(2*f(i,j)-5*f(i,j-1)+
     &              4*f(i,j-2)-f(i,j-3))/dy/dy
           end if

        else if (chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex) then
            f_yy=(f(i,j+1)-2*f(i,j)+f(i,j-1))/dy/dy 

        else
            if (first) then
               first=.false.
               write(*,*) 'df2_int: error in chr (F)'
               write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
               write(*,*) '    (first error only)'
            end if
            return
        end if

        if (ltrace) then
           write(*,*) 'df2_int for ',name
           write(*,*) ' f_tt=',f_tt
           write(*,*) ' f_tx=',f_tx
           write(*,*) ' f_ty=',f_ty
           write(*,*) ' f_xx=',f_xx
           write(*,*) ' f_xy=',f_xy
           write(*,*) ' f_yy=',f_yy
        end if

        return
        end

c----------------------------------------------------------------------
c Background AdS values ... with these new variables
c the AdS metric has been factored into the maple already
c----------------------------------------------------------------------
        subroutine init_ghb_ads(gb_tt,gb_tx,gb_ty,gb_xx,
     &                      gb_xy,gb_yy,psi,Hb_t,Hb_x,Hb_y,
     &                      L,x,y,chr,ex,Nx,Ny,regtype)
        implicit none
        integer Nx,Ny
        integer regtype
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny),gb_yy(Nx,Ny)
        real*8 psi(Nx,Ny),tfunction(Nx,Ny),Hb_t(Nx,Ny)
        real*8 Hb_x(Nx,Ny),Hb_y(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L,x(Nx),y(Ny)

        integer i,j

        logical ltrace
        parameter (ltrace=.false.)
 
        ! initialize fixed-size variables
        data i,j/0,0/
  
        !--------------------------------------------------------------
 
        do i=1,Nx
           do j=1,Ny
              gb_tt(i,j)=0
              gb_tx(i,j)=0
              gb_ty(i,j)=0
              gb_xx(i,j)=0
              gb_xy(i,j)=0
              gb_yy(i,j)=0

              Hb_t(i,j)=0
              Hb_x(i,j)=0
              Hb_y(i,j)=0

              psi(i,j)=0
           end do
        end do

        call axi_reg_g(gb_tt,gb_tx,gb_ty,
     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
     &                 L,x,y,Nx,Ny,regtype)

        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,Nx,Ny,regtype)

        return
        end

c----------------------------------------------------------------------
c initializes f with a 2d gaussian-like profile ... multiplication
c by (1-rho^2)^3 is for correct asymptotics for AdS 
c
c f = A * (1-rho^2) exp (- (r-r0)^2/delta^2) , r>r0
c   = A , r < r0
c
c where r = sqrt ( (1-ex^2)*(x)^2 + (1-ey^2)*(y)^2 )
c----------------------------------------------------------------------
        subroutine gauss2d(f,A,B,r0,delta,xu0,yu0,ex,ey,L,x,y,Nx,Ny,
     &                     rhoc,rhod,stype)
        implicit none
        integer Nx,Ny
        real*8 f(Nx,Ny),x(Nx),y(Ny)
        real*8 A,B,r0,delta,ex,ey,xu0,yu0,L

        integer i,j
        integer stype
        real*8 r,x0,y0,rho0,chi0,csr,xb,yb

        real*8 rhoc,rhod
        real*8 f1,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j/0,0/
        data r,x0,y0,rho0,csr,xb,yb/0.0,0.0,0.0,0.0,0.0,0.0,0.0/
  
        !--------------------------------------------------------------
 

        do i=1,Nx
           do j=1,Ny
              f(i,j)=0
              x0=x(i)
              y0=y(j)
              rho0=sqrt(x0**2+y0**2)
              chi0=atan2(y0,x0)
              if (chi0.lt.0) chi0=chi0+2*PI

              xb=x0-xu0
              yb=y0-yu0
              r=sqrt(xb**2+yb**2-ex**2*xb**2-ey**2*yb**2)

              f1=trans(rho0,rhoc,rhod)

              ! Gaussian phi=amp*exp(-r^2)/r^4 profile 
              ! remember that phi=phi1*(1-rho^2)^3
              ! NOTE: compactified coordinate here is: rho

              if (stype.eq.0) then

                 ! note that (1-x^2)^4=(1+2r)^4/(1+r)^8
                 ! = (1-rho^2)^4 (1+rho(4-rho))^4 / (1+rho(2-rho))^8

                 ! note that x=r/(1+r) and r=2rho/(1-rho^2)
                 ! so x=r/(1+r)=2rho/(1+rho*(2-rho))

                 if (rho0.ge.1) then
                    f(i,j)=0
                 else if (2*r/(1+r*(2-r)).gt.r0) then
                    f(i,j)=
     &              A*exp(-((2*r/(1+r*(2-r))-r0)/delta)**2)
     &            *(1-rho0**2)*(1+rho0*(4-rho0))**4/(1+rho0*(2-rho0))**8
     &             +B*cos(chi0)*4*f1*(1-f1)
                 else
                    f(i,j)=
     &              A
     &            *(1-rho0**2)*(1+rho0*(4-rho0))**4/(1+rho0*(2-rho0))**8
     &             +B*cos(chi0)*4*f1*(1-f1)
                 end if

              else if (stype.eq.1) then

                 if (rho0.ge.1) then
                    f(i,j)=0
                 else if (rho0.gt.r0) then
                    f(i,j)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 else
                    f(i,j)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 end if

              else if (stype.eq.2) then !annulus ID

                 if (rho0.ge.1) then
                    f(i,j)=0
                 else if (rho0.gt.r0) then
                    f(i,j)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 else
                    f(i,j)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 end if

              end if

           end do
        end do

        return
        end

c----------------------------------------------------------------------
c initializes complex s.f. phi with approximate q-ball profile
c
c phi = amp * exp (-rho^2/delta^2 + i*omega*t) * (1-rho^2)
c
c phi = d(phi)/dt = i*omega*phi
c
c where rgo = sqrt ((x-x0)^2 + (y-y0)^2)
c
c and profile is multiplied by (1-rho^2)^3 for correct asymptotics
c----------------------------------------------------------------------
        subroutine approx_qb(phi_r,phi_i,phi_r_dot,phi_i_dot,amp,xu0,
     &                       delta,omega,v0,L,x,y,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 phi_r(Nx,Ny),phi_i(Nx,Ny)
        real*8 phi_r_dot(Nx,Ny),phi_i_dot(Nx,Ny)
        real*8 x(Nx),y(Ny),L,v0(2)
        real*8 amp,omega,delta,xu0(2)

        integer i,j
        real*8 r,x0,y0,rho0,xb,yb

        ! initialize fixed-size variables
        data i,j/0,0/
        data r,x0,y0,rho0,xb,yb/0.0,0.0,0.0,0.0,0.0,0.0/
  
        !--------------------------------------------------------------
 
        do i=1,Nx
           do j=1,Ny
              x0=x(i)
              y0=y(j)
              rho0=sqrt(x0**2+y0**2)

              xb=x0-xu0(1)
              yb=y0-xu0(2)
              r=sqrt(xb**2+yb**2)

              if (rho0.ge.1) then
                 phi_r(i,j)=0
                 phi_i(i,j)=0
                 phi_r_dot(i,j)=0
                 phi_i_dot(i,j)=0
              else
                 phi_r(i,j)=amp*exp(-(r/delta)**2)*(1-rho0**2)**2
                 phi_i(i,j)=0
                 phi_r_dot(i,j)=0
                 phi_i_dot(i,j)=omega*phi_r(i,j)
              end if

          end do
        end do

        return
        end

c-----------------------------------------------------------------------
c for variables with lin_zero_bnd ... zeros residual there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx,Ny)
        implicit none
        integer Nx,Ny,all
        real*8 f(Nx,Ny)
        integer phys_bdy(4)

        integer i,j,k,is,ie,js,je

        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je/0,0,0,0,0,0,0/
  
        !--------------------------------------------------------------
 
        if (phys_bdy(1).eq.1.or.all.eq.1) then
           do j=2,Ny-1
              f(2,j)=0
           end do
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
           do j=2,Ny-1
              f(Nx-1,j)=0
           end do
        end if

        if (phys_bdy(3).eq.1.or.all.eq.1) then
           do i=3,Nx-2
              f(i,2)=0
           end do
        end if

        if (phys_bdy(4).eq.1.or.all.eq.1) then
           do i=3,Nx-2
              f(i,Ny-1)=0
           end do
        end if

        return
        end

c---------------------------------------------------------------------------
c applies a simple KO filter to f and modified
c near excision boundaries as follows:
c (or all boundaries if do_bdy=1, 
c  or neither if do_bdy=0 and internal flag no_uc_ex_bdy=.true.)
c
c with undivided operators Up(f) = f(i+1)-f(i)
c                          Um(f) = f(i)-f(i-1)
c
c Interior: eps*w4*(Up Um)^2
c left+1  : eps*w3*(Up Up Um)
c left    : eps*w2*(Up Um)
c right-1 : eps*w3*(Um Um Up)
c right   : eps*w2*(Up Um)
c
c NOTE: SPECIAL FLAG ... if do_ex>0, eps is treated as a scalar, 
c       else eps is treated as an array.
c
c same as dmdiss3d_ex(), but with 'kmax' flag, applying separate
c sweeps from kmax.eq.1
c------------------------------------------------------------------------
        subroutine dmdiss3d_ex_gen(f,work,eps,do_bdy,phys_bdy_type,even,
     &             odd,nx,ny,nz,chr,ex,do_ex,ind_sweeps,kmax)
        implicit none
        integer nx,ny,nz,do_bdy,phys_bdy_type(6),even,odd,do_ex,kmax
        real*8 f(nx,ny,nz),work(nx,ny,nz),eps(nx,ny,nz)
        real*8 chr(nx,ny,nz),ex
        logical ind_sweeps

        integer i,j,k,bo1,bo2
        integer pass,npass,kd
        real*8 eps_eff,f_hf,norm_f
        real*8 w4,w3,w2
        parameter (w4=1.0d0/16.0d0,w3=1.0d0/8.0d0,w2=1.0d0/4.0d0)
        logical no_uc_ex_diss
        parameter (no_uc_ex_diss=.true.)

        ! initialize fixed-size variables
        data i,j,k,bo1,bo2/0,0,0,0,0/
        data pass,npass,kd/0,0,0/
        data eps_eff,f_hf,norm_f/0.0,0.0,0.0/

        !--------------------------------------------------------------
        
        eps_eff=eps(1,1,1)

        if (do_bdy.eq.0) then
           bo1=0
           bo2=0
           if (no_uc_ex_diss) then
              bo1=-max(nx,ny,nz)
              bo2=bo1
           end if
        else
           bo1=1
           bo2=2
        end if

        do kd=kmax,1,-1

        npass=1
        if (ny.gt.1.and.ind_sweeps) npass=2
        if (nz.gt.1.and.ind_sweeps) npass=3

        do pass=1,npass

         do i=1,nx
           do j=1,ny
              do k=1,nz
                 work(i,j,k)=f(i,j,k)
              end do
           end do
         end do

         do i=1,nx
          do j=1,ny
           do k=1,nz
             if ((chr(i,j,k).ne.ex)) then

               if (do_ex.lt.0) eps_eff=eps(i,j,k)

               if (.not.ind_sweeps.or.pass.eq.1) then
                f_hf=0
                if (i.ge.(2*kd+1).and.i.le.(nx-2*kd).and.
     &           ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex
     &             .and.chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)))
     &          then
                   f_hf=w4*(
     &                 work(i-2*kd,j,k)+work(i+2*kd,j,k)
     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.ge.(kd+1).and.i.lt.(2*kd+1).and.
     &                   phys_bdy_type(1).eq.odd.and.
     &                ((chr(i-kd,j,k).ne.ex.and.
     &                 chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                (-work(i,j,k))+work(i+2*kd,j,k)
     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.lt.(kd+1).and.
     &                   phys_bdy_type(1).eq.odd.and.
     &               ((chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                (-work(i+2*kd,j,k))+work(i+2*kd,j,k)
     &          -4*((-work(i+1*kd,j,k))+work(i+1*kd,j,k))+6*work(i,j,k))
                else if (i.gt.(Nx-2*kd).and.i.le.(Nx-kd)
     &                   .and.phys_bdy_type(2).eq.odd.and.
     &              ((chr(i+kd,j,k).ne.ex.and.
     &                chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                 work(i-2*kd,j,k)+(-work(i,j,k))
     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.gt.(Nx-kd).and.phys_bdy_type(2).eq.odd.and.
     &             ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                  work(i-2*kd,j,k)+(-work(i-2*kd,j,k))
     &              -4*(work(i-kd,j,k)+(-work(i-kd,j,k)))+6*work(i,j,k))
                else if (i.ge.(kd+1).and.i.lt.(2*kd+1)
     &                   .and.phys_bdy_type(1).eq.even.and.
     &              ((chr(i-kd,j,k).ne.ex.and.
     &                 chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                    (work(i,j,k))+work(i+2*kd,j,k)
     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.lt.(kd+1).and.phys_bdy_type(1).eq.even.and.
     &               ((chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                   (work(i+2*kd,j,k))+work(i+2*kd,j,k)
     &               -4*((work(i+kd,j,k))+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.gt.(Nx-2*kd).and.i.le.(Nx-kd)
     *                   .and.phys_bdy_type(2).eq.even.and.
     &              ((chr(i+kd,j,k).ne.ex.and.
     &                 chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                     work(i-2*kd,j,k)+(work(i,j,k))
     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
                else if (i.gt.(Nx-kd).and.phys_bdy_type(2).eq.even.and.
     &             ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                   work(i-2*kd,j,k)+(work(i-2*kd,j,k))
     &               -4*(work(i-kd,j,k)+(work(i-kd,j,k)))+6*work(i,j,k))
                else if (i.gt.(2-bo1)*kd.and.i.le.(nx-2*kd).and.
     &                 (chr(i-kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex.and.
     &                  chr(i+2*kd,j,k).ne.ex)) 
     &          then
                   f_hf=w3*(
     &                    -work(i-kd,j,k)+3*work(i,j,k)
     &                  -3*work(i+kd,j,k)+work(i+2*kd,j,k))
                else if (i.gt.2*kd.and.i.le.(nx-2*kd+bo1*kd).and.
     &                 (chr(i-kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex.and.
     &                   chr(i-2*kd,j,k).ne.ex) ) 
     &          then
                   f_hf=w3*(
     &                    -work(i+kd,j,k)+3*work(i,j,k)
     &                  -3*work(i-kd,j,k)+work(i-2*kd,j,k))
                else if ((i.gt.(2-bo2)*kd.or.   
     &           (i.le.2*kd.and.chr(max(1,i-kd),j,k).eq.ex))
     &                       .and.i.le.(nx-2*kd).and.
     &                (chr(i+kd,j,k).ne.ex.and.chr(i+2*kd,j,k).ne.ex) )
     &          then
                   f_hf=w2*(
     &                   work(i,j,k)-2*work(i+kd,j,k)+work(i+2*kd,j,k))
                else if (i.gt.2*kd.and.(i.lt.(nx-kd+bo2*kd).or.
     &              (i.ge.(nx-kd).and.chr(min(nx,i+kd),j,k).eq.ex)).and.
     &                 (chr(i-kd,j,k).ne.ex.and.chr(i-2*kd,j,k).ne.ex) )
     &          then
                   f_hf=w2*(
     &                    work(i,j,k)-2*work(i-kd,j,k)+work(i-2*kd,j,k))
                end if

                f(i,j,k)=f(i,j,k)-eps_eff*f_hf
               end if

               if (.not.ind_sweeps.or.pass.eq.2) then
                f_hf=0
                if (j.ge.(2*kd+1).and.j.le.(ny-2*kd).and.
     &           ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex
     &             .and.chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)))
     &          then
                   f_hf=w4*(
     &                 work(i,j-2*kd,k)+work(i,j+2*kd,k)
     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.ge.(kd+1).and.j.lt.(2*kd+1).and.
     &                   phys_bdy_type(3).eq.odd.and.
     &                ((chr(i,j-kd,k).ne.ex.and.
     &                 chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                (-work(i,j,k))+work(i,j+2*kd,k)
     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.lt.(kd+1).and.
     &                   phys_bdy_type(3).eq.odd.and.
     &               ((chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                (-work(i,j+2*kd,k))+work(i,j+2*kd,k)
     &          -4*((-work(i,j+1*kd,k))+work(i,j+1*kd,k))+6*work(i,j,k))
                else if (j.gt.(ny-2*kd).and.j.le.(ny-kd)
     &                   .and.phys_bdy_type(4).eq.odd.and.
     &              ((chr(i,j+kd,k).ne.ex.and.
     &                chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                 work(i,j-2*kd,k)+(-work(i,j,k))
     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.gt.(Ny-kd).and.phys_bdy_type(4).eq.odd.and.
     &             ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                  work(i,j-2*kd,k)+(-work(i,j-2*kd,k))
     &              -4*(work(i,j-kd,k)+(-work(i,j-kd,k)))+6*work(i,j,k))
                else if (j.ge.(kd+1).and.j.lt.(2*kd+1)
     &                   .and.phys_bdy_type(3).eq.even.and.
     &              ((chr(i,j-kd,k).ne.ex.and.
     &                 chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                    (work(i,j,k))+work(i,j+2*kd,k)
     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.lt.(kd+1).and.phys_bdy_type(3).eq.even.and.
     &               ((chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                   (work(i,j+2*kd,k))+work(i,j+2*kd,k)
     &               -4*((work(i,j+kd,k))+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.gt.(Ny-2*kd).and.j.le.(Ny-kd)
     *                   .and.phys_bdy_type(4).eq.even.and.
     &              ((chr(i,j+kd,k).ne.ex.and.
     &                 chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                     work(i,j-2*kd,k)+(work(i,j,k))
     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
                else if (j.gt.(Ny-kd).and.phys_bdy_type(4).eq.even.and.
     &             ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
     &          then
                   f_hf=w4*(
     &                   work(i,j-2*kd,k)+(work(i,j-2*kd,k))
     &               -4*(work(i,j-kd,k)+(work(i,j-kd,k)))+6*work(i,j,k))
                else if (j.gt.(2-bo1)*kd.and.j.le.(ny-2*kd).and.
     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex.and.
     &                  chr(i,j+2*kd,k).ne.ex)) 
     &          then
                   f_hf=w3*(
     &                    -work(i,j-kd,k)+3*work(i,j,k)
     &                  -3*work(i,j+kd,k)+work(i,j+2*kd,k))
                else if (j.gt.2*kd.and.j.le.(ny-2*kd+bo1*kd).and.
     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex.and.
     &                   chr(i,j-2*kd,k).ne.ex) ) 
     &          then
                   f_hf=w3*(
     &                    -work(i,j+kd,k)+3*work(i,j,k)
     &                  -3*work(i,j-kd,k)+work(i,j-2*kd,k))
                else if ((j.gt.(2-bo2)*kd.or.   
     &           (j.le.2*kd.and.chr(i,max(1,j-kd),k).eq.ex))
     &                       .and.j.le.(ny-2*kd).and.
     &                (chr(i,j+kd,k).ne.ex.and.chr(i,j+2*kd,k).ne.ex) )
     &          then
                   f_hf=w2*(
     &                   work(i,j,k)-2*work(i,j+kd,k)+work(i,j+2*kd,k))
                else if (j.gt.2*kd.and.(j.lt.(ny-kd+bo2*kd).or.
     &              (j.ge.(ny-kd).and.chr(i,min(ny,j+kd),k).eq.ex)).and.
     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j-2*kd,k).ne.ex) )
     &          then
                   f_hf=w2*(
     &                    work(i,j,k)-2*work(i,j-kd,k)+work(i,j-2*kd,k))
                end if

                f(i,j,k)=f(i,j,k)-eps_eff*f_hf
               end if

               if ((.not.ind_sweeps.or.pass.eq.3).and.Nz.gt.1) then
                f_hf=0
                write(*,*) 'dmdiss3d_ex_gen: not yet updated for 3D'
                stop
               end if

             end if
            end do
           end do
         end do
        end do

        end do !kd

        return
        end

c----------------------------------------------------------------------
c initializes the metric to an exact black hole solution
c with radius parameter r0
c----------------------------------------------------------------------
        subroutine init_ads5d_bh(r0,L,gb_tt,gb_tx,gb_ty,gb_xx,
     &                         gb_xy,gb_yy,psi,
     &                         gb_tt_t,gb_tx_t,gb_ty_t,gb_xx_t,gb_xy_t,
     &                         gb_yy_t,psi_t,
     &                         Hb_t,Hb_x,Hb_y,
     &                         Hb_t_t,Hb_x_t,Hb_y_t,
     &                         phys_bdy,x,y,dt,chr,ex,Nx,Ny,regtype)
        implicit none

        integer Nx,Ny
        integer regtype
        integer phys_bdy(6)
        real*8 r0
        real*8 dt,ex,L
        real*8 chr(Nx,Ny)
        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny)
        real*8 Hb_y(Nx,Ny)
        real*8 Hb_t_t(Nx,Ny),Hb_x_t(Nx,Ny)
        real*8 Hb_y_t(Nx,Ny)
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny)
        real*8 gb_ty(Nx,Ny)
        real*8 gb_xx(Nx,Ny),gb_xy(Nx,Ny)
        real*8 gb_yy(Nx,Ny),psi(Nx,Ny)
        real*8 gb_tt_t(Nx,Ny),gb_tx_t(Nx,Ny)
        real*8 gb_ty_t(Nx,Ny),psi_t(Nx,Ny)
        real*8 gb_xx_t(Nx,Ny),gb_xy_t(Nx,Ny)
        real*8 gb_yy_t(Nx,Ny)
        real*8 x(Nx),y(Ny)

        integer n
        parameter (n=6)

        integer i,j

        real*8 rho0,f1,f0,cF0,C0,A0,B0,D0
        real*8 x0,y0
        real*8 r_h,rho_h
        real*8 small
        parameter (small=1d-10)

        real*8 tfunction(Nx,Ny)

        ! initialize fixed-size variables
        data i,j/0,0/

        data rho0,f1,f0,cF0/0.0,0.0,0.0,0.0/
        data C0,A0,B0,D0/0.0,0.0,0.0,0.0/
        data x0,y0/0.0,0.0/
        data r_h,rho_h/0.0,0.0/

        !--------------------------------------------------------------

        ! compute horizon global radius r_h and corresponding compactified rho_h
        !(NOTE ... here the x0,y0,rho0 are compactified (CODE)
        ! versions of the coordinates)
        r_h=(-2*3**0.3333333333333333*L**2 + 
     &      2**0.3333333333333333*(9*L**2*r0 + Sqrt(12*L**6 +
     &      81*L**4*r0**2))**
     &      0.6666666666666666)/
     &      (6**0.6666666666666666*(9*L**2*r0 + Sqrt(12*L**6 +
     &      81*L**4*r0**2))**
     &      0.3333333333333333)
        rho_h=(-1 + sqrt(1 + r_h**2))/r_h

        ! initialize metric 
        do i=1,Nx
           do j=1,Ny
              gb_tt_t(i,j)=0
              gb_tx_t(i,j)=0
              gb_ty_t(i,j)=0
              gb_xx_t(i,j)=0
              gb_xy_t(i,j)=0
              gb_yy_t(i,j)=0
              psi_t(i,j)=0
              Hb_t_t(i,j)=0
              Hb_x_t(i,j)=0
              Hb_y_t(i,j)=0
              if (chr(i,j).eq.ex) then
                 gb_tt(i,j)=0
                 gb_tx(i,j)=0
                 gb_ty(i,j)=0
                 gb_xx(i,j)=0
                 gb_xy(i,j)=0
                 gb_yy(i,j)=0
                 psi(i,j)=0
                 Hb_t(i,j)=0
                 Hb_x(i,j)=0
                 Hb_y(i,j)=0
              else
                 x0=x(i)
                 y0=y(j)
                 rho0=sqrt(x0**2+y0**2)

                 ! EF-like-near-horizon Schwarzschild-like-near-bdy coordinates
                 ! TODO: add AdS_L dependence; currently assumes AdS_L=1!
                 gb_tt(i,j)=-(r0*(-1 + rho0**2))/(2*rho0)
                 gb_tx(i,j)=(2*(-1 + rho0)**2*(1 + rho0**2)*x0)/
     &                      (rho0*(1 + rho0)**2*(-1 + rho_h)**4)
                 gb_ty(i,j)=(2*(-1 + rho0)**2*(1 + rho0**2)*y0)/
     &                      (rho0*(1 + rho0)**2*(-1 + rho_h)**4)
                 gb_xx(i,j)=(4*(-(-1 + rho0**2)**2 - (2*(-1 +
     &                      rho0**4)**2*(rho0 - rho_h)*(-2+rho0+rho_h)*
     &                      (2 + (-2 + rho0)*rho0 + (-2 + rho_h)*rho_h)*
     &                      (2 + (-2 + rho0)*rho0*(2 +(-2+rho0)*rho0) + 
     &                      (-2+rho_h)*rho_h*
     &                      (2+(-2+rho_h)*rho_h))*x0**2)/
     &                      (rho0*(r0*(-1 + rho0**2)**3 + 
     &                      2*rho0*(1 + rho0**2)**2)*(-1 + rho_h)**8) + 
     &                      ((-1 + rho0**2)**2*y0**2)/rho0**2))/
     &                      (-1 + rho0**2)**4
                 gb_xy(i,j)=(4*(-1 + (rho0*(2 - 2*r0*rho0 + 4*rho0**2 +
     &                      7*r0*rho0**3 + 2*rho0**4 - 
     &                      9*r0*rho0**5 + 5*r0*rho0**7 - r0*rho0**9 - 
     &                      (2*(-1 + rho0)**10*(1 + rho0 + rho0**2 +
     &                      rho0**3)**2)/(-1 + rho_h)**8))/
     &                      (r0*(-1 + rho0**2)**3 + 
     &                      2*rho0*(1 + rho0**2)**2))*x0*y0)/
     &                      (rho0**2*(-1 + rho0**2)**4)
                 gb_yy(i,j)=(4*(-rho0**2 + x0**2 - (2*rho0*(1 +
     &                      rho0**2)**2*(rho0 - rho_h)*(-2+rho0+rho_h)*
     &                      (2 + (-2 + rho0)*rho0 + (-2 + rho_h)*rho_h)*
     &                      (2 + (-2 + rho0)*rho0*(2 +(-2+rho0)*rho0) + 
     &                      (-2+rho_h)*rho_h*
     &                      (2+(-2+rho_h)*rho_h))*y0**2)/
     &                      ((r0*(-1+rho0**2)**3+2*rho0*(1+rho0**2)**2)*
     &                      (-1 +rho_h)**8)))/
     &                      (rho0 - rho0**3)**2

!                 ! (Schw coordinates)!
!                 ! TODO: add AdS_L dependence; currently assumes AdS_L=1!
!                 gb_tt(i,j)=0
!                 gb_xx(i,j)=0
!                 gb_xy(i,j)=0
!                 gb_yy(i,j)=0

              end if
           end do
        end do

        ! y=0 axis regularization
        call axi_reg_g(gb_tt,gb_tx,gb_ty,
     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
     &                 L,x,y,Nx,Ny,regtype)

        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,Nx,Ny,regtype)

        return
        end

c----------------------------------------------------------------------
c calculates inverse of a 4x4 matrix with components g11,...,g44
c----------------------------------------------------------------------
        subroutine calc_g0uu(g11,g12,g13,g14,g22,g23,g24,
     &                       g33,g34,g44,g0_uu)
        implicit none
        real*8 g11,g12,g13,g14
        real*8 g22,g23,g24
        real*8 g33,g34
        real*8 g44
        real*8 g0_uu(4,4)                        

        real*8 invdenominator

        !--------------------------------------------------------------

        invdenominator=
     &      g14**2*g23**2 - 2*g13*g14*g23*g24 + g13**2*g24**2 - 
     &      g14**2*g22*g33 + 2*g12*g14*g24*g33 - g11*g24**2*g33 + 
     &      2*g13*g14*g22*g34 - 2*g12*g14*g23*g34 - 2*g12*g13*g24*g34 + 
     &      2*g11*g23*g24*g34 + g12**2*g34**2 - g11*g22*g34**2 - 
     &      g13**2*g22*g44 + 2*g12*g13*g23*g44 - g11*g23**2*g44 - 
     &      g12**2*g33*g44 + g11*g22*g33*g44

        g0_uu(1,1)=
     &      (
     &      -(g24**2*g33) + 2*g23*g24*g34 - g22*g34**2 - g23**2*g44 
     &      + g22*g33*g44
     &      )
     &      /invdenominator
        g0_uu(1,2)=
     &      (
     &      g14*g24*g33 - g14*g23*g34 - g13*g24*g34 + g12*g34**2 + 
     &      g13*g23*g44 - g12*g33*g44
     &      )
     &      /invdenominator
        g0_uu(1,3)=
     &      (
     &      -(g14*g23*g24) + g13*g24**2 + g14*g22*g34 - g12*g24*g34 - 
     &      g13*g22*g44 + g12*g23*g44
     &      )
     &      /invdenominator
        g0_uu(1,4)=
     &      (
     &      g14*g23**2 - g13*g23*g24 - g14*g22*g33 + g12*g24*g33 + 
     &      g13*g22*g34 - g12*g23*g34
     &      )
     &      /invdenominator
        g0_uu(2,2)=
     &      (
     &      -(g14**2*g33) + 2*g13*g14*g34 - g11*g34**2 - g13**2*g44 
     &      + g11*g33*g44
     &      )
     &      /invdenominator
        g0_uu(2,3)=
     &      (
     &      g14**2*g23 - g13*g14*g24 - g12*g14*g34 + g11*g24*g34 + 
     &      g12*g13*g44 - g11*g23*g44
     &      )
     &      /invdenominator
        g0_uu(2,4)=
     &      (
     &      -(g13*g14*g23) + g13**2*g24 + g12*g14*g33 - g11*g24*g33 - 
     &      g12*g13*g34 + g11*g23*g34
     &      )
     &      /invdenominator
        g0_uu(3,3)=
     &      (
     &      -(g14**2*g22) + 2*g12*g14*g24 - g11*g24**2 - g12**2*g44 
     &      + g11*g22*g44
     &      )
     &      /invdenominator
        g0_uu(3,4)=
     &      (
     &      g13*g14*g22 - g12*g14*g23 - g12*g13*g24 + g11*g23*g24 + 
     &      g12**2*g34 - g11*g22*g34
     &      )
     &      /invdenominator
        g0_uu(4,4)=
     &      (
     &      -(g13**2*g22) + 2*g12*g13*g23 - g11*g23**2 - g12**2*g33 
     &      + g11*g22*g33
     &      )
     &      /invdenominator
        

        return
        end  

c----------------------------------------------------------------------
c calculates all the tensorial objects in x,y coordinates, at point i,j
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  A_l,A_l_x,Hads_l,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  einstein_ll,set_ll,
     &                  phi10_x,phi10_xx,
     &                  x,y,dt,chr,L,ex,Nx,Ny,i,j)
        implicit none

        integer Nx,Ny
        integer i,j,k

        real*8 chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),dt,L

        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)

        integer a,b,c,d,e,f,g,h
        real*8 dx,dy
        real*8 x0,y0
        real*8 rho0
        real*8 f0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 grad_phi1_sq

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t, gb_tt_x, gb_tt_y
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty
        real*8 gb_tt_xx,gb_tt_xy
        real*8 gb_tt_yy
        real*8 gb_tx_t, gb_tx_x, gb_tx_y
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty
        real*8 gb_tx_xx,gb_tx_xy
        real*8 gb_tx_yy
        real*8 gb_ty_t, gb_ty_x, gb_ty_y
        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty
        real*8 gb_ty_xx,gb_ty_xy
        real*8 gb_ty_yy
        real*8 gb_xx_t, gb_xx_x, gb_xx_y
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty
        real*8 gb_xx_xx,gb_xx_xy
        real*8 gb_xx_yy
        real*8 gb_xy_t, gb_xy_x, gb_xy_y
        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty
        real*8 gb_xy_xx,gb_xy_xy
        real*8 gb_xy_yy
        real*8 gb_yy_t, gb_yy_x, gb_yy_y
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty
        real*8 gb_yy_xx,gb_yy_xy
        real*8 gb_yy_yy
        real*8 psi_t, psi_x, psi_y
        real*8 psi_tt,psi_tx,psi_ty
        real*8 psi_xx,psi_xy
        real*8 psi_yy
        real*8 phi1_t, phi1_x, phi1_y
        real*8 phi1_tt,phi1_tx,phi1_ty
        real*8 phi1_xx,phi1_xy
        real*8 phi1_yy

        real*8 gb_tt0,gb_tx0,gb_ty0
        real*8 gb_xx0,gb_xy0
        real*8 gb_yy0
        real*8 psi0
        real*8 phi10

        real*8 g0_tt_ads_x,g0_tt_ads_xx,g0_tt_ads_xy
        real*8 g0_tt_ads_y,g0_tt_ads_yy
        real*8 g0_xx_ads_x,g0_xx_ads_xx,g0_xx_ads_xy
        real*8 g0_xx_ads_y,g0_xx_ads_yy
        real*8 g0_xy_ads_x,g0_xy_ads_xx,g0_xy_ads_xy
        real*8 g0_xy_ads_y,g0_xy_ads_yy
        real*8 g0_yy_ads_x,g0_yy_ads_xx,g0_yy_ads_xy
        real*8 g0_yy_ads_y,g0_yy_ads_yy
        real*8 g0_psi_ads_x,g0_psi_ads_xx,g0_psi_ads_xy
        real*8 g0_psi_ads_y,g0_psi_ads_yy

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_psi_ads0

        real*8 g0u_tt_ads0,g0u_xx_ads0
        real*8 g0u_xy_ads0,g0u_yy_ads0,g0u_psi_ads0

        real*8 Hb_t_t,Hb_t_x,Hb_t_y
        real*8 Hb_x_t,Hb_x_x,Hb_x_y
        real*8 Hb_y_t,Hb_y_x,Hb_y_y

        real*8 Hb_t0,Hb_x0,Hb_y0

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))

        x0=x(i)
        y0=y(j)
        rho0=sqrt(x0**2+y0**2)
       
        f0=(1-rho0**2)**2+4*rho0**2/L**2

        ! set gads values
        g0_tt_ads0 =-f0
     &             /(1-rho0**2)**2
        g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *4
        g0_xy_ads0 =((1+rho0**2)**2/f0-1)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *x0*y0
     &             *4
        g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
     &             /(1-rho0**2)**2
     &             /rho0**2
     &             *4
        g0_psi_ads0=(y0**2)
     &             /(1-rho0**2)**2
     &             *4

        g0u_tt_ads0 =1/g0_tt_ads0
        g0u_xx_ads0 =g0_yy_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
        g0u_xy_ads0 =-g0_xy_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
        g0u_yy_ads0 =g0_xx_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
        g0u_psi_ads0=1/g0_psi_ads0

        ! set gbar values
        gb_tt0=gb_tt_n(i,j)
        gb_tx0=gb_tx_n(i,j)
        gb_ty0=gb_ty_n(i,j)
        gb_xx0=gb_xx_n(i,j)
        gb_xy0=gb_xy_n(i,j)
        gb_yy0=gb_yy_n(i,j)
        psi0=psi_n(i,j)

        ! set hbar values
        Hb_t0=Hb_t_n(i,j)
        Hb_x0=Hb_x_n(i,j)
        Hb_y0=Hb_y_n(i,j)

        ! set phi1 value
        phi10=phi1_n(i,j)

        ! ASSUMES L=1
        g0_tt_ads_x  =(8*x0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_y  =(8*y0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_xx =(-8*(1 + 3*x0**4 - y0**4 + 2*x0**2*(4 + y0**2)))/
     &  (L**2*(-1 + x0**2 + y0**2)**4)
        g0_tt_ads_xy =(-32*x0*y0*(2 + x0**2 + y0**2))/
     &  (L**2*(-1 + x0**2 + y0**2)**4)
        g0_tt_ads_yy =(8*(-1 + x0**4 - 2*(4 + x0**2)*y0**2 - 3*y0**4))/
     &  (L**2*(-1 + x0**2 + y0**2)**4)

        g0_xx_ads_x  =(-16*x0*(8*y0**2*(-1 + 3*x0**2 + 3*y0**2) + 
     &                L**4*(-1 + rho0**2)**2*
     &                (3 + x0**4 - 4*y0**2 + y0**4 + 2*x0**2*(2+y0**2))+
     &                2*L**2*(-1 + x0**6 + 11*y0**2 +5*y0**4*(-3+y0**2)+
     &                x0**4*(5 + 7*y0**2) + 
     &                x0**2*(3 - 10*y0**2 + 11*y0**4))))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_xx_ads_y  =(-16*y0*(8*(-x0**4 + 2*y0**4 + x0**2*(1+y0**2))+
     &                L**4*(-1 + rho0**2)**2*
     &                (x0**4 + (-1 + y0**2)**2 + 2*x0**2*(3 + y0**2)) + 
     &                8*L**2*(y0**2*(-1 + y0**2)**2 + x0**4*(3 + y0**2)+
     &                x0**2*(-1 + y0**2 + 2*y0**4))))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_xx_ads_xx =(16*(32*y0**2*(3*(x0**2 - 4*x0**4 + 7*x0**6) + 
     &                (-1 - 8*x0**2 + 39*x0**4)*y0**2 + 
     &                (4 + 15*x0**2)*y0**4 - 3*y0**6) + 
     &                L**6*(-1 + x0**2 + y0**2)**4*
     &                (5*x0**6 - (-3 + y0**2)*(-1 + y0**2)**2 + 
     &                x0**4*(33 + 9*y0**2) + 
     &                3*x0**2*(13 - 14*y0**2 + y0**4)) + 
     &                2*L**4*(-1 + x0**2 + y0**2)**2*
     &                (9*x0**8 + x0**6*(78 + 60*y0**2) - 
     &                (-1 + y0**2)**2*(1 - 16*y0**2 + 7*y0**4) + 
     &                2*x0**4*(40 - 67*y0**2 + 43*y0**4) + 
     &                2*x0**2*(-11 + 88*y0**2 - 91*y0**4 + 14*y0**6)) + 
     &                8*L**2*(3*x0**10 + 4*x0**8*(7 + 16*y0**2) - 
     &                2*(y0 - y0**3)**2*(1 - 7*y0**2 + 4*y0**4) + 
     &                2*x0**6*(13 - 63*y0**2 + 83*y0**4) + 
     &                6*x0**4*(-2 + 31*y0**2 - 51*y0**4 + 24*y0**6) + 
     &                x0**2*(-1 + y0)*(1 + y0)*
     &                (-3 + 31*y0**2 - 91*y0**4 + 31*y0**6))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_xx_ads_xy =(32*x0*y0*(8*L**2*(1 + 6*x0**2-28*x0**4+50*x0**6 -
     &                5*x0**8 + 2*(-7 + 5*x0**2*(3 + 3*x0**2 + x0**4))*
     &                y0**2 + 2*(29 - 45*x0**2 + 30*x0**4)*y0**4 + 
     &                70*(-1 + x0**2)*y0**6 + 25*y0**8) + 
     &                L**6*(-1 + x0**2 + y0**2)**4*
     &                (11 + 3*x0**4 - 14*y0**2 + 3*y0**4 + 
     &                x0**2*(26 + 6*y0**2)) - 
     &                32*(3*x0**6 - y0**2 + 4*y0**4 - 9*y0**6 - 
     &                x0**4*(4 + 3*y0**2) + x0**2*(1 - 15*y0**4)) + 
     &                4*L**4*(-1 + x0**2 + y0**2)**2*
     &                (-4 + x0**6 + 31*y0**2 - 38*y0**4 + 11*y0**6 + 
     &                x0**4*(42 + 13*y0**2) + 
     &                x0**2*(-3 + 4*y0**2 + 23*y0**4))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_xx_ads_yy =(-16*(L**6*(-1 + x0**2 + y0**2)**4*
     &                (-1 - 5*x0**2 + 5*x0**4 + x0**6 - 
     &                3*(1 + 22*x0**2 + x0**4)*y0**2 - 
     &                9*(-1 + x0**2)*y0**4 - 5*y0**6) + 
     &                4*L**4*(-1 + x0**2 + y0**2)**2*
     &                (x0**8 - 3*y0**2*(-1 + y0**2)**2*(1 + 5*y0**2) + 
     &                x0**6*(11 + 8*y0**2) - 
     &                x0**4*(13 + 111*y0**2 + 2*y0**4) + 
     &                x0**2*(1 + 46*y0**2 - 95*y0**4 - 24*y0**6)) - 
     &                32*(x0**8 - x0**6*(2 + 11*y0**2) + 
     &                x0**4*(1 + 14*y0**2 - 15*y0**4) + 
     &                x0**2*y0**2*(-3 + 18*y0**2 + 7*y0**4) + 
     &                2*(y0**6 + 5*y0**8)) - 
     &                8*L**2*(x0**10 + 
     &                6*y0**4*(-1 + y0**2)**2*(1 + 5*y0**2) - 
     &                2*x0**8*(8 + 13*y0**2) + 
     &                x0**6*(22 + 138*y0**2 - 54*y0**4) + 
     &                2*x0**4*(-4 - 55*y0**2 + 135*y0**4 + 2*y0**6) + 
     &                x0**2*(1 + 38*y0**2 - 114*y0**4 + 62*y0**6 + 
     &                61*y0**8))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)

        g0_xy_ads_x  =(-16*(-1 + L)*(1 + L)*y0*
     &                (L**2*(1 + 7*x0**2 - y0**2)*(-1 + rho0**2)**2 + 
     &                4*(5*x0**4 + y0**2 - y0**4 +x0**2*(-1+4*y0**2))))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_xy_ads_y  =(16*(-1 + L)*(1 + L)*x0*
     &                (L**2*(-1 + x0**2 - 7*y0**2)*(-1 + rho0**2)**2 + 
     &                4*(x0**4 + y0**2 - 5*y0**4 - x0**2*(1+4*y0**2))))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_xy_ads_xx =(-128*x0*y0*(4*(x0**2 - 3*y0**2) - 
     &                L**6*(3+7*x0**2-3*y0**2)*(-1+ x0**2 + y0**2)**4 + 
     &                4*(x0**2 + y0**2)*
     &                (x0**2*(-4 + 15*x0**2) + 6*(2 + x0**2)*y0**2 - 
     &                9*y0**4) + L**4*(-1 + x0**2 + y0**2)**2*
     &                (6 + x0**2 - 50*x0**4 + 7*x0**6 + 
     &                (-33 - 20*x0**2 + 11*x0**4)*y0**2 + 
     &                (30 + x0**2)*y0**4 - 3*y0**6) + 
     &                L**2*(-3+2*x0**2+52*x0**4- 138*x0**6 + 39*x0**8 + 
     &                2*(21 - 34*x0**2 - 87*x0**4 + 48*x0**6)*y0**2 + 
     &                6*(-20 + 11*x0**2 + 9*x0**4)*y0**4 + 
     &                6*(17 - 4*x0**2)*y0**6 - 21*y0**8)))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_xy_ads_xy =(-16*(-(L**4*(-1 + x0**2 + y0**2)**2*
     &                (-1 - 4*x0**2 + 66*x0**4 - 68*x0**6 + 7*x0**8 - 
     &                4*(1 + 35*x0**2 - 109*x0**4 + 13*x0**6)*y0**2 + 
     &                2*(33 + 218*x0**2 - 59*x0**4)*y0**4 - 
     &                4*(17 + 13*x0**2)*y0**6 + 7*y0**8)) + 
     &                L**6*(-1 + x0**2 + y0**2)**4*
     &                (-1 + 7*x0**4 - 6*y0**2 + 7*y0**4 - 
     &                6*x0**2*(1 + 11*y0**2)) + 
     &                16*(-5*x0**8 - y0**4 + 6*y0**6 - 5*y0**8 + 
     &                x0**6*(6 + 28*y0**2) + 
     &                2*x0**2*y0**2*(3 - 7*y0**2 + 14*y0**4) + 
     &                x0**4*(-1 - 14*y0**2 + 66*y0**4)) + 
     &                16*L**2*(-3*x0**10 + x0**8*(14 + 15*y0**2) + 
     &                x0**6*(-15 - 64*y0**2 + 60*y0**4) + 
     &                y0**4*(4 - 15*y0**2 + 14*y0**4 - 3*y0**6) + 
     &                x0**2*y0**2*
     &                (-12 + 41*y0**2 - 64*y0**4 + 15*y0**6) + 
     &                x0**4*(4 + 41*y0**2 - 156*y0**4 + 60*y0**6))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_xy_ads_yy =(-128*x0*y0*(4*(-3*x0**2 + y0**2) + 
     &                L**6*(-3 + 3*x0**2 - 7*y0**2)*
     &                (-1 + x0**2 + y0**2)**4 - 
     &                4*(x0**2 + y0**2)*
     &                (9*x0**4+4*y0**2 - 15*y0**4 - 6*x0**2*(2 + y0**2))
     &                - L**4*(-1 + x0**2 + y0**2)**2*
     &                (-6 + 3*x0**6 - y0**2 + 50*y0**4 - 7*y0**6 - 
     &                x0**4*(30 + y0**2) + 
     &                x0**2*(33 + 20*y0**2 - 11*y0**4)) + 
     &                L**2*(-3-21*x0**8+2*y0**2+ 52*y0**4 - 138*y0**6 + 
     &                39*y0**8 + 6*x0**6*(17 - 4*y0**2) + 
     &                6*x0**4*(-20 + 11*y0**2 + 9*y0**4) + 
     &                2*x0**2*(21 - 34*y0**2 - 87*y0**4 + 48*y0**6))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)

        g0_yy_ads_x  =(-16*x0*(16*x0**4 + 8*(1 + x0**2)*y0**2 -8*y0**4+ 
     &                L**4*(-1 + rho0**2)**2*
     &                ((-1 + x0**2)**2 + 2*(3 + x0**2)*y0**2 + y0**4) + 
     &                8*L**2*((x0 - x0**3)**2 + 
     &                (-1 + x0**2 + 2*x0**4)*y0**2 + (3+x0**2)*y0**4)))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_yy_ads_y  =(-16*y0*(8*x0**2*(-1 + 3*x0**2 + 3*y0**2) + 
     &                L**4*(-1 + rho0**2)**2*
     &                (3 + x0**4 + 4*y0**2 + y0**4 +2*x0**2*(-2+y0**2))
     &                + 2*L**2*(-1 + 5*x0**6 + 3*y0**2 + 5*y0**4 +y0**6+
     &                x0**4*(-15 + 11*y0**2) + 
     &                x0**2*(11 - 10*y0**2 + 7*y0**4))))/
     &                ((-1 + rho0**2)**3*
     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
        g0_yy_ads_xx =(16*(4*L**4*(-1 + x0**2 + y0**2)**2*
     &                (3*x0**2*(-1 + x0**2)**2*(1 + 5*x0**2) + 
     &                (-1 - 46*x0**2 + 95*x0**4 + 24*x0**6)*y0**2 + 
     &                (13 + 111*x0**2 + 2*x0**4)*y0**4 - 
     &                (11 + 8*x0**2)*y0**6 - y0**8) + 
     &                8*L**2*(6*x0**4*(-1 + x0**2)**2*(1 + 5*x0**2) + 
     &                (1 + 38*x0**2 - 114*x0**4 + 62*x0**6 + 61*x0**8)*
     &                y0**2 + 2*(-4 - 55*x0**2 + 135*x0**4 + 2*x0**6)*
     &                y0**4 + 2*(11 + 69*x0**2 - 27*x0**4)*y0**6 - 
     &                2*(8 + 13*x0**2)*y0**8 + y0**10) + 
     &                32*(10*x0**8 + y0**4 - 2*y0**6 + y0**8 + 
     &                3*x0**4*y0**2*(6-5*y0**2) + x0**6*(2 + 7*y0**2) + 
     &                x0**2*y0**2*(-3 + 14*y0**2 - 11*y0**4)) + 
     &                L**6*(-1 + x0**2 + y0**2)**4*
     &                (1 + 5*x0**6 + 5*y0**2 + 9*x0**4*(-1 + y0**2) - 
     &                y0**4*(5+y0**2) + 3*x0**2*(1 + 22*y0**2 + y0**4)))
     &                )/((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_yy_ads_xy =(32*x0*y0*(32*(x0**2 - 4*x0**4 + 9*x0**6 + 
     &                (-1 + 15*x0**4)*y0**2 + (4 + 3*x0**2)*y0**4 - 
     &                3*y0**6) + L**6*(-1 + x0**2 + y0**2)**4*
     &                (11 + 3*x0**4 + 26*y0**2 + 3*y0**4 + 
     &                2*x0**2*(-7 + 3*y0**2)) + 
     &                4*L**4*(-1 + x0**2 + y0**2)**2*
     &                (-4 + 11*x0**6 - 3*y0**2 + 42*y0**4 + y0**6 + 
     &                x0**4*(-38 + 23*y0**2) + 
     &                x0**2*(31 + 4*y0**2 + 13*y0**4)) + 
     &                8*L**2*(1+25*x0**8+6*y0**2- 28*y0**4 + 50*y0**6 - 
     &                5*y0**8 + 70*x0**6*(-1 + y0**2) + 
     &                x0**4*(58 - 90*y0**2 + 60*y0**4) + 
     &                2*x0**2*(-7 + 5*y0**2*(3 + 3*y0**2 + y0**4)))))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
        g0_yy_ads_yy =(-16*(32*x0**2*(x0**2 - 4*x0**4 + 3*x0**6 + 
     &                (-3 + 8*x0**2 - 15*x0**4)*y0**2 + 
     &                3*(4 - 13*x0**2)*y0**4 - 21*y0**6) + 
     &                L**6*(-1 + x0**2 + y0**2)**4*
     &                ((-3 + x0**2)*(-1 + x0**2)**2 - 
     &                3*(13 - 14*x0**2 + x0**4)*y0**2 - 
     &                3*(11 + 3*x0**2)*y0**4 - 5*y0**6) + 
     &                2*L**4*(-1 + x0**2 + y0**2)**2*
     &                ((-1 + x0**2)**2*(1 - 16*x0**2 + 7*x0**4) + 
     &                2*(11 - 88*x0**2 + 91*x0**4 - 14*x0**6)*y0**2 - 
     &                2*(40 - 67*x0**2 + 43*x0**4)*y0**4 - 
     &                6*(13 + 10*x0**2)*y0**6 - 9*y0**8) + 
     &                8*L**2*(2*(x0-x0**3)**2*(1 - 7*x0**2 + 4*x0**4) + 
     &                (-3+ 34*x0**2 - 122*x0**4 + 122*x0**6 - 31*x0**8)*
     &                y0**2 - 6*(-2 + 31*x0**2 - 51*x0**4 + 24*x0**6)*
     &                y0**4 - 2*(13 - 63*x0**2 + 83*x0**4)*y0**6 - 
     &                4*(7 + 16*x0**2)*y0**8 - 3*y0**10)))/
     &                ((-1 + x0**2 + y0**2)**4*
     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)

        g0_psi_ads_x  =(-16*x0*y0**2)/(-1 + rho0**2)**3
        g0_psi_ads_y  =(-8*(y0 - x0**2*y0 + y0**3))/(-1 + rho0**2)**3
        g0_psi_ads_xx =(-16*y0**2*(-1 - 5*x0**2 + y0**2))/
     &                 (-1 + x0**2 + y0**2)**4
        g0_psi_ads_xy =(-32*x0*y0*(-1 + x0**2 - 2*y0**2))/
     &                 (-1 + x0**2 + y0**2)**4
        g0_psi_ads_yy =(8*((-1+x0**2)**2-8*(-1+x0**2)*y0**2+3*y0**4))/
     &                 (-1 + x0**2 + y0**2)**4

        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t,
     &       gb_tt_x,gb_tt_y,gb_tt_tt,gb_tt_tx,gb_tt_ty,
     &       gb_tt_xx,gb_tt_xy,gb_tt_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,'gb_tt')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t,
     &       gb_tx_x,gb_tx_y,gb_tx_tt,gb_tx_tx,gb_tx_ty,
     &       gb_tx_xx,gb_tx_xy,gb_tx_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,'gb_tx')
        call df2_int(gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t,
     &       gb_ty_x,gb_ty_y,gb_ty_tt,gb_ty_tx,gb_ty_ty,
     &       gb_ty_xx,gb_ty_xy,gb_ty_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,'gb_ty')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t,
     &       gb_xx_x,gb_xx_y,gb_xx_tt,gb_xx_tx,gb_xx_ty,
     &       gb_xx_xx,gb_xx_xy,gb_xx_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,'gb_xx')
        call df2_int(gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t,
     &       gb_xy_x,gb_xy_y,gb_xy_tt,gb_xy_tx,gb_xy_ty,
     &       gb_xy_xx,gb_xy_xy,gb_xy_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,'gb_xy')
        call df2_int(gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t,
     &       gb_yy_x,gb_yy_y,gb_yy_tt,gb_yy_tx,gb_yy_ty,
     &       gb_yy_xx,gb_yy_xy,gb_yy_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,"gb_yy")
        call df2_int(psi_np1,psi_n,psi_nm1,psi_t,psi_x,
     &       psi_y,psi_tt,psi_tx,psi_ty,psi_xx,psi_xy,psi_yy,
     &       x,y,dt,i,j,chr,ex,Nx,Ny,"psi")

        ! calculate hbar derivatives
        call df1_int(Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t,Hb_t_x,
     &       Hb_t_y,x,y,dt,i,j,
     &       chr,ex,Nx,Ny,'Hb_t')
        call df1_int(Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t,Hb_x_x,
     &       Hb_x_y,x,y,dt,i,j,
     &       chr,ex,Nx,Ny,'Hb_x')
        call df1_int(Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t,Hb_y_x,
     &       Hb_y_y,x,y,dt,i,j,
     &       chr,ex,Nx,Ny,'Hb_y')

        ! calculate phi1 derivatives
        call df2_int(phi1_np1,phi1_n,phi1_nm1,phi1_t,phi1_x,
     &       phi1_y,phi1_tt,phi1_tx,phi1_ty,phi1_xx,
     &       phi1_xy,phi1_yy,x,y,dt,i,j,
     &       chr,ex,Nx,Ny,'phi1')

        ! give values to the metric
        g0_ll(1,1)=g0_tt_ads0+gb_tt0
        g0_ll(1,2)=           gb_tx0
        g0_ll(1,3)=           gb_ty0
        g0_ll(2,2)=g0_xx_ads0+gb_xx0
        g0_ll(2,3)=g0_xy_ads0+gb_xy0
        g0_ll(3,3)=g0_yy_ads0+gb_yy0
        g0_ll(4,4)=g0_psi_ads0+psi0*y0**2

        g0_ll_x(1,1,1)   =0
     &                   +gb_tt_t
        g0_ll_x(1,1,2)   =g0_tt_ads_x
     &                   +gb_tt_x
        g0_ll_x(1,1,3)   =g0_tt_ads_y
     &                   +gb_tt_y
        g0_ll_xx(1,1,1,1)=0
     &                   +gb_tt_tt
        g0_ll_xx(1,1,1,2)=0
     &                   +gb_tt_tx
        g0_ll_xx(1,1,1,3)=0
     &                   +gb_tt_ty
        g0_ll_xx(1,1,2,2)=g0_tt_ads_xx
     &                   +gb_tt_xx
        g0_ll_xx(1,1,2,3)=g0_tt_ads_xy
     &                   +gb_tt_xy
        g0_ll_xx(1,1,3,3)=g0_tt_ads_yy
     &                   +gb_tt_yy

        g0_ll_x(1,2,1)   =0
     &                   +gb_tx_t
        g0_ll_x(1,2,2)   =
     &                   +gb_tx_x
        g0_ll_x(1,2,3)   =
     &                   +gb_tx_y
        g0_ll_xx(1,2,1,1)=0
     &                   +gb_tx_tt
        g0_ll_xx(1,2,1,2)=0
     &                   +gb_tx_tx
        g0_ll_xx(1,2,1,3)=0
     &                   +gb_tx_ty
        g0_ll_xx(1,2,2,2)=
     &                   +gb_tx_xx
        g0_ll_xx(1,2,2,3)=
     &                   +gb_tx_xy
        g0_ll_xx(1,2,3,3)=
     &                   +gb_tx_yy

        g0_ll_x(1,3,1)   =0
     &                   +gb_ty_t
        g0_ll_x(1,3,2)   =
     &                   +gb_ty_x
        g0_ll_x(1,3,3)   =
     &                   +gb_ty_y
        g0_ll_xx(1,3,1,1)=0
     &                   +gb_ty_tt
        g0_ll_xx(1,3,1,2)=0
     &                   +gb_ty_tx
        g0_ll_xx(1,3,1,3)=0
     &                   +gb_ty_ty
        g0_ll_xx(1,3,2,2)=
     &                   +gb_ty_xx
        g0_ll_xx(1,3,2,3)=
     &                   +gb_ty_xy
        g0_ll_xx(1,3,3,3)=
     &                   +gb_ty_yy

        g0_ll_x(2,2,1)   =0
     &                   +gb_xx_t
        g0_ll_x(2,2,2)   =g0_xx_ads_x
     &                   +gb_xx_x
        g0_ll_x(2,2,3)   =g0_xx_ads_y
     &                   +gb_xx_y
        g0_ll_xx(2,2,1,1)=0
     &                   +gb_xx_tt
        g0_ll_xx(2,2,1,2)=0
     &                   +gb_xx_tx
        g0_ll_xx(2,2,1,3)=0
     &                   +gb_xx_ty
        g0_ll_xx(2,2,2,2)=g0_xx_ads_xx
     &                   +gb_xx_xx
        g0_ll_xx(2,2,2,3)=g0_xx_ads_xy
     &                   +gb_xx_xy
        g0_ll_xx(2,2,3,3)=g0_xx_ads_yy
     &                   +gb_xx_yy

        g0_ll_x(2,3,1)   =0
     &                   +gb_xy_t
        g0_ll_x(2,3,2)   =g0_xy_ads_x
     &                   +gb_xy_x
        g0_ll_x(2,3,3)   =g0_xy_ads_y
     &                   +gb_xy_y
        g0_ll_xx(2,3,1,1)=0
     &                   +gb_xy_tt
        g0_ll_xx(2,3,1,2)=0
     &                   +gb_xy_tx
        g0_ll_xx(2,3,1,3)=0
     &                   +gb_xy_ty
        g0_ll_xx(2,3,2,2)=g0_xy_ads_xx
     &                   +gb_xy_xx
        g0_ll_xx(2,3,2,3)=g0_xy_ads_xy
     &                   +gb_xy_xy
        g0_ll_xx(2,3,3,3)=g0_xy_ads_yy
     &                   +gb_xy_yy

        g0_ll_x(3,3,1)   =0
     &                   +gb_yy_t
        g0_ll_x(3,3,2)   =g0_yy_ads_x
     &                   +gb_yy_x
        g0_ll_x(3,3,3)   =g0_yy_ads_y
     &                   +gb_yy_y
        g0_ll_xx(3,3,1,1)=0
     &                   +gb_yy_tt
        g0_ll_xx(3,3,1,2)=0
     &                   +gb_yy_tx
        g0_ll_xx(3,3,1,3)=0
     &                   +gb_yy_ty
        g0_ll_xx(3,3,2,2)=g0_yy_ads_xx
     &                   +gb_yy_xx
        g0_ll_xx(3,3,2,3)=g0_yy_ads_xy
     &                   +gb_yy_xy
        g0_ll_xx(3,3,3,3)=g0_yy_ads_yy
     &                   +gb_yy_yy

        g0_ll_x(4,4,1)   =0
     &                   +psi_t*y0**2
        g0_ll_x(4,4,2)   =g0_psi_ads_x
     &                   +psi_x*y0**2
        g0_ll_x(4,4,3)   =g0_psi_ads_y
     &                   +psi_y*y0**2
     &                   +psi0*(2*y0)
        g0_ll_xx(4,4,1,1)=0
     &                   +psi_tt*y0**2
        g0_ll_xx(4,4,1,2)=0
     &                   +psi_tx*y0**2
        g0_ll_xx(4,4,1,3)=0
     &                   +psi_ty*y0**2
     &                   +psi_t*(2*y0)
        g0_ll_xx(4,4,2,2)=g0_psi_ads_xx
     &                   +psi_xx*y0**2
        g0_ll_xx(4,4,2,3)=g0_psi_ads_xy
     &                   +psi_xy*y0**2
     &                   +psi_x*(2*y0)
        g0_ll_xx(4,4,3,3)=g0_psi_ads_yy
     &                   +psi_yy*y0**2
     &                   +psi_y*(2*y0)
     &                   +psi_y*(2*y0)
     &                   +psi0*(2)

        ! give values to the metric inverse
        call calc_g0uu(g0_ll(1,1),g0_ll(1,2),g0_ll(1,3),g0_ll(1,4),
     &         g0_ll(2,2),g0_ll(2,3),g0_ll(2,4),
     &         g0_ll(3,3),g0_ll(3,4),g0_ll(4,4),g0_uu)

        do a=1,3
          do b=a+1,4
            g0_ll(b,a)=g0_ll(a,b)
            g0_uu(b,a)=g0_uu(a,b) 
            do c=1,4
              g0_ll_x(b,a,c)=g0_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              g0_uu_x(a,b,c)=0
              do d=1,4
                do e=1,4
                  g0_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                          -g0_ll_x(d,e,c)
     &                           *g0_uu(a,d)*g0_uu(b,e)
                end do
     &  
              end do
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                g0_ll_xx(a,b,c,d)=
     &             g0_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

        ! give values to the Christoffel symbols
        do a=1,4
          do b=1,4
            do c=1,4
              gamma_ull(a,b,c)=0
              do d=1,4
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g0_uu(a,d)
     &                                *(g0_ll_x(c,d,b)
     &                                 -g0_ll_x(b,c,d)
     &                                 +g0_ll_x(d,b,c))
              end do
            end do
          end do
        end do

        ! calculate Christoffel symbol derivatives at point i,j
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        do a=1,4
          do b=1,4
            do c=1,4
              do e=1,4
                gamma_ull_x(a,b,c,e)=0
                do d=1,4
                  gamma_ull_x(a,b,c,e)=gamma_ull_x(a,b,c,e)
     &              +0.5d0*g0_uu_x(a,d,e)*(g0_ll_x(b,d,c)+
     &                     g0_ll_x(c,d,b)-g0_ll_x(b,c,d))
     &              +0.5d0*g0_uu(a,d)*(g0_ll_xx(b,d,c,e)+
     &                     g0_ll_xx(c,d,b,e)-g0_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do

        ! give values to the AdS metric
        gads_ll(1,1)=g0_tt_ads0
        gads_ll(2,2)=g0_xx_ads0
        gads_ll(2,3)=g0_xy_ads0
        gads_ll(3,3)=g0_yy_ads0
        gads_ll(4,4)=g0_psi_ads0

        gads_uu(1,1)=g0u_tt_ads0
        gads_uu(2,2)=g0u_xx_ads0
        gads_uu(2,3)=g0u_xy_ads0
        gads_uu(3,3)=g0u_yy_ads0
        gads_uu(4,4)=g0u_psi_ads0

        gads_ll_x(1,1,2)   =g0_tt_ads_x
        gads_ll_x(1,1,3)   =g0_tt_ads_y
        gads_ll_xx(1,1,2,2)=g0_tt_ads_xx
        gads_ll_xx(1,1,2,3)=g0_tt_ads_xy
        gads_ll_xx(1,1,3,3)=g0_tt_ads_yy
        gads_ll_x(2,2,2)   =g0_xx_ads_x
        gads_ll_x(2,2,3)   =g0_xx_ads_y
        gads_ll_xx(2,2,2,2)=g0_xx_ads_xx
        gads_ll_xx(2,2,2,3)=g0_xx_ads_xy
        gads_ll_xx(2,2,3,3)=g0_xx_ads_yy
        gads_ll_x(2,3,2)   =g0_xy_ads_x
        gads_ll_x(2,3,3)   =g0_xy_ads_y
        gads_ll_xx(2,3,2,2)=g0_xy_ads_xx
        gads_ll_xx(2,3,2,3)=g0_xy_ads_xy
        gads_ll_xx(2,3,3,3)=g0_xy_ads_yy
        gads_ll_x(3,3,2)   =g0_yy_ads_x
        gads_ll_x(3,3,3)   =g0_yy_ads_y
        gads_ll_xx(3,3,2,2)=g0_yy_ads_xx
        gads_ll_xx(3,3,2,3)=g0_yy_ads_xy
        gads_ll_xx(3,3,3,3)=g0_yy_ads_yy
        gads_ll_x(4,4,2)   =g0_psi_ads_x
        gads_ll_x(4,4,3)   =g0_psi_ads_y
        gads_ll_xx(4,4,2,2)=g0_psi_ads_xx
        gads_ll_xx(4,4,2,3)=g0_psi_ads_xy
        gads_ll_xx(4,4,3,3)=g0_psi_ads_yy

        do a=1,3
          do b=a+1,4
            gads_ll(b,a)=gads_ll(a,b)
            gads_uu(b,a)=gads_uu(a,b)
            do c=1,4
              gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              gads_uu_x(a,b,c)=
     &              -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &              -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                               +gads_uu(a,2)*gads_uu(b,1))
     &              -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,1))
     &              -gads_ll_x(1,4,c)*(gads_uu(a,1)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,1))
     &              -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &              -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,2))
     &              -gads_ll_x(2,4,c)*(gads_uu(a,2)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,2))
     &              -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
     &              -gads_ll_x(3,4,c)*(gads_uu(a,3)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,3))
     &              -gads_ll_x(4,4,c)*gads_uu(a,4)*gads_uu(b,4)
            end do
          end do
        end do

        ! give values to the metric deviation
        h0_ll(1,1)=gb_tt0 
        h0_ll(1,2)=gb_tx0
        h0_ll(1,3)=gb_ty0
        h0_ll(2,2)=gb_xx0
        h0_ll(2,3)=gb_xy0
        h0_ll(3,3)=gb_yy0
        h0_ll(4,4)=psi0*y0**2
  
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(1,3)=g0_uu(1,3)
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(2,3)=g0_uu(2,3)-gads_uu(2,3)
        h0_uu(3,3)=g0_uu(3,3)-gads_uu(3,3)
        h0_uu(4,4)=g0_uu(4,4)-gads_uu(4,4)

        h0_ll_x(1,1,1)   =g0_ll_x(1,1,1)
        h0_ll_x(1,1,2)   =g0_ll_x(1,1,2)-gads_ll_x(1,1,2)
        h0_ll_x(1,1,3)   =g0_ll_x(1,1,3)-gads_ll_x(1,1,3)
        h0_ll_xx(1,1,1,1)=g0_ll_xx(1,1,1,1)
        h0_ll_xx(1,1,1,2)=g0_ll_xx(1,1,1,2)
        h0_ll_xx(1,1,1,3)=g0_ll_xx(1,1,1,3)
        h0_ll_xx(1,1,2,2)=g0_ll_xx(1,1,2,2)-gads_ll_xx(1,1,2,2)
        h0_ll_xx(1,1,2,3)=g0_ll_xx(1,1,2,3)-gads_ll_xx(1,1,2,3)
        h0_ll_xx(1,1,3,3)=g0_ll_xx(1,1,3,3)-gads_ll_xx(1,1,3,3)

        h0_ll_x(1,2,1)   =g0_ll_x(1,2,1)
        h0_ll_x(1,2,2)   =g0_ll_x(1,2,2)
        h0_ll_x(1,2,3)   =g0_ll_x(1,2,3)
        h0_ll_xx(1,2,1,1)=g0_ll_xx(1,2,1,1)
        h0_ll_xx(1,2,1,2)=g0_ll_xx(1,2,1,2)
        h0_ll_xx(1,2,1,3)=g0_ll_xx(1,2,1,3)
        h0_ll_xx(1,2,2,2)=g0_ll_xx(1,2,2,2)
        h0_ll_xx(1,2,2,3)=g0_ll_xx(1,2,2,3)
        h0_ll_xx(1,2,3,3)=g0_ll_xx(1,2,3,3)

        h0_ll_x(1,3,1)   =g0_ll_x(1,3,1)
        h0_ll_x(1,3,2)   =g0_ll_x(1,3,2)
        h0_ll_x(1,3,3)   =g0_ll_x(1,3,3)
        h0_ll_xx(1,3,1,1)=g0_ll_xx(1,3,1,1)
        h0_ll_xx(1,3,1,2)=g0_ll_xx(1,3,1,2)
        h0_ll_xx(1,3,1,3)=g0_ll_xx(1,3,1,3)
        h0_ll_xx(1,3,2,2)=g0_ll_xx(1,3,2,2)
        h0_ll_xx(1,3,2,3)=g0_ll_xx(1,3,2,3)
        h0_ll_xx(1,3,3,3)=g0_ll_xx(1,3,3,3)

        h0_ll_x(2,2,1)   =g0_ll_x(2,2,1)
        h0_ll_x(2,2,2)   =g0_ll_x(2,2,2)-gads_ll_x(2,2,2)
        h0_ll_x(2,2,3)   =g0_ll_x(2,2,3)-gads_ll_x(2,2,3)
        h0_ll_xx(2,2,1,1)=g0_ll_xx(2,2,1,1)
        h0_ll_xx(2,2,1,2)=g0_ll_xx(2,2,1,2)
        h0_ll_xx(2,2,1,3)=g0_ll_xx(2,2,1,3)
        h0_ll_xx(2,2,2,2)=g0_ll_xx(2,2,2,2)-gads_ll_xx(2,2,2,2)
        h0_ll_xx(2,2,2,3)=g0_ll_xx(2,2,2,3)-gads_ll_xx(2,2,2,3)
        h0_ll_xx(2,2,3,3)=g0_ll_xx(2,2,3,3)-gads_ll_xx(2,2,3,3)

        h0_ll_x(2,3,1)   =g0_ll_x(2,3,1)
        h0_ll_x(2,3,2)   =g0_ll_x(2,3,2)-gads_ll_x(2,3,2)
        h0_ll_x(2,3,3)   =g0_ll_x(2,3,3)-gads_ll_x(2,3,3)
        h0_ll_xx(2,3,1,1)=g0_ll_xx(2,3,1,1)
        h0_ll_xx(2,3,1,2)=g0_ll_xx(2,3,1,2)
        h0_ll_xx(2,3,1,3)=g0_ll_xx(2,3,1,3)
        h0_ll_xx(2,3,2,2)=g0_ll_xx(2,3,2,2)-gads_ll_xx(2,3,2,2)
        h0_ll_xx(2,3,2,3)=g0_ll_xx(2,3,2,3)-gads_ll_xx(2,3,2,3)
        h0_ll_xx(2,3,3,3)=g0_ll_xx(2,3,3,3)-gads_ll_xx(2,3,3,3)

        h0_ll_x(3,3,1)   =g0_ll_x(3,3,1)
        h0_ll_x(3,3,2)   =g0_ll_x(3,3,2)-gads_ll_x(3,3,2)
        h0_ll_x(3,3,3)   =g0_ll_x(3,3,3)-gads_ll_x(3,3,3)
        h0_ll_xx(3,3,1,1)=g0_ll_xx(3,3,1,1)
        h0_ll_xx(3,3,1,2)=g0_ll_xx(3,3,1,2)
        h0_ll_xx(3,3,1,3)=g0_ll_xx(3,3,1,3)
        h0_ll_xx(3,3,2,2)=g0_ll_xx(3,3,2,2)-gads_ll_xx(3,3,2,2)
        h0_ll_xx(3,3,2,3)=g0_ll_xx(3,3,2,3)-gads_ll_xx(3,3,2,3)
        h0_ll_xx(3,3,3,3)=g0_ll_xx(3,3,3,3)-gads_ll_xx(3,3,3,3)

        h0_ll_x(4,4,1)   =g0_ll_x(4,4,1)
        h0_ll_x(4,4,2)   =g0_ll_x(4,4,2)-gads_ll_x(4,4,2)
        h0_ll_x(4,4,3)   =g0_ll_x(4,4,3)-gads_ll_x(4,4,3)
        h0_ll_xx(4,4,1,1)=g0_ll_xx(4,4,1,1)
        h0_ll_xx(4,4,1,2)=g0_ll_xx(4,4,1,2)
        h0_ll_xx(4,4,1,3)=g0_ll_xx(4,4,1,3)
        h0_ll_xx(4,4,2,2)=g0_ll_xx(4,4,2,2)-gads_ll_xx(4,4,2,2)
        h0_ll_xx(4,4,2,3)=g0_ll_xx(4,4,2,3)-gads_ll_xx(4,4,2,3)
        h0_ll_xx(4,4,3,3)=g0_ll_xx(4,4,3,3)-gads_ll_xx(4,4,3,3)

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)
        h0_uu_x(1,1,3)=g0_uu_x(1,1,3)-gads_uu_x(1,1,3)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)
        h0_uu_x(1,2,3)=g0_uu_x(1,2,3)

        h0_uu_x(1,3,1)=g0_uu_x(1,3,1) 
        h0_uu_x(1,3,2)=g0_uu_x(1,3,2) 
        h0_uu_x(1,3,3)=g0_uu_x(1,3,3) 

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)
        h0_uu_x(2,2,3)=g0_uu_x(2,2,3)-gads_uu_x(2,2,3)

        h0_uu_x(2,3,1)=g0_uu_x(2,3,1)
        h0_uu_x(2,3,2)=g0_uu_x(2,3,2)-gads_uu_x(2,3,2)
        h0_uu_x(2,3,3)=g0_uu_x(2,3,3)-gads_uu_x(2,3,3)

        h0_uu_x(3,3,1)=g0_uu_x(3,3,1)
        h0_uu_x(3,3,2)=g0_uu_x(3,3,2)-gads_uu_x(3,3,2)
        h0_uu_x(3,3,3)=g0_uu_x(3,3,3)-gads_uu_x(3,3,3)

        h0_uu_x(4,4,1)=g0_uu_x(4,4,1)
        h0_uu_x(4,4,2)=g0_uu_x(4,4,2)-gads_uu_x(4,4,2)
        h0_uu_x(4,4,3)=g0_uu_x(4,4,3)-gads_uu_x(4,4,3)

        do a=1,3
          do b=a+1,4
            h0_ll(b,a)=h0_ll(a,b)
            h0_uu(b,a)=h0_uu(a,b)
            do c=1,4
              h0_ll_x(b,a,c)=h0_ll_x(a,b,c)
              h0_uu_x(b,a,c)=h0_uu_x(a,b,c)
            end do
          end do
        end do

        ! give values to the source functions
        A_l(1)=Hb_t0*(1-rho0**2)
        A_l(2)=Hb_x0*(1-rho0**2)
        A_l(3)=Hb_y0*(1-rho0**2)

        A_l_x(1,1)=Hb_t_t*(1-rho0**2)
        A_l_x(1,2)=Hb_t_x*(1-rho0**2)
     &            -Hb_t0*2*x0
        A_l_x(1,3)=Hb_t_y*(1-rho0**2)
     &            -Hb_t0*2*y0

        A_l_x(2,1)=Hb_x_t*(1-rho0**2)
        A_l_x(2,2)=Hb_x_x*(1-rho0**2)
     &            -Hb_x0*2*x0
        A_l_x(2,3)=Hb_x_y*(1-rho0**2)
     &            -Hb_x0*2*y0

        A_l_x(3,1)=Hb_y_t*(1-rho0**2)
        A_l_x(3,2)=Hb_y_x*(1-rho0**2)
     &            -Hb_y0*2*x0
        A_l_x(3,3)=Hb_y_y*(1-rho0**2)
     &            -Hb_y0*2*y0

        ! give values to the AdS source functions
        Hads_l(1)=0
        Hads_l(2)=(-2*x0*(3 + rho0**2))/(-1 + (rho0**2)**2)
        Hads_l(3)=1/y0 - (2*y0*(3 + rho0**2))/(-1 + (rho0**2)**2)
        Hads_l(4)=0

        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-rho0**2)**2
        phi10_x(2)=phi1_x*(1-rho0**2)**2
     &            +phi10*(-4*x0)*(1-rho0**2)
        phi10_x(3)=phi1_y*(1-rho0**2)**2
     &            +phi10*(-4*y0)*(1-rho0**2)

        phi10_xx(1,1)=phi1_tt*(1-rho0**2)**2
        phi10_xx(1,2)=phi1_tx*(1-rho0**2)**2
     &               +phi1_t*(-4*x0)*(1-rho0**2)
        phi10_xx(1,3)=phi1_ty*(1-rho0**2)**2
     &               +phi1_t*(-4*y0)*(1-rho0**2)
        phi10_xx(2,2)=phi1_xx*(1-rho0**2)**2
     &               +phi1_x*(2)*(-4*x0)*(1-rho0**2)
     &               +phi10*(-4*(1-rho0**2)+8*x0**2)
        phi10_xx(2,3)=phi1_xy*(1-rho0**2)**2
     &               +phi1_x*(-4*y0)*(1-rho0**2)
     &               +phi1_y*(-4*x0)*(1-rho0**2)
     &               +phi10*(-4*x0)*(-2*y0)
        phi10_xx(3,3)=phi1_yy*(1-rho0**2)**2
     &               +phi1_y*(2)*(-4*y0)*(1-rho0**2)
     &               +phi10*(-4*(1-rho0**2)+8*y0**2)

        do a=1,3
          do b=a+1,4
            phi10_xx(b,a)=phi10_xx(a,b)
          end do
        end do

        ! calculate Riemann tensor at point i,j
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                riemann_ulll(a,b,c,d)=
     &                gamma_ull_x(a,b,d,c)-gamma_ull_x(a,b,c,d)
                do e=1,4
                   riemann_ulll(a,b,c,d)=riemann_ulll(a,b,c,d)
     &               +gamma_ull(a,c,e)*gamma_ull(e,b,d)
     &               -gamma_ull(a,d,e)*gamma_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i,j
        !(R_bd = R^a_bad)
        do b=1,4
          do d=1,4
            ricci_ll(b,d)=0
            do a=1,4
              ricci_ll(b,d)=ricci_ll(b,d)+riemann_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i,j
        !(R_a^b = R_ad g^db)
        do a=1,4
          do b=1,4
            ricci_lu(a,b)=0
            do d=1,4
              ricci_lu(a,b)=ricci_lu(a,b)+ricci_ll(a,d)*g0_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        ricci=0
        do a=1,4
          ricci=ricci+ricci_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i,j
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,4
          do b=1,4
            einstein_ll(a,b)=ricci_ll(a,b)-0.5d0*ricci*g0_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i,j 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1_sq=0
        do a=1,4
          do b=1,4
            grad_phi1_sq=grad_phi1_sq
     &                  +phi10_x(a)*phi10_x(b)*g0_uu(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            set_ll(a,b)=
     &            phi10_x(a)*phi10_x(b)
     &           -g0_ll(a,b)*(grad_phi1_sq/2)
          end do
        end do

        return
        end

