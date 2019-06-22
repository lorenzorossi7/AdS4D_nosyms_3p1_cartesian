c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1. The tensor components are given in
c spherical polar coordinates.
c----------------------------------------------------------------------
!
!        subroutine quasiset(
!     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
!     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
!     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
!     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
!     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
!     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
!     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
!     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
!     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
!     &                  psi_np1,psi_n,psi_nm1,
!     &                  quasiset_tt,quasiset_tchi,quasiset_txi,
!     &                  quasiset_chichi,quasiset_chixi,
!     &                  quasiset_xixi,quasiset_mass,
!     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)
!
!        implicit none
!
!        integer Nx,Ny,Nz
!        integer phys_bdy(6),ghost_width(6)
!        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
!        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
!        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
!        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
!        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
!        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
!        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
!        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
!        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
!        real*8 psi_np1(Nx,Ny,Nz),psi_n(Nx,Ny,Nz),psi_nm1(Nx,Ny,Nz)
!        real*8 L
!        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex
!
!        real*8 chrql(Nx,Ny,Nz)
!
!        real*8 quasiset_tt(Nx,Ny,Nz),quasiset_tchi(Nx,Ny,Nz)
!        real*8 quasiset_txi(Nx,Ny,Nz),quasiset_chichi(Nx,Ny,Nz)
!        real*8 quasiset_chixi(Nx,Ny,Nz),quasiset_xixi(Nx,Ny,Nz)
!        real*8 quasiset_mass(Nx,Ny,Nz)
!
!        integer i,j,k,is,ie,js,je,ks,ke
!
!        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
!        integer ix,jy,kz
!
!        real*8 x0,y0,z0,rho0,q
!
!        real*8 dx,dy,dz
!        integer numbdypoints
!
!        real*8 PI
!        parameter (PI=3.141592653589793d0)
!
!
!
!        !--------------------------------------------------------------
!        ! the following are first and second time derivatives of *n*
!        ! level variables, and as these are the only derivatives we
!        ! use we drop any _n identifier
!        !--------------------------------------------------------------
!        real*8 gb_tt_t, gb_tt_x, gb_tt_y
!        real*8 gb_tt_z
!        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty
!        real*8 gb_tt_tz
!        real*8 gb_tt_xx,gb_tt_xy
!        real*8 gb_tt_xz
!        real*8 gb_tt_yy
!        real*8 gb_tt_yz
!        real*8 gb_tt_zz
!
!        real*8 gb_tx_t, gb_tx_x, gb_tx_y
!        real*8 gb_tx_z
!        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty
!        real*8 gb_tx_tz
!        real*8 gb_tx_xx,gb_tx_xy
!        real*8 gb_tx_xz
!        real*8 gb_tx_yy
!        real*8 gb_tx_yz
!        real*8 gb_tx_zz
!
!        real*8 gb_ty_t, gb_ty_x, gb_ty_y
!        real*8 gb_ty_z
!        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty
!        real*8 gb_ty_tz
!        real*8 gb_ty_xx,gb_ty_xy
!        real*8 gb_ty_xz
!        real*8 gb_ty_yy
!        real*8 gb_ty_yz
!        real*8 gb_ty_zz
!
!        real*8 gb_tz_t, gb_tz_x, gb_tz_y
!        real*8 gb_tz_z
!        real*8 gb_tz_tt,gb_tz_tx,gb_tz_ty
!        real*8 gb_tz_tz
!        real*8 gb_tz_xx,gb_tz_xy
!        real*8 gb_tz_xz
!        real*8 gb_tz_yy
!        real*8 gb_tz_yz
!        real*8 gb_tz_zz
!
!        real*8 gb_xx_t, gb_xx_x, gb_xx_y
!        real*8 gb_xx_z
!        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty
!        real*8 gb_xx_tz
!        real*8 gb_xx_xx,gb_xx_xy
!        real*8 gb_xx_xz
!        real*8 gb_xx_yy
!        real*8 gb_xx_yz
!        real*8 gb_xx_zz
!
!        real*8 gb_xy_t, gb_xy_x, gb_xy_y
!        real*8 gb_xy_z
!        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty
!        real*8 gb_xy_tz
!        real*8 gb_xy_xx,gb_xy_xy
!        real*8 gb_xy_xz
!        real*8 gb_xy_yy
!        real*8 gb_xy_yz
!        real*8 gb_xy_zz
!
!        real*8 gb_xz_t, gb_xz_x, gb_xz_y
!        real*8 gb_xz_z
!        real*8 gb_xz_tt,gb_xz_tx,gb_xz_ty
!        real*8 gb_xz_tz
!        real*8 gb_xz_xx,gb_xz_xy
!        real*8 gb_xz_xz
!        real*8 gb_xz_yy
!        real*8 gb_xz_yz
!        real*8 gb_xz_zz
!
!        real*8 gb_yy_t, gb_yy_x, gb_yy_y
!        real*8 gb_yy_z
!        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty
!        real*8 gb_yy_tz
!        real*8 gb_yy_xx,gb_yy_xy
!        real*8 gb_yy_xz
!        real*8 gb_yy_yy
!        real*8 gb_yy_yz
!        real*8 gb_yy_zz
!
!        real*8 gb_yz_t, gb_yz_x, gb_yz_y
!        real*8 gb_yz_z
!        real*8 gb_yz_tt,gb_yz_tx,gb_yz_ty
!        real*8 gb_yz_tz
!        real*8 gb_yz_xx,gb_yz_xy
!        real*8 gb_yz_xz
!        real*8 gb_yz_yy
!        real*8 gb_yz_yz
!        real*8 gb_yz_zz
!
!        real*8 psi_t, psi_x, psi_y
!        real*8 psi_z
!        real*8 psi_tt,psi_tx,psi_ty
!        real*8 psi_tz
!        real*8 psi_xx,psi_xy
!        real*8 psi_xz
!        real*8 psi_yy
!        real*8 psi_yz
!        real*8 psi_zz
!
!!----------------------------------------------------------------------
!
!        dx=(x(2)-x(1))
!        dy=(y(2)-y(1))
!        dz=(z(2)-z(1))
!
!        ! set index bounds for main loop
!        is=2
!        ie=Nx-1
!        js=2
!        je=Ny-1
!        ks=2
!        ke=Nz-1
!
!        ! adjust index bounds to compensate for ghost_width
!        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
!        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
!        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
!        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
!        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
!        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)
!            if ( chr(i,j,k).ne.ex) then
!
!
!                ! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
!                ! identified with the bdy CFT stress-energy tensor one-point function
!                !these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor
!
!             q=1-rho0
!
!             quasiset_tt(i,j,k)=(2*(1+3*(z0**2+y0**2)
!     &                          -(-x0**2+z0**2+y0**2)/rho0**2)
!     &                           *(gb_xx_n(i,j,k)/q)
!     &                          +(4*x0*z0*(-3+2/rho0**2)
!     &                           -(3*rho0**2*y0)/sqrt(z0**2+y0**2))
!     &                           *(gb_xz_n(i,j,k)/q)
!     &                          -12*x0*y0*(gb_xy_n(i,j,k)/q)
!     &                          +(8*x0*y0
!     &                           *(gb_xy_n(i,j,k)/q))/rho0**2
!     &                          +(3*z0*rho0**2
!     &                           *(gb_xy_n(i,j,k)/q))/sqrt(z0**2+y0**2)
!     &                          +(4*z0**2*(psi_n(i,j,k)/q))/rho0**2
!     &                          +(6*x0**2*z0**2
!     &                           *(psi_n(i,j,k)/q))/(z0**2+y0**2)
!     &                          +(3*x0*z0*rho0**2*y0
!     &                           *(psi_n(i,j,k)/q))
!     &                           /(sqrt(-x0**2+rho0**2)
!     &                           *(z0**2+y0**2))
!     &                          +(8*z0*y0
!     &                           *(gb_yz_n(i,j,k)/q))/rho0**2
!     &                          -(3*x0*z0**2*rho0**2
!     &                           *(gb_yz_n(i,j,k)/q))
!     &                           /(sqrt(-x0**2+rho0**2)*(z0**2+y0**2))
!     &                          +(12*x0**2*z0*y0
!     &                           *(gb_yz_n(i,j,k)/q))/(z0**2+y0**2)
!     &                          +(3*x0*rho0**2*y0**2
!     &                           *(gb_yz_n(i,j,k)/q))
!     &                           /(sqrt(-x0**2+rho0**2)*(z0**2+y0**2))
!     &                          +(y0*(-((3*x0*z0*rho0**2)
!     &                           /sqrt(-x0**2+rho0**2))
!     &                           +6*x0**2*y0
!     &                           +(4*y0*(z0**2+y0**2))/rho0**2)
!     &                           *(gb_yy_n(i,j,k)/q))
!     &                           /(z0**2+y0**2))/(32*PI)
!
!             quasiset_tchi(i,j,k)  = (3*(-(z0**2+y0**2)
!     &                              *(gb_tx_n(i,j,k)/q)
!     &                              +x0*(z0*(gb_tz_n(i,j,k)/q)
!     &                               +y0*(gb_ty_n(i,j,k)/q))))
!     &                              /(16*sqrt(z0**2+y0**2))
!
!             quasiset_txi(i,j,k)   =(3.0d0/8.0d0)*(y0*(gb_tz_n(i,j,k)/q)
!     &                               -z0*(gb_ty_n(i,j,k)/q))
!
!             quasiset_chichi(i,j,k)=(1.0d0/64.0d0)*PI
!     &                              *(12*(gb_tt_n(i,j,k)/q)
!     &                              +(2*(-4*x0**2*(gb_xx_n(i,j,k)/q)
!     &                              +(-8*x0*z0+(3*rho0**4*y0)
!     &                               /sqrt(z0**2+y0**2))
!     &                               *(gb_xz_n(i,j,k)/q)
!     &                              -8*x0*y0*(gb_xy_n(i,j,k)/q)
!     &                              -(3*z0*rho0**4*(gb_xy_n(i,j,k)/q))
!     &                               /sqrt(z0**2+y0**2)
!     &                              -4*z0**2*(psi_n(i,j,k)/q)
!     &                              -(3*x0*z0*rho0**4*y0
!     &                               *(psi_n(i,j,k)/q))
!     &                               /(sqrt(-x0**2+rho0**2)
!     &                                *(z0**2+y0**2))
!     &                              -8*z0*y0*(gb_yz_n(i,j,k)/q)
!     &                              +(3*x0*z0**2*rho0**4
!     &                               *(gb_yz_n(i,j,k)/q))
!     &                               /(sqrt(-x0**2+rho0**2)
!     &                                *(z0**2+y0**2))
!     &                              -(3*x0*rho0**4*y0**2
!     &                               *(gb_yz_n(i,j,k)/q))
!     &                               /(sqrt(-x0**2+rho0**2)
!     &                                *(z0**2+y0**2))
!     &                              -4*y0**2*(gb_yy_n(i,j,k)/q)
!     &                              +(3*x0*z0*rho0**4*y0
!     &                               *(gb_yy_n(i,j,k)/q))
!     &                               /(sqrt(-x0**2+rho0**2)
!     &                                *(z0**2+y0**2))))/rho0**2)
!
!
!             quasiset_chixi(i,j,k) =(3*PI*(-y0*(z0**2+y0**2)
!     &                               *(gb_xz_n(i,j,k)/q)
!     &                              +z0*(z0**2+y0**2)
!     &                               *(gb_xy_n(i,j,k)/q)
!     &                              +x0*(z0*y0*(psi_n(i,j,k)/q)
!     &                              +(-z0**2+y0**2)*(gb_yz_n(i,j,k)/q)
!     &                              -z0*y0*(gb_yy_n(i,j,k)/q))))
!     &                               /(8*sqrt(z0**2+y0**2))
!
!             quasiset_xixi(i,j,k)  =-((PI*(z0**2+y0**2)
!     &                              *(-6*(gb_tt_n(i,j,k)/q)
!     &                              +(2*(x0**2+rho0**2
!     &                               +z0**2*(-1+3*rho0**2)
!     &                               +(-1+3*rho0**2)*y0**2)
!     &                               *(gb_xx_n(i,j,k)/q))/rho0**2
!     &                              +(4*x0*z0*(2-3*rho0**2)
!     &                               *(gb_xz_n(i,j,k)/q))/rho0**2
!     &                              +(4*x0*(2-3*rho0**2)*y0
!     &                               *(gb_xy_n(i,j,k)/q))/rho0**2
!     &                              +((2+3*rho0**2)*(z0**2
!     &                               *(psi_n(i,j,k)/q)
!     &                              +y0*(2*z0*(gb_yz_n(i,j,k)/q)
!     &                               +y0*(gb_yy_n(i,j,k)/q))))
!     &                               /(z0**2+y0**2)
!     &                              +((-2+3*rho0**2)
!     &                              *(x0**2-z0**2-y0**2)
!     &                              *(z0**2*(psi_n(i,j,k)/q)
!     &                               +y0*(2*z0* (gb_yz_n(i,j,k)/q)
!     &                               +y0*(gb_yy_n(i,j,k)/q))))
!     &                               /(rho0**2*(z0**2+y0**2))))
!     &                              /(8*rho0**2))
!
!             quasiset_mass(i,j,k)  =gb_xx_np1(i,j,k)
!
!            else
!
!             quasiset_tt(i,j,k)=0
!             quasiset_tchi(i,j,k)=0
!             quasiset_txi(i,j,k)=0
!             quasiset_chichi(i,j,k)=0
!             quasiset_chixi(i,j,k)=0
!             quasiset_xixi(i,j,k)=0
!
!
!            end if
