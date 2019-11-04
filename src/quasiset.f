c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 at points near the boundary. 
c The tensor components are given in spherical polar coordinates.
c----------------------------------------------------------------------

        subroutine quasiset_ll(
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
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_massdensityll,
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
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

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

        real*8 dxcar_dxsph(4,4)
        real*8 hcar_n(4,4)
        real*8 hsph_n(4,4),gbsph_n(4,4)

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

!!!!ghost_width not needed as lons as we don't use derivatives
!        ! adjust index bounds to compensate for ghost_width
!        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
!        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
!        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
!        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
!        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
!        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        do i=is,ie
         do j=js,je
          do k=ks,ke
            if ( chr(i,j,k).ne.ex) then


                ! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
                ! identified with the bdy CFT stress-energy tensor one-point function
                !these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor
             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0

!calculate regularized metric components in spherical coordinates in terms of regularized metric components in Cartesian coordinates
! we use the following coordinate transformation (notice that the angles are rescaled w.r.t. the usual spherical coordinates): x=rho*cos(PI*chi),y=rho*sin(PI*chi)*cos(2*PI*xi),z=rho*sin(PI*chi)*sin(2*PI*xi)

!transformation matrix

        dtdt=1
        dtdrho=0
        dtdchi=0
        dtdxi=0
        dxdt=0
        dxdrho=x0/rho0
        dxdchi=-PI*sqrt(y0**2+z0**2)
        dxdxi=0
        dydt=0
        dydrho=y0/rho0
        dydchi=PI*x0*y0/(sqrt(y0**2+z0**2))
        dydxi=-2*PI*z0
        dzdt=0
        dzdrho=z0/rho0
        dzdchi=PI*x0*z0/(sqrt(y0**2+z0**2))
        dzdxi=2*PI*y0

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

       do a=1,3
          do b=a+1,4
            gbsph_n(b,a)=gbsph_n(a,b)
          end do
        end do

             if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
             quasiset_tt_ll(i,j,k)=(12*(gbsph_n(3,3)/q)
     &                        + 8*PI**2*(gbsph_n(2,2)/q)
     &                        +  (3*rho0**2*(gbsph_n(4,4)/q))
     &                          /(z0**2+y0**2)
     &                        )/(64*PI**3)


             quasiset_tchi_ll(i,j,k)  = (3*(gbsph_n(1,3)/q))/(16*PI)


             quasiset_txi_ll(i,j,k)   = (3*(gbsph_n(1,4)/q))/(16*PI)

             quasiset_chichi_ll(i,j,k)=(3.0d0/16.0d0)*PI
     &                                  *(gbsph_n(1,1)/q)
     &                                 -(1.0d0/8.0d0)*PI
     &                                  *(gbsph_n(2,2)/q)
     &                                 -(3*rho0**2*(gbsph_n(4,4)/q))
     &                                  /(64*PI*(z0**2+y0**2))


             quasiset_chixi_ll(i,j,k) =(3*(gbsph_n(3,4)/q))/(16*PI)


             quasiset_xixi_ll(i,j,k)  =( (z0**2+y0**2)*(-3
     &                                  *(gbsph_n(3,3)/q)
     &                                 +PI**2*(3*(gbsph_n(1,1)/q)
     &                                 -2*(gbsph_n(2,2)/q)))
     &                                 )/(4*PI*rho0**2)

             quasiset_massdensityll(i,j,k)=(sqrt(y0**2+z0**2)
     &                                   *(12*(gbsph_n(3,3)/q)
     &                                   +8*PI**2*(gbsph_n(2,2)/q)
     &                                   +(3*rho0**2*(gbsph_n(4,4)/q))
     &                                   /(y0**2+z0**2)))
     &                                   /(32*PI*rho0)
             else
             quasiset_tt_ll(i,j,k)=0

             quasiset_tchi_ll(i,j,k)  =0


             quasiset_txi_ll(i,j,k)   =0

             quasiset_chichi_ll(i,j,k)=0

             quasiset_chixi_ll(i,j,k) =0

             quasiset_xixi_ll(i,j,k)  =0

             quasiset_massdensityll(i,j,k)=0
             end if


            else

             quasiset_tt_ll(i,j,k)=0
             quasiset_tchi_ll(i,j,k)=0
             quasiset_txi_ll(i,j,k)=0
             quasiset_chichi_ll(i,j,k)=0
             quasiset_chixi_ll(i,j,k)=0
             quasiset_xixi_ll(i,j,k)=0
             quasiset_massdensityll(i,j,k)=0


            end if

          end do
         end do
        end do


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
     &                  quasiset_tt,quasiset_tchi,quasiset_txi,
     &                  quasiset_chichi,quasiset_chixi,
     &                  quasiset_xixi,
     &                  quasiset_massdensity,AdS_mass,
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

        real*8 quasiset_tt_p1,quasiset_tchi_p1
        real*8 quasiset_txi_p1,quasiset_chichi_p1
        real*8 quasiset_chixi_p1,quasiset_xixi_p1
        real*8 quasiset_massdensity_p1

        real*8 quasiset_tt_p2,quasiset_tchi_p2
        real*8 quasiset_txi_p2,quasiset_chichi_p2
        real*8 quasiset_chixi_p2,quasiset_xixi_p2
        real*8 quasiset_massdensity_p2

        real*8 quasiset_tt(numbdypoints),quasiset_tchi(numbdypoints)
        real*8 quasiset_txi(numbdypoints),quasiset_chichi(numbdypoints)
        real*8 quasiset_chixi(numbdypoints),quasiset_xixi(numbdypoints)
        real*8 quasiset_massdensity(numbdypoints)

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

        real*8 xp1,yp1,zp1,xp2,yp2,zp2
        real*8 xex,yex,zex
        real*8 maxxyzp1

        real*8 extrapalongx,extrapalongy,extrapalongz

        real*8 dp1p2

        real*8 AdS_mass

!!!!!!!!!!!!!!!!

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        call quasiset_ll(
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
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_massdensityll,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)


        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

!        ! adjust index bounds to compensate for ghost_width
!        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
!        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
!        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
!        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
!        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
!        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)
       
!!!!!!!!!DEBUG!!!!!!!!
!        do lind=1,numbdypoints
!         quasiset_tt(lind)=3
!        end do
!!!!!!!!!!!!!!!!!!
 
        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           quasiset_tt_p1=quasiset_tt_ll(i,j,k)
           quasiset_tchi_p1=quasiset_tchi_ll(i,j,k)
           quasiset_txi_p1=quasiset_txi_ll(i,j,k)
           quasiset_chichi_p1=quasiset_chichi_ll(i,j,k)
           quasiset_chixi_p1=quasiset_chixi_ll(i,j,k)
           quasiset_xixi_p1=quasiset_xixi_ll(i,j,k)
           quasiset_massdensity_p1=quasiset_massdensityll(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1
              if (maxxyzp1.eq.abs(xp1)) then
               if (xp1.gt.0) then
 
                  xextrap(lind)=sqrt(1-yp1**2-zp1**2)
                  yextrap(lind)=yp1
                  zextrap(lind)=zp1
                  xp2=x(i-1)
                  quasiset_tt_p2=quasiset_tt_ll(i-1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i-1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i-1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i-1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i-1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i-1,j,k)
                  quasiset_massdensity_p2=
     &                           quasiset_massdensityll(i-1,j,k)
                  xex=xextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongx(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,xp1,xp2,xex)
 
 
              else
                  xextrap(lind)=-sqrt(1-yp1**2-zp1**2)
                  yextrap(lind)=yp1
                  zextrap(lind)=zp1
                  xp2=x(i+1)
                  quasiset_tt_p2=quasiset_tt_ll(i+1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i+1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i+1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i+1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i+1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i+1,j,k)
                  quasiset_massdensity_p2=
     &                                 quasiset_massdensityll(i+1,j,k)
                  xex=xextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongx(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,xp1,xp2,xex)
 
 
              end if
             else if (maxxyzp1.eq.abs(yp1)) then
              if (yp1.gt.0) then
                  yextrap(lind)=sqrt(1-xp1**2-zp1**2)
                  xextrap(lind)=xp1
                  zextrap(lind)=zp1
                  yp2=y(j-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j-1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j-1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j-1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j-1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j-1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j-1,k)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j-1,k)
                  yex=yextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongy(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,yp1,yp2,yex)   
 
 
              else
                  yextrap(lind)=-sqrt(1-xp1**2-zp1**2)
                  xextrap(lind)=xp1
                  zextrap(lind)=zp1
                  yp2=y(j+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j+1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j+1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j+1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j+1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j+1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j+1,k)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j+1,k)
                  yex=yextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongy(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,yp1,yp2,yex)
 
 
              end if             
             else
                 if (zp1.gt.0) then
                  zextrap(lind)=sqrt(1-yp1**2-xp1**2)
                  yextrap(lind)=yp1
                  xextrap(lind)=xp1
                  zp2=z(k-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k-1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k-1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k-1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k-1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k-1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k-1)
                quasiset_massdensity_p2=quasiset_massdensityll(i,j,k-1)
                  zex=zextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongz(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,zp1,zp2,zex)
 
              else
                  zextrap(lind)=-sqrt(1-yp1**2-xp1**2)
                  yextrap(lind)=yp1
                  xextrap(lind)=xp1
                  zp2=z(k+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k+1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k+1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k+1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k+1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k+1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k+1)
                 quasiset_massdensity_p2=quasiset_massdensityll(i,j,k+1)
                  zex=zextrap(lind)
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
          quasiset_massdensity(lind)=
     &                         extrapalongz(quasiset_massdensity_p1
     &                         ,quasiset_massdensity_p2,zp1,zp2,zex)
 
               end if
              end if
            end if
          end do
         end do
        end do

!calculate the value of angular coordinates chi and xi at the boundary points where we extrapolate the quasi-local stress energy tensor
        rhoextrap=1.0d0
        do e=1,numbdypoints
         chiextrap(e)=(1/PI)*acos(xextrap(e)/rhoextrap)
         if (zextrap(e).lt.0) then
             xiextrap(e)=(1/(2*PI))*(atan2(zextrap(e),yextrap(e))+2*PI)
         else
             xiextrap(e)=(1/(2*PI))*atan2(zextrap(e),yextrap(e))
         end if
        end do

!calculate the number of different values taken by chi and xi at the boundary points
        call bdy_N(
     &                  numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi)
!        write(*,*) "bdy_Nchi,bdy_Nxi=",bdy_Nchi,bdy_Nxi

!compute the double integral in chi and xi of the massdensity
        AdS_mass=0.0d0
        call doubleintegralonsphere(quasiset_massdensity,
     &                  xextrap,yextrap,zextrap,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi,
     &                  AdS_mass)

!       write(*,*) "AdS_mass=",AdS_mass


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

c-----------------------------------------------------------------------------------
c The following routine performs the double integral on the S^2 boundary. We approximate the integrals via the trapezoidal rule
c-----------------------------------------------------------------------------------

        subroutine doubleintegralonsphere(f,
     &                  xextrap,yextrap,zextrap,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi,
     &                  integral)

        implicit none
        integer numbdypoints
        real*8 f(numbdypoints)
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

        do i=1,bdy_Nchi
            chibdy(i)=0
        end do

        do j=1,bdy_Nxi
            xibdy(j)=0
        end do

        rhobdy=1.0d0

!obtain the arrays chibdy(bdy_Nchi) and xibdy(bdy_Nxi) containing the values of chi and xi at the boundary points, in increasing order
        call chibdy_xibdy(xextrap,yextrap,zextrap,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi,
     &                  chibdy,xibdy)


        additions=0
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

!          additions=additions+1

              integral=integral+
     &              (chibdy(i+1)-chibdy(i))/2 * (xibdy(j+1)-xibdy(j))/2
     &              *(f(lind_chipxip)+f(lind_chipxipp1)
     &              +f(lind_chipp1xip)+f(lind_chipp1xipp1))

         end do
        end do
 
!        write(*,*) "additions=",additions


        return
        end
!----------------------------------------------------------------------

c---------------------------------------------------------------------------------------------------------------------------------------------------------
c The following routine calculate the number of different values taken by chi and xi at the boundary points where the stress-energy tensor is extrapolated
c---------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine bdy_N(
     &                  numbdypoints,
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
     &                  xextrap,yextrap,zextrap,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi,
     &                  chibdy,xibdy)

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
!            write(*,*) "i,chibdy(i)=",i,chibdy(i)
        end do

        do j=1,bdy_Nxi 
            xiextrap_min=minval(xiextrap,mask=xiextrap.gt.xiextrap_min)
            xibdy(j)=xiextrap_min
!            write(*,*) "j,xibdy(j)=",j,xibdy(j)
        end do

        return
        end
!----------------------------------------------------------------------
