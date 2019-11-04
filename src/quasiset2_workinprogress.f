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
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if


!calculate regularized metric components in spherical coordinates in terms of regularized metric components in Cartesian coordinates
! we use the following coordinate transformation (notice that the angles are rescaled w.r.t. the usual spherical coordinates): x=rho*cos(PI*chi),y=rho*sin(PI*chi)*cos(2*PI*xi),z=rho*sin(PI*chi)*sin(2*PI*xi)

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
        dzdchi=PI*rho*cos(PI*chi0)*sin(2*PI*xi0)
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
     &                        +  (3*(gbsph_n(4,4)/q)/(sin(PI*chi0))**2)
     &                        )/(64*PI**3)


             quasiset_tchi_ll(i,j,k)  = (3*(gbsph_n(1,3)/q))/(16*PI)


             quasiset_txi_ll(i,j,k)   = (3*(gbsph_n(1,4)/q))/(16*PI)

             quasiset_chichi_ll(i,j,k)=(3.0d0/16.0d0)*PI
     &                                  *(gbsph_n(1,1)/q)
     &                                 -(1.0d0/8.0d0)*PI
     &                                  *(gbsph_n(2,2)/q)
     &                          -(3*(gbsph_n(4,4)/q)/(sin(PI*chi0))**2)
     &                                  /(64*PI)


             quasiset_chixi_ll(i,j,k) =(3*(gbsph_n(3,4)/q))/(16*PI)


             quasiset_xixi_ll(i,j,k)  =((sin(PI*chi0))**2*(-3
     &                                  *(gbsph_n(3,3)/q)
     &                                 +PI**2*(3*(gbsph_n(1,1)/q)
     &                                 -2*(gbsph_n(2,2)/q)))
     &                                 )/(4*PI)

             quasiset_massdensityll(i,j,k)=((sin(PI*chi0))**2
     &                                   *(12*(gbsph_n(3,3)/q)
     &                                   +8*PI**2*(gbsph_n(2,2)/q)
     &                                   +(3*(gbsph_n(4,4)/q))
     &                                    /(sin(PI*chi0))**2 ))
     &                                   /(32*PI)

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
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if


           quasiset_tt_p1=quasiset_tt_ll(i,j,k)
           quasiset_tchi_p1=quasiset_tchi_ll(i,j,k)
           quasiset_txi_p1=quasiset_txi_ll(i,j,k)
           quasiset_chichi_p1=quasiset_chichi_ll(i,j,k)
           quasiset_chixi_p1=quasiset_chixi_ll(i,j,k)
           quasiset_xixi_p1=quasiset_xixi_ll(i,j,k)
           quasiset_massdensity_p1=quasiset_massdensityll(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))
           maxxyp1=max(abs(xp1),abs(yp1))
           maxxzp1=max(abs(xp1),abs(zp1))
           maxyzp1=max(abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1
              if (maxxyzp1.eq.abs(xp1)) then
               if (xp1.gt.0) then



