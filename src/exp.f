c----------------------------------------------------------------------
c in cartesian coordinates t,x,y,z for x,y,z in [-1,1]
c
c support routines for apph.c
c
c (fixed using correct g=gads+(1-x^2)gb factorization)
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c returns true if a given point is not excised, and can be updated if
c in the *ghost region* of the grid
c----------------------------------------------------------------------
        logical function can_calc_ex(chr,i,j,k,Nx,Ny,Nz,ex)
        implicit none
        integer i,j,k,Nx,Ny,Nz
        real*8 chr(Nx,Ny,Nz),ex

        integer i1,j1,k1

        ! initialize fixed-size variables
        data i1,j1,k1/0,0,0/

        !--------------------------------------------------------------

        can_calc_ex=.true.

        if (chr(i,j,k).ne.ex.and.((i.gt.2.and.i.lt.Nx-1).and.
     &                            (j.gt.2.and.j.lt.Ny-1).and.
     &                            (k.gt.2.and.k.lt.Nz-1))
     &                         ) return

        do i1=max(1,min(Nx-2,i-1)),min(Nx,max(3,i+1))
           do j1=max(1,min(Ny-2,j-1)),min(Ny,max(3,j+1))
             do k1=max(1,min(Nz-2,k-1)),min(Nz,max(3,k+1))
              if (chr(i1,j1,k1).eq.ex) then
                 can_calc_ex=.false.
                 return
              end if
             end do
           end do
        end do

        return 
        end

c----------------------------------------------------------------------
c the following computes the outgoing null expansion (theta) 
c
c f is the level-surface function, and [is0..ie0,js0..je0,ks0..ke0] specifies
c the region overwhich theta should be computed
c
c is_ex is set to 1 if couldn't compute expansion an some point due
c to closeness of an excised region.
c----------------------------------------------------------------------
        subroutine calc_exp(theta,f,is0,ie0,js0,je0,ks0,ke0,is_ex,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    L,x,y,z,dt,chr,ex,do_ex,
     &                    Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,is0,ie0,js0,je0,ks0,ke0,do_ex
        integer is_ex
        real*8 theta(Nx,Ny,Nz),f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

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

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d,e

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 x0,y0,z0
        real*8 dx,dy,dz

        real*8 zeros(Nx,Ny,Nz)

        real*8 f_x,f_y,f_z,f_xx,f_xy,f_xz,f_yy,f_yz,f_zz
        real*8 tmp1,tmp2,tmp3,tmp4,tmp5

        real*8 n_l(4),s_l(4)
        real*8 n_u(4),s_u(4)
        real*8 n_l_x(4,4),n_u_x(4,4),s_l_x(4,4)
        real*8 f0_x(4) 
        real*8 f0_xx(4,4)
        real*8 gam_uu(4,4),sig_uu(4,4)
        real*8 gam_uu_x(4,4,4)
        real*8 normsusq

        real*8 g0gamfx(4)
        real*8 nufx,nuxfx(4),gamxfxfx(4)

        real*8 theta_ads

        integer do_ex_n
        logical ltrace,any_ex,can_calc_ex
        parameter (ltrace=.false.,do_ex_n=0)

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
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data i,j,k/0,0,0/
        data is,ie,js,je,ks,ke/0,0,0,0,0,0/

        data x0,y0,z0/0.0,0.0,0.0/
        data dx,dy,dz/0.0,0.0,0.0/

        data f_x,f_y,f_z/0.0,0.0,0.0/
        data f_xx,f_xy,f_xz/0.0,0.0,0.0/
        data f_yy,f_yz,f_zz/0.0,0.0,0.0/
        data tmp1,tmp2,tmp3,tmp4,tmp5/0.0,0.0,0.0,0.0,0.0/

        data n_l,s_l/4*0.0,4*0.0/
        data n_u,s_u/4*0.0,4*0.0/
        data n_l_x,n_u_x,s_l_x/16*0.0,16*0.0,16*0.0/
        data f0_x/4*0.0/
        data f0_xx/16*0.0/
        data gam_uu,sig_uu/16*0.0,16*0.0/
        data gam_uu_x/64*0.0/
        data normsusq/0.0/

        data g0gamfx/4*0.0/
        data nufx,nuxfx,gamxfxfx/0.0,4*0.0,4*0.0/

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data gads_ll,gads_uu/16*0.0,16*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/

        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/

        data g0_ll_xx/256*0.0/
        data gads_ll_xx/256*0.0/
        data h0_ll_xx/256*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data riemann_ulll/256*0.0/

        data A_l,Hads_l/4*0.0,4*0.0/
        data A_l_x/16*0.0/

        data phi10_x/4*0.0/
        data phi10_xx/16*0.0/

!----------------------------------------------------------------------

        if (ltrace) write(*,*) 'calc_exp ... N=',Nx,Ny,Nz

        is_ex=0

        do i=1,Nx
          do j=1,Ny
           do k=1,Nz
            zeros(i,j,k)=0
           end do
          end do
        end do

        do i=is0,ie0
          do j=js0,je0
           do k=ks0,ke0 
           theta(i,j,k)=1.0d0
           end do
          end do
        end do

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        is=max(2,is0)
        ie=min(Nx-1,ie0)
        js=max(2,js0)
        je=min(Ny-1,je0)
        ks=max(2,ks0)
        ke=min(Nz-1,ke0)

        do i=is,ie
          do j=js,je
           do k=ks,ke 
            if (ltrace) write(*,*) 'i,j,k:',i,j,k

            any_ex=.false.

            if (.not.(do_ex.eq.0.or.
     &          can_calc_ex(chr,i,j,k,Nx,Ny,Nz,ex))) then
                write(*,*) ' can_calc_ex: i,j,k has excised neighbors,'
                write(*,*) '              so cannot be updated '
                write(*,*) ' i,j,k=',i,j,k
                write(*,*) ' x(i),y(j),z(k)=',x(i),y(j),z(k)
                write(*,*) ' chr(i,j,k)=',chr(i,j,k)
                write(*,*) ' Nx,Ny,Nz=',Nx,Ny,Nz
                write(*,*) ' dx,dy,dz=',dx,dy,dz
               any_ex=.true.
               is_ex=1
            else
               x0=x(i)
               y0=y(j)              
               z0=z(k)
 
               ! computes tensors at point i,j 
               call tensor_init(
     &                 gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                 gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                 gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                 gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                 gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                 gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                 gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                 gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                 gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                 psi_np1,psi_n,psi_nm1,
     &                 zeros,zeros,zeros,
     &                 zeros,zeros,zeros,
     &                 zeros,zeros,zeros,
     &                 zeros,zeros,zeros,
     &                 zeros,zeros,zeros,
     &                 g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                 gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                 h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                 A_l,A_l_x,Hads_l,
     &                 gamma_ull,gamma_ull_x,
     &                 riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                 einstein_ll,set_ll,
     &                 phi10_x,phi10_xx,
     &                 x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k)

            end if

            if (.not.any_ex) then

              ! computes needed derivatives of f=r-R(chi,phi)
              call df2_int(f,f,f,
     &             tmp1, f_x, f_y, f_z,
     &             tmp2,tmp3,tmp4,tmp5,
     &             f_xx,f_xy,f_xz,
     &             f_yy,f_yz,
     &             f_zz,
     &             x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'f')

              ! define unit time-like vector n, normal to t=const
              ! surfaces
              n_l(1)=-1/sqrt(-g0_uu(1,1))
              do a=1,4
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)
              end do
              do b=1,4
                n_l_x(1,b)=-1/2.0d0/sqrt(-g0_uu(1,1))**3*g0_uu_x(1,1,b)
              end do
              do a=1,4
                do b=1,4
                  n_u_x(a,b)=n_l_x(1,b)*g0_uu(a,1)+
     &                       n_l_x(2,b)*g0_uu(a,2)+
     &                       n_l_x(3,b)*g0_uu(a,3)+
     &                       n_l_x(4,b)*g0_uu(a,4)+
     &                       n_l(1)*g0_uu_x(a,1,b)+
     &                       n_l(2)*g0_uu_x(a,2,b)+
     &                       n_l(3)*g0_uu_x(a,3,b)+
     &                       n_l(4)*g0_uu_x(a,4,b)
                end do
              end do

              ! define gradients of the flow field f=r-AH_R(chi,phi)
              f0_x(1)=0
              f0_x(2)=f_x
              f0_x(3)=f_y
              f0_x(4)=f_z
              f0_xx(1,1)=0
              f0_xx(1,2)=0
              f0_xx(1,3)=0
              f0_xx(1,4)=0
              f0_xx(2,2)=f_xx
              f0_xx(2,3)=f_xy
              f0_xx(2,4)=f_xz
              f0_xx(3,3)=f_yy
              f0_xx(3,4)=f_yz
              f0_xx(4,4)=f_zz

              do a=1,3
                do b=a+1,4
                  f0_xx(b,a)=f0_xx(a,b)
                end do
              end do

              ! define metric on codimension-1 surfaces 
              do a=1,4
                do b=1,4
                  gam_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)
                end do
              end do
              do a=1,4
                do b=1,4
                  do c=1,4
                    gam_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                             +n_u_x(a,c)*n_u(b)
     &                             +n_u(a)*n_u_x(b,c)
                  end do
                end do
              end do

              ! define unit space-like vector s, orthogonal to n and
              ! projected gradient of the flow field f
              do a=1,4
                s_u(a)=0.0d0
                do b=1,4
                  s_u(a)=s_u(a)+gam_uu(a,b)*f0_x(b)
                end do
              end do
              normsusq=0.0d0
              do a=1,4
                do b=1,4
                  normsusq=normsusq+gam_uu(a,b)*f0_x(a)*f0_x(b)
                end do
              end do
              do a=1,4
                s_u(a)=s_u(a)/sqrt(normsusq)
              end do
              do a=1,4
                s_l(a)=s_u(1)*g0_ll(a,1)+
     &                 s_u(2)*g0_ll(a,2)+
     &                 s_u(3)*g0_ll(a,3)+
     &                 s_u(4)*g0_ll(a,4)
              end do

              nufx=0
              do a=1,4
                nufx=nufx
     &              +n_u(a)*f0_x(a)
                nuxfx(a)=0
                gamxfxfx(a)=0
                do c=1,4
                  nuxfx(a)=nuxfx(a)
     &                    +n_u_x(c,a)*f0_x(c)
                  do d=1,4
                    gamxfxfx(a)=gamxfxfx(a)
     &                        +gam_uu_x(c,d,a)*f0_x(c)*f0_x(d)
     &                        +gam_uu(c,d)*f0_xx(c,a)*f0_x(d)
     &                        +gam_uu(c,d)*f0_x(c)*f0_xx(d,a)
                  end do
                end do
              end do
              do a=1,4
                do b=1,4
                  s_l_x(a,b)=
     &                     (f0_xx(a,b)+n_l_x(a,b)*nufx+n_l(a)*nuxfx(b))
     &                     /sqrt(normsusq)
     &                    -(f0_x(a)+n_l(a)*nufx)*gamxfxfx(b)
     &                     /2.0d0/sqrt(normsusq)**3
                end do
              end do

              ! define metric on codimension-2 surfaces 
              do a=1,4
                do b=1,4
                  sig_uu(a,b)=g0_uu(a,b)+n_u(a)*n_u(b)-s_u(a)*s_u(b)
                end do
              end do

              ! for theta: outward null expansion 
              theta(i,j,k)=0.0d0
              do c=1,4
                do d=1,4
                  theta(i,j,k)=theta(i,j,k)
     &                   +sig_uu(c,d)*(n_l_x(c,d)+s_l_x(c,d))
                  do e=1,4
                    theta(i,j,k)=theta(i,j,k)
     &                   -sig_uu(c,d)*gamma_ull(e,c,d)*(n_l(e)+s_l(e))
                  end do
                end do
              end do

              ! if theta is nan, then set it to the following instead
              if (.not.(theta(i,j,k).le.abs(theta(i,j,k)))) then
                 theta(i,j,k)=0.0d0
              end if
 
            end if
           end do
          end do
        end do

        !--------------------------------------------------------------
        ! fill in the borders if desired via first order extrapolation
        !--------------------------------------------------------------
 
        do j=js,je
         do k=ks,ke 
          if (is0.eq.1) theta(1,j,k)=theta(2,j,k)
          if (ie0.eq.Nx) theta(Nx,j,k)=theta(Nx-1,j,k)
         end do
        end do

        do i=is,ie
         do k=ks,ke 
          if (js0.eq.1) theta(i,1,k)=theta(i,2,k)
          if (je0.eq.Ny) theta(i,Ny,k)=theta(i,Ny-1,k)
         end do
        end do

        do i=is,ie
         do j=js,je 
          if (ks0.eq.1) theta(i,j,1)=theta(i,j,2)
          if (ke0.eq.Nz) theta(i,j,Nz)=theta(i,j,Nz-1)
         end do
        end do

        return
        end

c-----------------------------------------------------------------------
c given AH_R,AH_xc, computes the corresponding coordinate center
c and principle axis radii for excision ex_r0,ex_xc0
c-----------------------------------------------------------------------
        subroutine fill_ex_params(AH_R,AH_xc,ex_r0,ex_xc0,
     &                            AH_Nchi,AH_Nphi,dx,dy,dz,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3),ex_r0(3),ex_xc0(3)
        
        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer i,j,k
        real*8 x,y,z,AH_chi,AH_phi,dahchi,dahphi
        real*8 dx,dy,dz
        real*8 xmin,xmax,ymin,ymax,zmin,zmax

        ! initialize fixed-size variables
        data i,j,k/0,0,0/

        data x,y,z,AH_chi,AH_phi/0.0,0.0,0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/
        data xmin,xmax,ymin,ymax,zmin,zmax/0.0,0.0,0.0,0.0,0.0,0.0/

        !--------------------------------------------------------------

        xmin=1
        xmax=0
        ymin=1
        ymax=0
        zmin=1
        zmax=0

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        do i=1,AH_Nchi
          do j=1,AH_Nphi

            AH_chi=(i-1)*dahchi
            AH_phi=(j-1)*dahphi

            ! AH_R,AH_chi: polar coordinates of point on AH, wrt center of AH
            ! AH_xc(1),AH_xc(2),AH_xc(3): cartesian coordinates of center of AH, wrt origin
            ! x,y,z cartesian coordinates of point on AH, wrt origin
!            x=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!            y=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!            z=AH_R(i,j)*cos(AH_chi)+AH_xc(3)

            x=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
            y=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
            z=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

            xmin=min(x,xmin)
            ymin=min(y,ymin)
            zmin=min(z,zmin)

            xmax=max(x,xmax)
            ymax=max(y,ymax)
            zmax=max(z,ymax)

          end do
        end do

        !sets excision coordinate center as AH coordinate center
        ex_xc0(1)=AH_xc(1)
        ex_xc0(2)=AH_xc(2)
        ex_xc0(3)=AH_xc(3)

        !sets excision principal axis radii as ellipse semi-axes
        ex_r0(1)=(xmax-xmin)/2
        ex_r0(2)=(ymax-ymin)/2
        ex_r0(3)=(zmax-zmin)/2

        return
        end

c-----------------------------------------------------------------------
c the follow checks whether the coordinate ellipse a, 
c given by [r_a,xc_a], is entirely contained in b, given by [r_b,xc_b]
c If so, a_in_b and a_int_b is set to 1 (else 0)
c If a intersects b (but does *not* entirely contain b), 
c then a_int_b is set to 1 (else 0)
c
c NOTE: only the 2*dim "corner" points along the principle axis
c are checked, and so there are situations where a_int_b might
c incorrectly be cleared. However, for the kinds of AH's we're
c dealing with, this should not happen
c-----------------------------------------------------------------------
        subroutine is_inside(a_in_b,a_int_b,r_a,xc_a,r_b,xc_b,dim)
        implicit none
        integer a_in_b,dim,a_int_b,side
        real*8 r_a(dim),xc_a(dim)
        real*8 r_b(dim),xc_b(dim)

        real*8 d2,d0
        logical ltrace

        ! initialize fixed-size variables
        data d2,d0/0.0,0.0/

        !--------------------------------------------------------------
        ! we only check points along the principle axis
        !--------------------------------------------------------------

        a_in_b=1
        a_int_b=0

        d0=(xc_a(2) - xc_b(2))**2/r_b(2)**2
        if (dim.gt.2) d0=d0+(xc_a(3) - xc_b(3))**2/r_b(3)**2
        d2=(xc_a(1)+r_a(1) - xc_b(1))**2/r_b(1)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1
        d2=(xc_a(1)-r_a(1) - xc_b(1))**2/r_b(1)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1

        d0=(xc_a(1) - xc_b(1))**2/r_b(1)**2
        if (dim.gt.2) d0=d0+(xc_a(3) - xc_b(3))**2/r_b(3)**2
        d2=(xc_a(2)+r_a(2) - xc_b(2))**2/r_b(2)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1
        d2=(xc_a(2)-r_a(2) - xc_b(2))**2/r_b(2)**2 + d0
        if (d2.gt.1) a_in_b=0
        if (d2.lt.1) a_int_b=1

        if (dim.gt.2) then
           d0=(xc_a(1) - xc_b(1))**2/r_b(1)**2
           d0=d0+(xc_a(2) - xc_b(2))**2/r_b(2)**2
           d2=(xc_a(3)+r_a(3) - xc_b(3))**2/r_b(3)**2 + d0
           if (d2.gt.1) a_in_b=0
           if (d2.lt.1) a_int_b=1
           d2=(xc_a(3)-r_a(3) - xc_b(3))**2/r_b(3)**2 + d0
           if (d2.gt.1) a_in_b=0
           if (d2.lt.1) a_int_b=1
        end if

        !--------------------------------------------------------------
        ! the above would have missed an intersection where a
        ! is large, and b is small and away from a's principle axis.
        ! to check for an intersection, all of b's points must
        ! be one side or another of a:
        !--------------------------------------------------------------

        side=0

        d0=(xc_b(2) - xc_a(2))**2/r_a(2)**2
        if (dim.gt.2) d0=d0+(xc_b(3) - xc_a(3))**2/r_a(3)**2
        d2=(xc_b(1)+r_b(1) - xc_a(1))**2/r_a(1)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1
        d2=(xc_b(1)-r_b(1) - xc_a(1))**2/r_a(1)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1

        d0=(xc_b(1) - xc_a(1))**2/r_a(1)**2
        if (dim.gt.2) d0=d0+(xc_b(3) - xc_a(3))**2/r_a(3)**2
        d2=(xc_b(2)+r_b(2) - xc_a(2))**2/r_a(2)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1
        d2=(xc_b(2)-r_b(2) - xc_a(2))**2/r_a(2)**2 + d0
        if (d2.gt.1) side=side+1
        if (d2.lt.1) side=side-1

        if (dim.gt.2) then
           d0=(xc_b(1) - xc_a(1))**2/r_a(1)**2
           d0=d0+(xc_b(2) - xc_a(2))**2/r_a(2)**2
           d2=(xc_b(3)+r_b(3) - xc_a(3))**2/r_a(3)**2 + d0
           if (d2.gt.1) side=side+1
           if (d2.lt.1) side=side-1
           d2=(xc_b(3)-r_b(3) - xc_a(3))**2/r_a(3)**2 + d0
           if (d2.gt.1) side=side+1
           if (d2.lt.1) side=side-1
        end if

        if (abs(side).ne.2*dim)  a_int_b=1

        return
        end

c-----------------------------------------------------------------------
c is a point on the hypersurface r=R(chi) (relative to the center
c AH_xc) interior to the cartesian bounding box?? 
c-----------------------------------------------------------------------
        subroutine ah_is_int(is_int,AH_R,AH_xc,i,j,bbox,dx,dy,dz,
     &                       AH_Nchi,AH_Nphi,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi,i,j,is_int
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3)
        real*8 dx,dy,dz,bbox(6)
        
        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 x,y,z,AH_chi,AH_phi,dahchi,dahphi

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables 
        data x,y,AH_chi,AH_phi/0.0,0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/

        !--------------------------------------------------------------

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        AH_chi=(i-1)*dahchi
        AH_phi=(j-1)*dahphi

        ! AH_R,AH_chi: polar coordinates of point on AH, wrt center of AH
        ! AH_xc(1),AH_xc(2),AH_xc(3): cartesian coordinates of center of AH, wrt origin
        ! x,y,z cartesian coordinates of point on AH, wrt origin
!        x=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!        y=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!        z=AH_R(i,j)*cos(AH_chi)+AH_xc(3)

        x=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
        y=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
        z=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

        if ((x-bbox(1)).gt.(2.5*dx).and.(bbox(2)-x).gt.(2.5*dx).and.
     &      (y-bbox(3)).gt.(2.5*dy).and.(bbox(4)-y).gt.(2.5*dy).and.
     &      (z-bbox(5)).gt.(2.5*dz).and.(bbox(6)-z).gt.(2.5*dz))
     &  then
           is_int=1
        else
           is_int=0
        end if

        if (ltrace) then
           write(*,*)
           write(*,*) 'is_int: ',is_int
           write(*,*) 'bbox(1),bbox(2)=',bbox(1),bbox(2)
           write(*,*) 'bbox(3),bbox(4)=',bbox(3),bbox(4)
           write(*,*) 'bbox(5),bbox(6)=',bbox(5),bbox(6)
           write(*,*) 'x,y,z=',x,y,z
           write(*,*) 'R,AH_chi,AH_phi',AH_R(i,j),AH_chi,AH_phi
           write(*,*) 'xc,yc,zc',AH_xc(1),AH_xc(2),AH_xc(3)
        end if

        return
        end

c-----------------------------------------------------------------------
c marks unowned points (-1) as owned by node=rank if expansion can be
c calculated via interior stencils (2 *cells* away from boundaries,
c as calculated above) on node rank ... fills in AH_lev(..)=L as well,
c note that physical boundaries are exceptions; see ah_is_int() 
c-----------------------------------------------------------------------
        subroutine ah_fill_own(AH_R,AH_xc,AH_own,AH_lev,bbox,dx,dy,dz,
     &                         rank,L,AH_Nchi,AH_Nphi,axisym)
        implicit none
        integer axisym
        integer rank,AH_Nchi,AH_Nphi,L
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3)
        integer AH_own(AH_Nchi,AH_Nphi),AH_lev(AH_Nchi,AH_Nphi)
        real*8 dx,dy,dz,bbox(6)
        
        integer i,j,is_int

        logical ltrace
        parameter (ltrace=.false.)

        ! initialize fixed-size variables
        data i,j,is_int/0,0,0/

        !--------------------------------------------------------------

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              if (AH_own(i,j).eq.-1) then
                 call ah_is_int(is_int,AH_R,AH_xc,i,j,bbox,dx,dy,dz,
     &                          AH_Nchi,AH_Nphi,axisym)
                 if (is_int.eq.1) then
                    AH_own(i,j)=rank
                    AH_lev(i,j)=L
                 end if
                 if (ltrace) then
                    write(*,*) 'ah_fill_own: bbox=',bbox
                    write(*,*) 'i,j=',i,j,' own=',AH_own(i,j)
                 end if
              end if
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c evaluates the hypersurface function f = r - R(chi,phi) in
c the range [is..ie,js..je,ks..ke], with
c
c x=AH_R(chi,phi) + AH_xc(1)
c y=chi/PI + AH_xc(2) 
c
c AH_chi goes from 0 ("x" axis) to PI ("x" axis)
c AH_phi goes from 0 ("z" axis) to PI ("z" axis) 
c
c AH_chi and AH_phi are measured relative to center of AH, which in
c cartesian coordinates is given by AH_xc
c-----------------------------------------------------------------------
        subroutine ah_fill_f(AH_R,AH_xc,f,is,ie,js,je,ks,ke,x,y,z,
     &                       AH_Nchi,AH_Nphi,Nx,Ny,Nz,axisym)
        implicit none
        integer axisym
        integer is,ie,js,je,ks,ke,AH_Nchi,AH_Nphi,Nx,Ny,Nz
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3)
        real*8 x(Nx),y(Ny),z(Nz),f(Nx,Ny,Nz)
        
        integer i,j,k,i0,j0,is0,ie0,js0,je0,ks0,ke0
        real*8 AH_chi,AH_phi,r,xb,yb,zb,dahchi,dahphi,ft,fp,rtp
        real*8 xb0,yb0,zb0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical ltrace
        parameter (ltrace=.false.)

        ! initialize fixed-size variables
        data i,j,k,i0,j0,is0,ie0,js0,je0,ks0,ke0/0,0,0,0,0,0,0,0,0,0,0/

        data AH_chi,AH_phi,r,xb,yb,zb/0.0,0.0,0.0,0.0,0.0,0.0/
        data dahchi,dahphi,ft,fp,rtp/0.0,0.0,0.0,0.0,0.0/
        data xb0,yb0,zb0/0.0,0.0,0.0/

        data dx,dy,dz/0.0,0.0,0.0/ 
        !--------------------------------------------------------------

        if (ltrace) then
           write(*,*) 'ah_fill_f:'
           write(*,*) 'is,ie,js,je,ks,ke:',is,ie,js,je,ks,ke
        end if

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        xb0=AH_xc(1)
        yb0=AH_xc(2)
        zb0=AH_xc(3)

        is0=is
        js0=js
        ks0=ks
        ie0=ie
        je0=je
        ke0=ke

        do i=is0,ie0
           xb=x(i)
           do j=js0,je0
              yb=y(j)
              do k=ks0,ke0 
                 zb=z(k)

                 !xb,yb,zb: cartesian coordinates of some point "not on, but near" AH, wrt origin
                 !AH_xc(1),AH_xc(2),AH_xc(3): cartesian coordinates of center of AH, wrt origin
                 !r,AH_chi,AH_phi polar coordinates of the point "not on, but near" AH, wrt center of AH
                 !AH_R,AH_chi,AH_phi: polar coordinates of point "projected" on AH, wrt center of AH)
                 r=sqrt( (xb-AH_xc(1))**2
     &                  +(yb-AH_xc(2))**2 
     &                  +(zb-AH_xc(3))**2 )
!                 AH_chi=acos( zb-AH_xc(3)/r )
!                 AH_phi=atan2( yb-AH_xc(2),xb-AH_xc(1) )

                 AH_chi=acos( xb-AH_xc(1)/r )
                 AH_phi=atan2( zb-AH_xc(3),yb-AH_xc(2) )

                 if (AH_phi.lt.0) AH_phi=AH_phi+2*PI

                 !-----------------------------------------------------
                 ! bi-linear interpolate in AH_chi and AH_phi
                 !-----------------------------------------------------
                 i0=min(int(AH_chi/dahchi+1),AH_Nchi-1)
                 j0=min(int(AH_phi/dahphi+1),AH_Nphi-1)

                 if (r.eq.0) then
!                   write(*,*) 'WARNING!!! i0 being set to 1'
!                   write(*,*) 'center of AH is sampled by ah_fill_f()'
!                   write(*,*) '(AH points too close to center of AH)'
!                   write(*,*) '(... increase AH_r0 or resolution)'
!                    write(*,*) ' i0,j0=',i0,j0
!                    write(*,*) 'xb,yb=',xb,yb
!                    write(*,*) 'is,ie=',is,ie
!                    write(*,*) 'js,je=',js,je
!                    write(*,*) 'ks,ke=',ks,ke
!                    write(*,*) 'i,j=',i,j
!                    write(*,*) 'AH_R(i0,j0)',AH_R(i0,j0)
!                    write(*,*) 'AH_chi,AH_phi=',AH_chi,AH_phi
!                    write(*,*) 'r=',r
!                    write(*,*) 'rtp=',rtp
!                    write(*,*) '---------------------------'
                   i0=1
                 end if

                 if (ltrace) write(*,*) 'i0,j0=',i0,j0, ' Nchi,Nphi=',
     &              AH_Nchi,AH_Nphi

                 ft=((AH_chi-(i0-1)*dahchi))/dahchi
                 fp=((AH_phi-(j0-1)*dahphi))/dahphi

                 rtp=AH_R(i0,j0)*(1-ft)*(1-fp)+
     &               AH_R(i0+1,j0)*(ft)*(1-fp)+
     &               AH_R(i0,j0+1)*(1-ft)*(fp)+
     &               AH_R(i0+1,j0+1)*(ft)*(fp)

                 f(i,j,k)=r-rtp

                 if (ltrace) then
                    write(*,*) 'xb,yb,zb=',xb,yb,zb
                    write(*,*) 'AH_R(i0,j0)',AH_R(i0,j0)
                    write(*,*) 'AH_chi,AH_phi=',AH_chi,AH_phi
                    write(*,*) 'r=',r
                    write(*,*) 'rtp=',rtp
                    write(*,*) '---------------------------'
                 end if
             end do
           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c the following calculates the expansion at the single point (i0,j0)
c on the hypersurface f=r-AH_R(chi,phi) (with offset AH_xc)
c
c also returns the area element da at (i0,j0);
c and if (i0,j0) corresponds to chi=Pi/2, the equatorial circumference
c element d_ceq, and if (i0,j0) corresponds to phi=0, the polar
c circumference element d_cp 
c
c NOTE: if Nchi is even, will miss equatorial plane by dahchi
c
c area element is set to -1000 if point was too close to excision boundary
c to compute.
c
c NOTE: this routine modifies f and chi, but only in the vicinity of (i0,j0)
c
c if use_AH_new!=0 (in AH_new.inc), then computes chi using new tensor routines
c-----------------------------------------------------------------------
        subroutine calc_exp0(AH_R,AH_xc,AH_theta,
     &                    i0,j0,AH_Nchi,
     &                    AH_Nphi,theta,f,da,d_ceq,d_cp,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    L,x,y,z,dt,chr,ex,do_ex,
     &                    Nx,Ny,Nz,axisym)
        implicit none
        integer axisym
        integer Nx,Ny,Nz,i0,j0,AH_Nchi,AH_Nphi,do_ex
        real*8 theta(Nx,Ny,Nz),f(Nx,Ny,Nz),AH_xc(3),da,d_ceq,d_cp
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_theta(AH_Nchi,AH_Nphi)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

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

        real*8 cosx(Nx),cosy(Ny),cosz(Nz)
        real*8 sinx(Nx),siny(Ny),sinz(Nz)

        integer i,j,k,i1,j1,k1,is,ie,js,je,ks,ke,is_ex

        real*8 x0,y0,z0
        real*8 PI
        parameter (PI=3.141592653589793d0)
        real*8 dx,dy,dz,dahchi,dahphi,AH_chi,AH_phi,fx,fy,fz
        real*8 r_t,r_p,st,sp,ct,cp,r
        real*8 dxb_dt,dyb_dt,dzb_dt
        real*8 dxb_dp,dyb_dp,dzb_dp
        real*8 dtt,dpp,dpt
        real*8 gb_xx0,gb_xy0,gb_yy0,psi0

        real*8 drh_dch,dx_dch,dy_dch
        real*8 dx_dch_rh,dy_dch_rh
        real*8 g_xx0,g_xy0,g_yy0,g_phph0
        real*8 g_chch0,g_const0

        real*8 rho0,f0

        real*8 g0_xx_ads0,g0_xy_ads0,g0_yy_ads0,g0_psi_ads0

        integer is_bad,fill_later

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables 
        data i,j,k,i1,j1,k1/0.0,0.0,0.0,0.0,0.0,0.0/
        data is,ie,js,je,ks,ke/0.0,0.0,0.0,0.0,0.0,0.0/
        data is_ex/0.0/

        data x0,y0,z0/0.0,0.0,0.0/
        data dx,dy,dz/0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/
        data AH_chi,AH_phi,fx,fy,fz/0.0,0.0,0.0,0.0,0.0/
        data r_t,r_p,st,sp,ct,cp,r/0.0,0.0,0.0,0.0,0.0,0.0,0.0/
        data dxb_dt,dyb_dt,dzb_dt/0.0,0.0,0.0/
        data dxb_dp,dyb_dp,dzb_dp/0.0,0.0,0.0/
        data dtt,dpp,dpt/0.0,0.0,0.0/
        data gb_xx0,gb_xy0,gb_yy0/0.0,0.0,0.0/
        data psi0/0.0/

        data drh_dch,dx_dch,dy_dch/0.0,0.0,0.0/
        data dx_dch_rh,dy_dch_rh/0.0,0.0/
        data g_xx0,g_xy0,g_yy0,g_phph0/0.0,0.0,0.0,0.0/
        data g_chch0,g_const0/0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        AH_chi=(i0-1)*dahchi
        AH_phi=(j0-1)*dahphi

        !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
        !AH_xc(1),AH_xc(2),AH_xc(3) cartesian coordinates of center of AH, wrt origin
        !x0,y0,z0 cartesian coordinates of point on AH, wrt origin
!        x0=AH_R(i0,j0)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!        y0=AH_R(i0,j0)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!        z0=AH_R(i0,j0)*cos(AH_chi)+AH_xc(3)

        x0=AH_R(i0,j0)*cos(AH_chi)+AH_xc(1)
        y0=AH_R(i0,j0)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
        z0=AH_R(i0,j0)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

        da=0
        d_ceq=0
        d_cp=0

        ! extract (i,j,z) cartesian grid point closest to x0,y0,z0
        i=(x0-x(1))/dx+1
        j=(y0-y(1))/dy+1
        k=(z0-z(1))/dz+1

        ! for bilinear interpolation of theta from 
        ! (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1),(i,j+1,k+1),(i+1,j+1,k+1)
        !(NOTE: since, x0,y0,z0 does not necessarily lie on a grid point these fx,fy,fz are *not* identically zero)
        fx=(x0-((i-1)*dx+x(1)))/dx
        fy=(y0-((j-1)*dy+y(1)))/dy
        fz=(z0-((k-1)*dz+z(1)))/dz

        if (ltrace) then
           write(*,*) 'calc_exp0: i0,j0=',i0,j0
           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
           write(*,*) ' i,j,k=',i,j,k
           write(*,*) ' x0,y0,z0=',x0,y0,z0
           write(*,*) ' x(1),y(1),z(1)=',x(1),y(1),z(1)
           write(*,*) ' dx,dy,dz=',dx,dy,dz
           write(*,*) '----------------------------'
        end if

        !------------------------------------------------------------
        ! we will use bilinear interpolation to calculate theta,
        ! hence we need f in a sufficient box around to caculate theta
        ! in a 2^2 box adjacent to (i,j)
        !------------------------------------------------------------

        is_bad=0
        fill_later=0
        if (i.lt.2.or.j.lt.2.or.k.lt.2.or.i.ge.Nx.or.j.ge.Ny.or.k.ge.Nz)
     &  then
          if (i0.eq.1.or.i0.eq.AH_Nchi.or.j0.eq.1.or.j0.eq.AH_Nphi)
     &    then !these are allowed ... reg_ah_r() will fill it in 
            fill_later=1
          else
            is_bad=1 
          end if
        end if

        if (AH_R(i0,j0).lt.(1.5d0*dx)) is_bad=1 !stops ahfinder when AH_R gets too small

        if (is_bad.eq.1) then
           write(*,*) '--------------------------------------------'
           write(*,*) 'WARNING!!! theta set to 100 at i0,j0=',i0,j0
           write(*,*) 'and regularity will *not* take care of it!'
           write(*,*) '(AH points too bunched up near i0 or j0 endpts)'
           write(*,*) '(... reduce AH_Nchi,AH_Nphi or increase res)'
!           write(*,*) ' i0,j0=',i0,j0
           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
!           write(*,*) ' xc(1),xc(2),xc(3)=',AH_xc(1),AH_xc(2),AH_xc(3)
           write(*,*) ' i,j,k=',i,j,k
           write(*,*) ' x0,y0,z0=',x0,y0,z0
!           write(*,*) ' x(1),y(1),z(1)=',x(1),y(1),z(1)
           write(*,*) ' dx,dy,dz=',dx,dy,dz
!           write(*,*) ' dahchi,dahphi=',dahchi,dahphi
!           write(*,*) '--------------------------------------------'
           AH_theta(i0,j0)=100
           return
        end if

        if (fill_later.eq.1) then
            AH_theta(i0,j0)=0
            return
        else

          is=max(1,i-2)
          ie=min(Nx,i+3)
          js=max(1,j-2)
          je=min(Ny,j+3)
          ks=max(1,k-2)
          ke=min(Nz,k+3)

          call ah_fill_f(AH_R,AH_xc,f,is,ie,js,je,ks,ke,x,y,z,
     &                AH_Nchi,AH_Nphi,Nx,Ny,Nz,axisym)

          is=i
          ie=i+1
          js=j
          je=j+1
          ks=k
          ke=k+1

          if (do_ex.ne.0) then
             do i1=is,ie
                do j1=js,je
                   do k1=ks,ke
                      if (chr(i1,j1,k1).eq.ex) then
                         write(*,*) ' calc_exp0: pt i1,j1,k1 is excised'
                         write(*,*) ' i0,j0=',i0,j0
                         write(*,*) ' AH_Nchi,AH_Nphi=',AH_Nchi,AH_Nphi
                         write(*,*) ' AH_chi,AH_phi=',AH_chi,AH_phi
                         write(*,*) ' x0,y0,z0=',x0,y0,z0
                         write(*,*) ' i,j,k=',i,j,k
                         write(*,*) ' xi,yj,zk=',x(i),y(j),z(k)
                         write(*,*) ' chr(i,j,k)=',chr(i,j,k)
                         write(*,*) ' i1,j1,k1=',i1,j1,k1
                         write(*,*) ' xi1,yj1,zk1,AH_R(i0,j0)='
     &                               ,x(i1),y(j1),z(k1),AH_R(i0,j0)
                         write(*,*) ' chr(i1,j1,k1)=',chr(i1,j1,k1)
                    write(*,*) ' rho1=',sqrt(x(i1)**2+y(j1)**2+z(k1)**2)
                         write(*,*) ' Nx,Ny,Nz=',Nx,Ny,Nz
                         da=-10000
                         AH_theta(i0,j0)=0
                         return
                      end if
                   end do
                end do
             end do
          end if

          call calc_exp(theta,f,is,ie,js,je,ks,ke,is_ex,
     &               gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &               gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &               gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &               gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &               gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &               gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &               gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &               gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &               gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &               psi_np1,psi_n,psi_nm1,
     &               L,x,y,z,dt,chr,ex,do_ex,Nx,Ny,Nz)

          ! interpolate theta from 
          ! (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1),(i,j+1,k+1),(i+1,j+1,k+1)
          AH_theta(i0,j0)=(1-fx)*(1-fy)*(1-fz)*theta(i,j,k)+
     &                    (  fx)*(1-fy)*(1-fz)*theta(i+1,j,k)+
     &                    (1-fx)*(  fy)*(1-fz)*theta(i,j+1,k)+
     &                    (1-fx)*(1-fy)*(  fz)*theta(i,j,k+1)+
     &                    (  fx)*(  fy)*(1-fz)*theta(i+1,j+1,k)+
     &                    (  fx)*(1-fy)*(  fz)*theta(i+1,j,k+1)+
     &                    (1-fx)*(  fy)*(  fz)*theta(i,j+1,k+1)+
     &                    (  fx)*(  fy)*(  fz)*theta(i+1,j+1,k+1)
        end if

        !--------------------------------------------------------------
        ! proper area element
        ! NOTE: use AdS5Ds to fill this out
        !--------------------------------------------------------------
!
!        rho0=sqrt(x0**2+y0**2)
!        f0=(1-rho0**2)**2+4*rho0**2/L**2
!
!        drh_dch=(AH_R(i0+1,j0)-AH_R(i0-1,j0))/2/dahchi
!
!        dx_dch=cos(AH_chi)*drh_dch-rho0*sin(AH_chi)  !when rho=rho(chi)
!        dy_dch=sin(AH_chi)*drh_dch+rho0*cos(AH_chi)  
!
!        dx_dch_rh=rho0*sin(AH_chi)  !when rho=const
!        dy_dch_rh=rho0*cos(AH_chi)
!
!        g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *4
!        g0_xy_ads0 =((1+rho0**2)**2/f0-1)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *x0*y0
!     &             *4
!        g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *4
!        g0_psi_ads0=(y0**2)
!     &             /(1-rho0**2)**2
!     &             *4
!
!        gb_xx0=((1-fx)*(1-fy)*gb_xx_n(i,j,k)+
!     &       (  fx)*(1-fy)*gb_xx_n(i+1,j,k)+
!     &       (1-fx)*(  fy)*gb_xx_n(i,j+1,k)+
!     &       (  fx)*(  fy)*gb_xx_n(i+1,j+1,k))
!
!        gb_xy0=((1-fx)*(1-fy)*gb_xy_n(i,j,k)+
!     &       (  fx)*(1-fy)*gb_xy_n(i+1,j,k)+
!     &       (1-fx)*(  fy)*gb_xy_n(i,j+1,k)+
!     &       (  fx)*(  fy)*gb_xy_n(i+1,j+1,k))
!
!        gb_yy0=((1-fx)*(1-fy)*gb_yy_n(i,j,k)+
!     &       (  fx)*(1-fy)*gb_yy_n(i+1,j,k)+
!     &       (1-fx)*(  fy)*gb_yy_n(i,j+1,k)+
!     &       (  fx)*(  fy)*gb_yy_n(i+1,j+1,k))
!
!        psi0=((1-fx)*(1-fy)*psi_n(i,j,k)+
!     &       (  fx)*(1-fy)*psi_n(i+1,j,k)+
!     &       (1-fx)*(  fy)*psi_n(i,j+1,k)+
!     &       (  fx)*(  fy)*psi_n(i+1,j+1,k))
!
!        g_xx0   =g0_xx_ads0+gb_xx0
!        g_xy0   =g0_xy_ads0+gb_xy0
!        g_yy0   =g0_yy_ads0+gb_yy0
!        g_phph0 =g0_psi_ads0+psi0*y0**2
!
!        g_chch0 =g_xx0*dx_dch_rh**2+2*g_xy0*dx_dch_rh*dy_dch_rh
!     &          +g_yy0*dy_dch_rh**2
!        g_const0=g_xx0*dx_dch**2+2*g_xy0*dx_dch*dy_dch+g_yy0*dy_dch**2
!
!        da=2*PI*sqrt(g_const0*g_phph0)*dahchi
!        if (is_ex.eq.1) then
!           da=-10000
!        end if
!
!        if (i0.eq.((AH_Nchi-1)/2)+1) then
!           d_ceq=2*PI*sqrt(g_phph0)
!        end if
!        d_cp=2*sqrt(g_chch0)*dahchi

        return
        end


c---------------------------------------------------------------------------------------------------
c the following calculates the metric components in Cartesian coordinates at the single point (i0,j0)
c on the hypersurface f=r-AH_R(chi,phi) (with offset AH_xc)
c---------------------------------------------------------------------------------------------------

        subroutine calc_ahmetric0(AH_R,AH_xc,
     &                    AH_g0_tt,
     &                    AH_g0_tx,AH_g0_ty,AH_g0_tz,
     &                    AH_g0_xx,AH_g0_xy,AH_g0_xz,
     &                    AH_g0_yy,AH_g0_yz,AH_g0_psi,
     &                    i0,j0,AH_Nchi,
     &                    AH_Nphi,
     &                    gb_tt_n,
     &                    gb_tx_n,
     &                    gb_ty_n,
     &                    gb_tz_n,
     &                    gb_xx_n,
     &                    gb_xy_n,
     &                    gb_xz_n,
     &                    gb_yy_n,
     &                    gb_yz_n,
     &                    psi_n,
     &                    L,x,y,z,dt,chr,ex,do_ex,
     &                    Nx,Ny,Nz,axisym)

        implicit none
        integer axisym
        integer Nx,Ny,Nz,i0,j0,AH_Nchi,AH_Nphi,do_ex
        real*8 theta(Nx,Ny,Nz),f(Nx,Ny,Nz),AH_xc(3),da,d_ceq,d_cp
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_theta(AH_Nchi,AH_Nphi)
        real*8 AH_g0_tt(AH_Nchi,AH_Nphi)
        real*8 AH_g0_tx(AH_Nchi,AH_Nphi),AH_g0_ty(AH_Nchi,AH_Nphi)
        real*8 AH_g0_tz(AH_Nchi,AH_Nphi)
        real*8 AH_g0_xx(AH_Nchi,AH_Nphi),AH_g0_xy(AH_Nchi,AH_Nphi)
        real*8 AH_g0_xz(AH_Nchi,AH_Nphi)
        real*8 AH_g0_yy(AH_Nchi,AH_Nphi),AH_g0_yz(AH_Nchi,AH_Nphi)
        real*8 AH_g0_psi(AH_Nchi,AH_Nphi)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

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

        real*8 cosx(Nx),cosy(Ny),cosz(Nz)
        real*8 sinx(Nx),siny(Ny),sinz(Nz)

        integer i,j,k,i1,j1,k1,is,ie,js,je,ks,ke,is_ex

        real*8 x0,y0,z0
        real*8 xi,yj,zk
        real*8 xip1,yjp1,zkp1
        real*8 rhoijk,rhoip1jk,rhoijp1k,rhoip1jp1k,rhoijkp1
        real*8 rhoip1jkp1,rhoijp1kp1,rhoip1jp1kp1
        real*8 PI
        parameter (PI=3.141592653589793d0)
        real*8 dx,dy,dz,dahchi,dahphi,AH_chi,AH_phi,fx,fy,fz
        real*8 r_t,r_p,st,sp,ct,cp,r
        real*8 dxb_dt,dyb_dt,dzb_dt
        real*8 dxb_dp,dyb_dp,dzb_dp
        real*8 dtt,dpp,dpt
        real*8 gb_xx0,gb_xy0,gb_yy0,psi0

        real*8 drh_dch,dx_dch,dy_dch
        real*8 dx_dch_rh,dy_dch_rh
        real*8 g_xx0,g_xy0,g_yy0,g_phph0
        real*8 g_chch0,g_const0

        real*8 rho0,f0

        real*8 g0_tt_adsijk
        real*8 g0_tt_adsip1jk
        real*8 g0_tt_adsijp1k
        real*8 g0_tt_adsip1jp1k
        real*8 g0_tt_adsijkp1
        real*8 g0_tt_adsip1jkp1
        real*8 g0_tt_adsijp1kp1
        real*8 g0_tt_adsip1jp1kp1

        real*8 g0_tx_adsijk
        real*8 g0_tx_adsip1jk
        real*8 g0_tx_adsijp1k
        real*8 g0_tx_adsip1jp1k
        real*8 g0_tx_adsijkp1
        real*8 g0_tx_adsip1jkp1
        real*8 g0_tx_adsijp1kp1
        real*8 g0_tx_adsip1jp1kp1

        real*8 g0_ty_adsijk
        real*8 g0_ty_adsip1jk
        real*8 g0_ty_adsijp1k
        real*8 g0_ty_adsip1jp1k
        real*8 g0_ty_adsijkp1
        real*8 g0_ty_adsip1jkp1
        real*8 g0_ty_adsijp1kp1
        real*8 g0_ty_adsip1jp1kp1

        real*8 g0_tz_adsijk
        real*8 g0_tz_adsip1jk
        real*8 g0_tz_adsijp1k
        real*8 g0_tz_adsip1jp1k
        real*8 g0_tz_adsijkp1
        real*8 g0_tz_adsip1jkp1
        real*8 g0_tz_adsijp1kp1
        real*8 g0_tz_adsip1jp1kp1

        real*8 g0_xx_adsijk
        real*8 g0_xx_adsip1jk
        real*8 g0_xx_adsijp1k
        real*8 g0_xx_adsip1jp1k
        real*8 g0_xx_adsijkp1
        real*8 g0_xx_adsip1jkp1
        real*8 g0_xx_adsijp1kp1
        real*8 g0_xx_adsip1jp1kp1

        real*8 g0_xy_adsijk
        real*8 g0_xy_adsip1jk
        real*8 g0_xy_adsijp1k
        real*8 g0_xy_adsip1jp1k
        real*8 g0_xy_adsijkp1
        real*8 g0_xy_adsip1jkp1
        real*8 g0_xy_adsijp1kp1
        real*8 g0_xy_adsip1jp1kp1

        real*8 g0_xz_adsijk
        real*8 g0_xz_adsip1jk
        real*8 g0_xz_adsijp1k
        real*8 g0_xz_adsip1jp1k
        real*8 g0_xz_adsijkp1
        real*8 g0_xz_adsip1jkp1
        real*8 g0_xz_adsijp1kp1
        real*8 g0_xz_adsip1jp1kp1

        real*8 g0_yy_adsijk
        real*8 g0_yy_adsip1jk
        real*8 g0_yy_adsijp1k
        real*8 g0_yy_adsip1jp1k
        real*8 g0_yy_adsijkp1
        real*8 g0_yy_adsip1jkp1
        real*8 g0_yy_adsijp1kp1
        real*8 g0_yy_adsip1jp1kp1

        real*8 g0_yz_adsijk
        real*8 g0_yz_adsip1jk
        real*8 g0_yz_adsijp1k
        real*8 g0_yz_adsip1jp1k
        real*8 g0_yz_adsijkp1
        real*8 g0_yz_adsip1jkp1
        real*8 g0_yz_adsijp1kp1
        real*8 g0_yz_adsip1jp1kp1

        real*8 g0_psi_adsijk
        real*8 g0_psi_adsip1jk
        real*8 g0_psi_adsijp1k
        real*8 g0_psi_adsip1jp1k
        real*8 g0_psi_adsijkp1
        real*8 g0_psi_adsip1jkp1
        real*8 g0_psi_adsijp1kp1
        real*8 g0_psi_adsip1jp1kp1

        real*8 g0_tt_ijk
        real*8 g0_tt_ip1jk
        real*8 g0_tt_ijp1k
        real*8 g0_tt_ip1jp1k
        real*8 g0_tt_ijkp1
        real*8 g0_tt_ip1jkp1
        real*8 g0_tt_ijp1kp1
        real*8 g0_tt_ip1jp1kp1

        real*8 g0_tx_ijk
        real*8 g0_tx_ip1jk
        real*8 g0_tx_ijp1k
        real*8 g0_tx_ip1jp1k
        real*8 g0_tx_ijkp1
        real*8 g0_tx_ip1jkp1
        real*8 g0_tx_ijp1kp1
        real*8 g0_tx_ip1jp1kp1

        real*8 g0_ty_ijk
        real*8 g0_ty_ip1jk
        real*8 g0_ty_ijp1k
        real*8 g0_ty_ip1jp1k
        real*8 g0_ty_ijkp1
        real*8 g0_ty_ip1jkp1
        real*8 g0_ty_ijp1kp1
        real*8 g0_ty_ip1jp1kp1

        real*8 g0_tz_ijk
        real*8 g0_tz_ip1jk
        real*8 g0_tz_ijp1k
        real*8 g0_tz_ip1jp1k
        real*8 g0_tz_ijkp1
        real*8 g0_tz_ip1jkp1
        real*8 g0_tz_ijp1kp1
        real*8 g0_tz_ip1jp1kp1

        real*8 g0_xx_ijk
        real*8 g0_xx_ip1jk
        real*8 g0_xx_ijp1k
        real*8 g0_xx_ip1jp1k
        real*8 g0_xx_ijkp1
        real*8 g0_xx_ip1jkp1
        real*8 g0_xx_ijp1kp1
        real*8 g0_xx_ip1jp1kp1

        real*8 g0_xy_ijk
        real*8 g0_xy_ip1jk
        real*8 g0_xy_ijp1k
        real*8 g0_xy_ip1jp1k
        real*8 g0_xy_ijkp1
        real*8 g0_xy_ip1jkp1
        real*8 g0_xy_ijp1kp1
        real*8 g0_xy_ip1jp1kp1

        real*8 g0_xz_ijk
        real*8 g0_xz_ip1jk
        real*8 g0_xz_ijp1k
        real*8 g0_xz_ip1jp1k
        real*8 g0_xz_ijkp1
        real*8 g0_xz_ip1jkp1
        real*8 g0_xz_ijp1kp1
        real*8 g0_xz_ip1jp1kp1

        real*8 g0_yy_ijk
        real*8 g0_yy_ip1jk
        real*8 g0_yy_ijp1k
        real*8 g0_yy_ip1jp1k
        real*8 g0_yy_ijkp1
        real*8 g0_yy_ip1jkp1
        real*8 g0_yy_ijp1kp1
        real*8 g0_yy_ip1jp1kp1

        real*8 g0_yz_ijk
        real*8 g0_yz_ip1jk
        real*8 g0_yz_ijp1k
        real*8 g0_yz_ip1jp1k
        real*8 g0_yz_ijkp1
        real*8 g0_yz_ip1jkp1
        real*8 g0_yz_ijp1kp1
        real*8 g0_yz_ip1jp1kp1

        real*8 g0_psi_ijk
        real*8 g0_psi_ip1jk
        real*8 g0_psi_ijp1k
        real*8 g0_psi_ip1jp1k
        real*8 g0_psi_ijkp1
        real*8 g0_psi_ip1jkp1
        real*8 g0_psi_ijp1kp1
        real*8 g0_psi_ip1jp1kp1


        real*8 g0_xx_ads0,g0_xy_ads0,g0_yy_ads0,g0_psi_ads0

        integer is_bad,fill_later

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables 
        data i,j,k,i1,j1,k1/0.0,0.0,0.0,0.0,0.0,0.0/
        data is,ie,js,je,ks,ke/0.0,0.0,0.0,0.0,0.0,0.0/
        data is_ex/0.0/

        data x0,y0,z0/0.0,0.0,0.0/
        data dx,dy,dz/0.0,0.0,0.0/
        data dahchi,dahphi/0.0,0.0/
        data AH_chi,AH_phi,fx,fy,fz/0.0,0.0,0.0,0.0,0.0/
        data r_t,r_p,st,sp,ct,cp,r/0.0,0.0,0.0,0.0,0.0,0.0,0.0/
        data dxb_dt,dyb_dt,dzb_dt/0.0,0.0,0.0/
        data dxb_dp,dyb_dp,dzb_dp/0.0,0.0,0.0/
        data dtt,dpp,dpt/0.0,0.0,0.0/
        data gb_xx0,gb_xy0,gb_yy0/0.0,0.0,0.0/
        data psi0/0.0/

        data drh_dch,dx_dch,dy_dch/0.0,0.0,0.0/
        data dx_dch_rh,dy_dch_rh/0.0,0.0/
        data g_xx0,g_xy0,g_yy0,g_phph0/0.0,0.0,0.0,0.0/
        data g_chch0,g_const0/0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        AH_chi=(i0-1)*dahchi
        AH_phi=(j0-1)*dahphi

        !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
        !AH_xc(1),AH_xc(2),AH_xc(3) cartesian coordinates of center of AH, wrt origin
        !x0,y0,z0 cartesian coordinates of point on AH, wrt origin
!        x0=AH_R(i0,j0)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!        y0=AH_R(i0,j0)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!        z0=AH_R(i0,j0)*cos(AH_chi)+AH_xc(3)

        x0=AH_R(i0,j0)*cos(AH_chi)+AH_xc(1)
        y0=AH_R(i0,j0)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
        z0=AH_R(i0,j0)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

        ! extract (i,j,z) cartesian grid point closest to x0,y0,z0
        i=(x0-x(1))/dx+1
        j=(y0-y(1))/dy+1
        k=(z0-z(1))/dz+1

        ! for bilinear interpolation of theta from 
        ! (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1),(i,j+1,k+1),(i+1,j+1,k+1)
        !(NOTE: since, x0,y0,z0 does not necessarily lie on a grid point these fx,fy,fz are *not* identically zero)
        fx=(x0-((i-1)*dx+x(1)))/dx
        fy=(y0-((j-1)*dy+y(1)))/dy
        fz=(z0-((k-1)*dz+z(1)))/dz


        if (ltrace) then
           write(*,*) 'calc_AHmetric0: i0,j0=',i0,j0
           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
           write(*,*) ' i,j,k=',i,j,k
           write(*,*) ' x0,y0,z0=',x0,y0,z0
           write(*,*) ' x(1),y(1),z(1)=',x(1),y(1),z(1)
           write(*,*) ' dx,dy,dz=',dx,dy,dz
           write(*,*) '----------------------------'
        end if

        !------------------------------------------------------------
        ! we will use bilinear interpolation to calculate theta,
        ! hence we need f in a sufficient box around to caculate theta
        ! in a 2^2 box adjacent to (i,j)
        !------------------------------------------------------------

        is_bad=0
        fill_later=0
        if (i.lt.2.or.j.lt.2.or.k.lt.2.or.i.ge.Nx.or.j.ge.Ny.or.k.ge.Nz)
     &  then
          if (i0.eq.1.or.i0.eq.AH_Nchi.or.j0.eq.1.or.j0.eq.AH_Nphi)
     &    then !these are allowed ... reg_ah_r() will fill it in 
            fill_later=1
          else
            is_bad=1
          end if
        end if

        if (AH_R(i0,j0).lt.(1.5d0*dx)) is_bad=1 !stops ahfinder when AH_R gets too small

        if (is_bad.eq.1) then
           write(*,*) '--------------------------------------------'
         write(*,*) 'WARNING!!! metric components set to 1E10 at i0,j0='
     &               ,i0,j0
           write(*,*) 'and regularity will *not* take care of it!'
           write(*,*) '(AH points too bunched up near i0 or j0 endpts)'
           write(*,*) '(... reduce AH_Nchi,AH_Nphi or increase res)'
!           write(*,*) ' i0,j0=',i0,j0
!           write(*,*) ' AH_R,AH_chi,AH_phi=',AH_R(i0,j0),AH_chi,AH_phi
!           write(*,*) ' xc(1),xc(2),xc(3)=',AH_xc(1),AH_xc(2),AH_xc(3)
!           write(*,*) ' i,j,k=',i,j,k
!           write(*,*) ' x0,y0,z0=',x0,y0,z0
!           write(*,*) ' x(1),y(1),z(1)=',x(1),y(1),z(1)
!           write(*,*) ' dx,dy,dz=',dx,dy,dz
!           write(*,*) ' dahchi,dahphi=',dahchi,dahphi
!           write(*,*) '--------------------------------------------'
           AH_g0_tt(i0,j0)=1E10
           AH_g0_tx(i0,j0)=1E10
           AH_g0_ty(i0,j0)=1E10
           AH_g0_tz(i0,j0)=1E10
           AH_g0_xx(i0,j0)=1E10
           AH_g0_xy(i0,j0)=1E10
           AH_g0_xy(i0,j0)=1E10
           AH_g0_yy(i0,j0)=1E10
           AH_g0_yz(i0,j0)=1E10
           AH_g0_psi(i0,j0)=1E10
           return
        end if

        if (fill_later.eq.1) then
           AH_g0_tt(i0,j0)=0
           AH_g0_tx(i0,j0)=0
           AH_g0_ty(i0,j0)=0
           AH_g0_tz(i0,j0)=0
           AH_g0_xx(i0,j0)=0
           AH_g0_xy(i0,j0)=0
           AH_g0_xy(i0,j0)=0
           AH_g0_yy(i0,j0)=0
           AH_g0_yz(i0,j0)=0
           AH_g0_psi(i0,j0)=0
            return
        else


          is=max(1,i-2)
          ie=min(Nx,i+3)
          js=max(1,j-2)
          je=min(Ny,j+3)
          ks=max(1,k-2)
          ke=min(Nz,k+3)


          is=i
          ie=i+1
          js=j
          je=j+1
          ks=k
          ke=k+1

          if (do_ex.ne.0) then
             do i1=is,ie
                do j1=js,je
                   do k1=ks,ke
                      if (chr(i1,j1,k1).eq.ex) then
                    write(*,*) ' calc_AHmetric0: pt i1,j1,k1 is excised'
                         write(*,*) ' i0,j0=',i0,j0
                         write(*,*) ' AH_Nchi,AH_Nphi=',AH_Nchi,AH_Nphi
                         write(*,*) ' AH_chi,AH_phi=',AH_chi,AH_phi
                         write(*,*) ' x0,y0,z0=',x0,y0,z0
                         write(*,*) ' i,j,k=',i,j,k
                         write(*,*) ' xi,yj,zk=',x(i),y(j),z(k)
                         write(*,*) ' chr(i,j,k)=',chr(i,j,k)
                         write(*,*) ' i1,j1,k1=',i1,j1,k1
                         write(*,*) ' xi1,yj1,zk1,AH_R(i0,j0)='
     &                               ,x(i1),y(j1),z(k1),AH_R(i0,j0)
                         write(*,*) ' chr(i1,j1,k1)=',chr(i1,j1,k1)
                    write(*,*) ' rho1=',sqrt(x(i1)**2+y(j1)**2+z(k1)**2)
                         write(*,*) ' Nx,Ny,Nz=',Nx,Ny,Nz
                                    AH_g0_tt(i0,j0)=0
                                    AH_g0_tx(i0,j0)=0
                                    AH_g0_ty(i0,j0)=0
                                    AH_g0_tz(i0,j0)=0
                                    AH_g0_xx(i0,j0)=0
                                    AH_g0_xy(i0,j0)=0
                                    AH_g0_xy(i0,j0)=0
                                    AH_g0_yy(i0,j0)=0
                                    AH_g0_yz(i0,j0)=0
                                    AH_g0_psi(i0,j0)=0
                         return
                      end if
                   end do
                end do
             end do
          end if

             xi=x(i)
             yj=y(j)
             zk=z(k)
             xip1=x(i+1)
             yjp1=y(j+1)
             zkp1=z(k+1)
             rhoijk=sqrt(xi**2+yj**2+zk**2)
             rhoip1jk=sqrt(xip1**2+yj**2+zk**2)
             rhoijp1k=sqrt(xi**2+yjp1**2+zk**2)
             rhoip1jp1k=sqrt(xip1**2+yjp1**2+zk**2)
             rhoijkp1=sqrt(xi**2+yj**2+zkp1**2)
             rhoip1jkp1=sqrt(xip1**2+yj**2+zkp1**2)
             rhoijp1kp1=sqrt(xi**2+yjp1**2+zkp1**2)
             rhoip1jp1kp1=sqrt(xip1**2+yjp1**2+zkp1**2)


        !metric components of pure AdS in Cartesian coordinates
        g0_tt_adsijk =-(4*rhoijk**2+L**2*(-1+rhoijk**2)**2)
     &               /L**2/(-1+rhoijk**2)**2
        g0_tt_adsip1jk =-(4*rhoip1jk**2+L**2*(-1+rhoip1jk**2)**2)
     &               /L**2/(-1+rhoip1jk**2)**2
        g0_tt_adsijp1k =-(4*rhoijp1k**2+L**2*(-1+rhoijp1k**2)**2)
     &               /L**2/(-1+rhoijp1k**2)**2
        g0_tt_adsip1jp1k =-(4*rhoip1jp1k**2+L**2*(-1+rhoip1jp1k**2)**2)
     &               /L**2/(-1+rhoip1jp1k**2)**2
        g0_tt_adsijkp1 =-(4*rhoijkp1**2+L**2*(-1+rhoijkp1**2)**2)
     &               /L**2/(-1+rhoijkp1**2)**2
        g0_tt_adsip1jkp1 =-(4*rhoip1jkp1**2+L**2*(-1+rhoip1jkp1**2)**2)
     &               /L**2/(-1+rhoip1jkp1**2)**2
        g0_tt_adsijp1kp1 =-(4*rhoijp1kp1**2+L**2*(-1+rhoijp1kp1**2)**2)
     &               /L**2/(-1+rhoijp1kp1**2)**2
        g0_tt_adsip1jp1kp1 =-(4*rhoip1jp1kp1**2
     &                        +L**2*(-1+rhoip1jp1kp1**2)**2)
     &               /L**2/(-1+rhoip1jp1kp1**2)**2


        g0_tx_adsijk =0
        g0_tx_adsip1jk =0
        g0_tx_adsijp1k =0
        g0_tx_adsip1jp1k =0
        g0_tx_adsijkp1 =0
        g0_tx_adsip1jkp1 =0
        g0_tx_adsijp1kp1 =0
        g0_tx_adsip1jp1kp1 =0

        g0_ty_adsijk =0
        g0_ty_adsip1jk =0
        g0_ty_adsijp1k =0
        g0_ty_adsip1jp1k =0
        g0_ty_adsijkp1 =0
        g0_ty_adsip1jkp1 =0
        g0_ty_adsijp1kp1 =0
        g0_ty_adsip1jp1kp1 =0

        g0_tz_adsijk =0
        g0_tz_adsip1jk =0
        g0_tz_adsijp1k =0
        g0_tz_adsip1jp1k =0
        g0_tz_adsijkp1 =0
        g0_tz_adsip1jkp1 =0
        g0_tz_adsijp1kp1 =0
        g0_tz_adsip1jp1kp1 =0


        g0_xx_adsijk =(8*(-1+L**2)*(xi**2-yj**2-zk**2)
     &           +8*rhoijk**2+4*L**2*(1+rhoijk**4))
     &          /(-1+rhoijk**2)**2/(4*rhoijk**2+L**2*(-1+rhoijk**2)**2)
        g0_xx_adsip1jk =(8*(-1+L**2)*(xip1**2-yj**2-zk**2)
     &           +8*rhoip1jk**2+4*L**2*(1+rhoip1jk**4))
     &   /(-1+rhoip1jk**2)**2/(4*rhoip1jk**2+L**2*(-1+rhoip1jk**2)**2)
        g0_xx_adsijp1k =(8*(-1+L**2)*(xi**2-yjp1**2-zk**2)
     &           +8*rhoijp1k**2+4*L**2*(1+rhoijp1k**4))
     &   /(-1+rhoijp1k**2)**2/(4*rhoijp1k**2+L**2*(-1+rhoijp1k**2)**2)
        g0_xx_adsip1jp1k =(8*(-1+L**2)*(xip1**2-yjp1**2-zk**2)
     &           +8*rhoip1jp1k**2+4*L**2*(1+rhoip1jp1k**4))
     &   /(-1+rhoip1jp1k**2)**2/(4*rhoip1jp1k**2
     &                           +L**2*(-1+rhoip1jp1k**2)**2)
        g0_xx_adsijkp1 =(8*(-1+L**2)*(xi**2-yj**2-zkp1**2)
     &           +8*rhoijkp1**2+4*L**2*(1+rhoijkp1**4))
     &   /(-1+rhoijkp1**2)**2/(4*rhoijkp1**2+L**2*(-1+rhoijkp1**2)**2)
        g0_xx_adsip1jkp1 =(8*(-1+L**2)*(xip1**2-yj**2-zkp1**2)
     &           +8*rhoip1jkp1**2+4*L**2*(1+rhoip1jkp1**4))
     &   /(-1+rhoip1jkp1**2)**2/(4*rhoip1jkp1**2
     &                           +L**2*(-1+rhoip1jkp1**2)**2)
        g0_xx_adsijp1kp1 =(8*(-1+L**2)*(xi**2-yjp1**2-zkp1**2)
     &           +8*rhoijp1kp1**2+4*L**2*(1+rhoijp1kp1**4))
     &   /(-1+rhoijp1kp1**2)**2/(4*rhoijp1kp1**2
     &                           +L**2*(-1+rhoijp1kp1**2)**2)
        g0_xx_adsip1jp1kp1 =(8*(-1+L**2)*(xip1**2-yjp1**2-zkp1**2)
     &           +8*rhoip1jp1kp1**2+4*L**2*(1+rhoip1jp1kp1**4))
     &   /(-1+rhoip1jp1kp1**2)**2/(4*rhoip1jp1kp1**2
     &                           +L**2*(-1+rhoip1jp1kp1**2)**2)

        g0_xy_adsijk =(16 *(-1 + L**2) *xi* yj)
     &              /((-1 + rhoijk**2)**2
     &               *(4 *rhoijk**2 +L**2 *(-1 +rhoijk**2)**2))
        g0_xy_adsip1jk =(16 *(-1 + L**2) *xip1* yj)
     &              /((-1 + rhoip1jk**2)**2
     &               *(4 *rhoip1jk**2 +L**2 *(-1 +rhoip1jk**2)**2))
        g0_xy_adsijp1k =(16 *(-1 + L**2) *xi* yjp1)
     &              /((-1 + rhoijp1k**2)**2
     &               *(4 *rhoijp1k**2 +L**2 *(-1 +rhoijp1k**2)**2))
        g0_xy_adsip1jp1k =(16 *(-1 + L**2) *xip1* yjp1)
     &              /((-1 + rhoip1jp1k**2)**2
     &               *(4 *rhoip1jp1k**2 +L**2 *(-1 +rhoip1jp1k**2)**2))
        g0_xy_adsijkp1 =(16 *(-1 + L**2) *xi* yj)
     &              /((-1 + rhoijkp1**2)**2
     &               *(4 *rhoijkp1**2 +L**2 *(-1 +rhoijkp1**2)**2))
        g0_xy_adsip1jkp1 =(16 *(-1 + L**2) *xip1* yj)
     &              /((-1 + rhoip1jkp1**2)**2
     &               *(4 *rhoip1jkp1**2 +L**2 *(-1 +rhoip1jkp1**2)**2))
        g0_xy_adsijp1kp1 =(16 *(-1 + L**2) *xi* yjp1)
     &              /((-1 + rhoijp1kp1**2)**2
     &               *(4 *rhoijp1kp1**2 +L**2 *(-1 +rhoijp1kp1**2)**2))
        g0_xy_adsip1jp1kp1 =(16 *(-1 + L**2) *xip1* yjp1)
     &              /((-1 + rhoip1jp1kp1**2)**2
     &          *(4 *rhoip1jp1kp1**2 +L**2 *(-1 +rhoip1jp1kp1**2)**2))

        g0_xz_adsijk =(16 *(-1 + L**2) *xi* zk)
     &              /((-1 + rhoijk**2)**2
     &               *(4 *rhoijk**2 +L**2 *(-1 +rhoijk**2)**2))
        g0_xz_adsip1jk =(16 *(-1 + L**2) *xip1* zk)
     &              /((-1 + rhoip1jk**2)**2
     &               *(4 *rhoip1jk**2 +L**2 *(-1 +rhoip1jk**2)**2))
        g0_xz_adsijp1k =(16 *(-1 + L**2) *xi* zk)
     &              /((-1 + rhoijp1k**2)**2
     &               *(4 *rhoijp1k**2 +L**2 *(-1 +rhoijp1k**2)**2))
        g0_xz_adsip1jp1k =(16 *(-1 + L**2) *xip1* zk)
     &              /((-1 + rhoip1jp1k**2)**2
     &               *(4 *rhoip1jp1k**2 +L**2 *(-1 +rhoip1jp1k**2)**2))
        g0_xz_adsijkp1 =(16 *(-1 + L**2) *xi* zkp1)
     &              /((-1 + rhoijkp1**2)**2
     &               *(4 *rhoijkp1**2 +L**2 *(-1 +rhoijkp1**2)**2))
        g0_xz_adsip1jkp1 =(16 *(-1 + L**2) *xip1* zkp1)
     &              /((-1 + rhoip1jkp1**2)**2
     &               *(4 *rhoip1jkp1**2 +L**2 *(-1 +rhoip1jkp1**2)**2))
        g0_xz_adsijp1kp1 =(16 *(-1 + L**2) *xi* zkp1)
     &              /((-1 + rhoijp1kp1**2)**2
     &               *(4 *rhoijp1kp1**2 +L**2 *(-1 +rhoijp1kp1**2)**2))
        g0_xz_adsip1jp1kp1 =(16 *(-1 + L**2) *xip1* zkp1)
     &              /((-1 + rhoip1jp1kp1**2)**2
     &          *(4 *rhoip1jp1kp1**2 +L**2 *(-1 +rhoip1jp1kp1**2)**2))


        g0_yy_adsijk =(4*(4*(xi**2+zk**2)+L**2*(xi**4+(1+yj**2)**2
     &              +2*(-1+yj**2)*zk**2+zk**4
     &              +2*xi**2*(-1+yj**2+zk**2))))
     &              /(L**2*(-1+rhoijk**2)**4
     &                +4*(-1+rhoijk**2)**2*(rhoijk**2))
        g0_yy_adsip1jk =(4*(4*(xip1**2+zk**2)+L**2*(xip1**4+(1+yj**2)**2
     &              +2*(-1+yj**2)*zk**2+zk**4
     &              +2*xip1**2*(-1+yj**2+zk**2))))
     &              /(L**2*(-1+rhoip1jk**2)**4
     &                +4*(-1+rhoip1jk**2)**2*(rhoip1jk**2))
        g0_yy_adsijp1k =(4*(4*(xi**2+zk**2)+L**2*(xi**4+(1+yjp1**2)**2
     &              +2*(-1+yjp1**2)*zk**2+zk**4
     &              +2*xi**2*(-1+yjp1**2+zk**2))))
     &              /(L**2*(-1+rhoijp1k**2)**4
     &                +4*(-1+rhoijp1k**2)**2*(rhoijp1k**2))
        g0_yy_adsip1jp1k =(4*(4*(xip1**2+zk**2)
     &                       +L**2*(xip1**4+(1+yjp1**2)**2
     &              +2*(-1+yjp1**2)*zk**2+zk**4
     &              +2*xip1**2*(-1+yjp1**2+zk**2))))
     &              /(L**2*(-1+rhoip1jp1k**2)**4
     &                +4*(-1+rhoip1jp1k**2)**2*(rhoip1jp1k**2))
        g0_yy_adsijkp1 =(4*(4*(xi**2+zkp1**2)+L**2*(xi**4+(1+yj**2)**2
     &              +2*(-1+yj**2)*zkp1**2+zkp1**4
     &              +2*xi**2*(-1+yj**2+zkp1**2))))
     &              /(L**2*(-1+rhoijkp1**2)**4
     &                +4*(-1+rhoijkp1**2)**2*(rhoijkp1**2))
        g0_yy_adsip1jkp1 =(4*(4*(xip1**2+zkp1**2)
     &                      +L**2*(xip1**4+(1+yj**2)**2
     &              +2*(-1+yj**2)*zkp1**2+zkp1**4
     &              +2*xip1**2*(-1+yj**2+zkp1**2))))
     &              /(L**2*(-1+rhoip1jkp1**2)**4
     &                +4*(-1+rhoip1jkp1**2)**2*(rhoip1jkp1**2))
        g0_yy_adsijp1kp1 =(4*(4*(xi**2+zkp1**2)
     &                        +L**2*(xi**4+(1+yjp1**2)**2
     &              +2*(-1+yjp1**2)*zkp1**2+zkp1**4
     &              +2*xi**2*(-1+yjp1**2+zkp1**2))))
     &              /(L**2*(-1+rhoijp1kp1**2)**4
     &                +4*(-1+rhoijp1kp1**2)**2*(rhoijp1kp1**2))
        g0_yy_adsip1jp1kp1 =(4*(4*(xip1**2+zkp1**2)
     &                        +L**2*(xip1**4+(1+yjp1**2)**2
     &              +2*(-1+yjp1**2)*zkp1**2+zkp1**4
     &              +2*xip1**2*(-1+yjp1**2+zkp1**2))))
     &              /(L**2*(-1+rhoip1jp1kp1**2)**4
     &                +4*(-1+rhoip1jp1kp1**2)**2*(rhoip1jp1kp1**2))

        g0_yz_adsijk =(16 *(-1 + L**2) *yj* zk)
     &              /((-1 + rhoijk**2)**2
     &               *(4 *rhoijk**2 +L**2 *(-1 +rhoijk**2)**2))
        g0_yz_adsip1jk =(16 *(-1 + L**2) *yj* zk)
     &              /((-1 + rhoip1jk**2)**2
     &               *(4 *rhoip1jk**2 +L**2 *(-1 +rhoip1jk**2)**2))
        g0_yz_adsijp1k =(16 *(-1 + L**2) *yjp1* zk)
     &              /((-1 + rhoijp1k**2)**2
     &               *(4 *rhoijp1k**2 +L**2 *(-1 +rhoijp1k**2)**2))
        g0_yz_adsip1jp1k =(16 *(-1 + L**2) *yjp1* zk)
     &              /((-1 + rhoip1jp1k**2)**2
     &               *(4 *rhoip1jp1k**2 +L**2 *(-1 +rhoip1jp1k**2)**2))
        g0_yz_adsijkp1 =(16 *(-1 + L**2) *yj* zkp1)
     &              /((-1 + rhoijkp1**2)**2
     &               *(4 *rhoijkp1**2 +L**2 *(-1 +rhoijkp1**2)**2))
        g0_yz_adsip1jkp1 =(16 *(-1 + L**2) *yj* zkp1)
     &              /((-1 + rhoip1jkp1**2)**2
     &               *(4 *rhoip1jkp1**2 +L**2 *(-1 +rhoip1jkp1**2)**2))
        g0_yz_adsijp1kp1 =(16 *(-1 + L**2) *yjp1* zkp1)
     &              /((-1 + rhoijp1kp1**2)**2
     &               *(4 *rhoijp1kp1**2 +L**2 *(-1 +rhoijp1kp1**2)**2))
        g0_yz_adsip1jp1kp1 =(16 *(-1 + L**2) *yjp1* zkp1)
     &              /((-1 + rhoip1jp1kp1**2)**2
     &          *(4 *rhoip1jp1kp1**2 +L**2 *(-1 +rhoip1jp1kp1**2)**2))

        g0_psi_adsijk=(4*(4*(xi**2+yj**2)
     &                    +L**2*((-1+xi**2+yj**2)**2
     &              +2*(1+xi**2+yj**2)*zk**2+zk**4)))
     &              /(L**2*(-1+rhoijk**2)**4
     &              +4*(-1+rhoijk**2)**2*(rhoijk**2))
        g0_psi_adsip1jk=(4*(4*(xip1**2+yj**2)
     &                    +L**2*((-1+xip1**2+yj**2)**2
     &              +2*(1+xip1**2+yj**2)*zk**2+zk**4)))
     &              /(L**2*(-1+rhoip1jk**2)**4
     &              +4*(-1+rhoip1jk**2)**2*(rhoip1jk**2))
        g0_psi_adsijp1k=(4*(4*(xi**2+yjp1**2)
     &                    +L**2*((-1+xi**2+yjp1**2)**2
     &              +2*(1+xi**2+yjp1**2)*zk**2+zk**4)))
     &              /(L**2*(-1+rhoijp1k**2)**4
     &              +4*(-1+rhoijp1k**2)**2*(rhoijp1k**2))
        g0_psi_adsip1jp1k=(4*(4*(xip1**2+yjp1**2)
     &                    +L**2*((-1+xip1**2+yjp1**2)**2
     &              +2*(1+xip1**2+yjp1**2)*zk**2+zk**4)))
     &              /(L**2*(-1+rhoip1jp1k**2)**4
     &              +4*(-1+rhoip1jp1k**2)**2*(rhoip1jp1k**2))
        g0_psi_adsijkp1=(4*(4*(xi**2+yj**2)
     &                    +L**2*((-1+xi**2+yj**2)**2
     &              +2*(1+xi**2+yj**2)*zkp1**2+zkp1**4)))
     &              /(L**2*(-1+rhoijkp1**2)**4
     &              +4*(-1+rhoijkp1**2)**2*(rhoijkp1**2))
        g0_psi_adsip1jkp1=(4*(4*(xip1**2+yj**2)
     &                    +L**2*((-1+xip1**2+yj**2)**2
     &              +2*(1+xip1**2+yj**2)*zkp1**2+zkp1**4)))
     &              /(L**2*(-1+rhoip1jkp1**2)**4
     &              +4*(-1+rhoip1jkp1**2)**2*(rhoip1jkp1**2))
        g0_psi_adsijp1kp1=(4*(4*(xi**2+yjp1**2)
     &                    +L**2*((-1+xi**2+yjp1**2)**2
     &              +2*(1+xi**2+yjp1**2)*zkp1**2+zkp1**4)))
     &              /(L**2*(-1+rhoijp1kp1**2)**4
     &              +4*(-1+rhoijp1kp1**2)**2*(rhoijp1kp1**2))
        g0_psi_adsip1jp1kp1=(4*(4*(xip1**2+yjp1**2)
     &                    +L**2*((-1+xip1**2+yjp1**2)**2
     &              +2*(1+xip1**2+yjp1**2)*zkp1**2+zkp1**4)))
     &              /(L**2*(-1+rhoip1jp1kp1**2)**4
     &              +4*(-1+rhoip1jp1kp1**2)**2*(rhoip1jp1kp1**2))

       !full metric in Cartesian coordinates
        g0_tt_ijk=      g0_tt_adsijk       +gb_tt_n(i,j,k)
        g0_tt_ip1jk=    g0_tt_adsip1jk     +gb_tt_n(i+1,j,k)
        g0_tt_ijp1k=    g0_tt_adsijp1k     +gb_tt_n(i,j+1,k)
        g0_tt_ip1jp1k=  g0_tt_adsip1jp1k   +gb_tt_n(i+1,j+1,k)
        g0_tt_ijkp1=    g0_tt_adsijkp1     +gb_tt_n(i,j,k+1)
        g0_tt_ip1jkp1=  g0_tt_adsip1jkp1   +gb_tt_n(i+1,j,k+1)
        g0_tt_ijp1kp1=  g0_tt_adsijp1kp1   +gb_tt_n(i,j+1,k+1)
        g0_tt_ip1jp1kp1=g0_tt_adsip1jp1kp1 +gb_tt_n(i+1,j+1,k+1)

        g0_tx_ijk=      g0_tx_adsijk       +gb_tx_n(i,j,k)
        g0_tx_ip1jk=    g0_tx_adsip1jk     +gb_tx_n(i+1,j,k)
        g0_tx_ijp1k=    g0_tx_adsijp1k     +gb_tx_n(i,j+1,k)
        g0_tx_ip1jp1k=  g0_tx_adsip1jp1k   +gb_tx_n(i+1,j+1,k)
        g0_tx_ijkp1=    g0_tx_adsijkp1     +gb_tx_n(i,j,k+1)
        g0_tx_ip1jkp1=  g0_tx_adsip1jkp1   +gb_tx_n(i+1,j,k+1)
        g0_tx_ijp1kp1=  g0_tx_adsijp1kp1   +gb_tx_n(i,j+1,k+1)
        g0_tx_ip1jp1kp1=g0_tx_adsip1jp1kp1 +gb_tx_n(i+1,j+1,k+1)

        g0_ty_ijk=      g0_ty_adsijk       +gb_ty_n(i,j,k)
        g0_ty_ip1jk=    g0_ty_adsip1jk     +gb_ty_n(i+1,j,k)
        g0_ty_ijp1k=    g0_ty_adsijp1k     +gb_ty_n(i,j+1,k)
        g0_ty_ip1jp1k=  g0_ty_adsip1jp1k   +gb_ty_n(i+1,j+1,k)
        g0_ty_ijkp1=    g0_ty_adsijkp1     +gb_ty_n(i,j,k+1)
        g0_ty_ip1jkp1=  g0_ty_adsip1jkp1   +gb_ty_n(i+1,j,k+1)
        g0_ty_ijp1kp1=  g0_ty_adsijp1kp1   +gb_ty_n(i,j+1,k+1)
        g0_ty_ip1jp1kp1=g0_ty_adsip1jp1kp1 +gb_ty_n(i+1,j+1,k+1)

        g0_tz_ijk=      g0_tz_adsijk       +gb_tz_n(i,j,k)
        g0_tz_ip1jk=    g0_tz_adsip1jk     +gb_tz_n(i+1,j,k)
        g0_tz_ijp1k=    g0_tz_adsijp1k     +gb_tz_n(i,j+1,k)
        g0_tz_ip1jp1k=  g0_tz_adsip1jp1k   +gb_tz_n(i+1,j+1,k)
        g0_tz_ijkp1=    g0_tz_adsijkp1     +gb_tz_n(i,j,k+1)
        g0_tz_ip1jkp1=  g0_tz_adsip1jkp1   +gb_tz_n(i+1,j,k+1)
        g0_tz_ijp1kp1=  g0_tz_adsijp1kp1   +gb_tz_n(i,j+1,k+1)
        g0_tz_ip1jp1kp1=g0_tz_adsip1jp1kp1 +gb_tz_n(i+1,j+1,k+1)

        g0_xx_ijk=      g0_xx_adsijk       +gb_xx_n(i,j,k)
        g0_xx_ip1jk=    g0_xx_adsip1jk     +gb_xx_n(i+1,j,k)
        g0_xx_ijp1k=    g0_xx_adsijp1k     +gb_xx_n(i,j+1,k)
        g0_xx_ip1jp1k=  g0_xx_adsip1jp1k   +gb_xx_n(i+1,j+1,k)
        g0_xx_ijkp1=    g0_xx_adsijkp1     +gb_xx_n(i,j,k+1)
        g0_xx_ip1jkp1=  g0_xx_adsip1jkp1   +gb_xx_n(i+1,j,k+1)
        g0_xx_ijp1kp1=  g0_xx_adsijp1kp1   +gb_xx_n(i,j+1,k+1)
        g0_xx_ip1jp1kp1=g0_xx_adsip1jp1kp1 +gb_xx_n(i+1,j+1,k+1)

        g0_xy_ijk=      g0_xy_adsijk       +gb_xy_n(i,j,k)
        g0_xy_ip1jk=    g0_xy_adsip1jk     +gb_xy_n(i+1,j,k)
        g0_xy_ijp1k=    g0_xy_adsijp1k     +gb_xy_n(i,j+1,k)
        g0_xy_ip1jp1k=  g0_xy_adsip1jp1k   +gb_xy_n(i+1,j+1,k)
        g0_xy_ijkp1=    g0_xy_adsijkp1     +gb_xy_n(i,j,k+1)
        g0_xy_ip1jkp1=  g0_xy_adsip1jkp1   +gb_xy_n(i+1,j,k+1)
        g0_xy_ijp1kp1=  g0_xy_adsijp1kp1   +gb_xy_n(i,j+1,k+1)
        g0_xy_ip1jp1kp1=g0_xy_adsip1jp1kp1 +gb_xy_n(i+1,j+1,k+1)

        g0_xz_ijk=      g0_xz_adsijk       +gb_xz_n(i,j,k)
        g0_xz_ip1jk=    g0_xz_adsip1jk     +gb_xz_n(i+1,j,k)
        g0_xz_ijp1k=    g0_xz_adsijp1k     +gb_xz_n(i,j+1,k)
        g0_xz_ip1jp1k=  g0_xz_adsip1jp1k   +gb_xz_n(i+1,j+1,k)
        g0_xz_ijkp1=    g0_xz_adsijkp1     +gb_xz_n(i,j,k+1)
        g0_xz_ip1jkp1=  g0_xz_adsip1jkp1   +gb_xz_n(i+1,j,k+1)
        g0_xz_ijp1kp1=  g0_xz_adsijp1kp1   +gb_xz_n(i,j+1,k+1)
        g0_xz_ip1jp1kp1=g0_xz_adsip1jp1kp1 +gb_xz_n(i+1,j+1,k+1)

        g0_yy_ijk=      g0_yy_adsijk       +gb_yy_n(i,j,k)
        g0_yy_ip1jk=    g0_yy_adsip1jk     +gb_yy_n(i+1,j,k)
        g0_yy_ijp1k=    g0_yy_adsijp1k     +gb_yy_n(i,j+1,k)
        g0_yy_ip1jp1k=  g0_yy_adsip1jp1k   +gb_yy_n(i+1,j+1,k)
        g0_yy_ijkp1=    g0_yy_adsijkp1     +gb_yy_n(i,j,k+1)
        g0_yy_ip1jkp1=  g0_yy_adsip1jkp1   +gb_yy_n(i+1,j,k+1)
        g0_yy_ijp1kp1=  g0_yy_adsijp1kp1   +gb_yy_n(i,j+1,k+1)
        g0_yy_ip1jp1kp1=g0_yy_adsip1jp1kp1 +gb_yy_n(i+1,j+1,k+1)

        g0_yz_ijk=      g0_yz_adsijk       +gb_yz_n(i,j,k)
        g0_yz_ip1jk=    g0_yz_adsip1jk     +gb_yz_n(i+1,j,k)
        g0_yz_ijp1k=    g0_yz_adsijp1k     +gb_yz_n(i,j+1,k)
        g0_yz_ip1jp1k=  g0_yz_adsip1jp1k   +gb_yz_n(i+1,j+1,k)
        g0_yz_ijkp1=    g0_yz_adsijkp1     +gb_yz_n(i,j,k+1)
        g0_yz_ip1jkp1=  g0_yz_adsip1jkp1   +gb_yz_n(i+1,j,k+1)
        g0_yz_ijp1kp1=  g0_yz_adsijp1kp1   +gb_yz_n(i,j+1,k+1)
        g0_yz_ip1jp1kp1=g0_yz_adsip1jp1kp1 +gb_yz_n(i+1,j+1,k+1)

        g0_psi_ijk=      g0_psi_adsijk       +psi_n(i,j,k)
        g0_psi_ip1jk=    g0_psi_adsip1jk     +psi_n(i+1,j,k)
        g0_psi_ijp1k=    g0_psi_adsijp1k     +psi_n(i,j+1,k)
        g0_psi_ip1jp1k=  g0_psi_adsip1jp1k   +psi_n(i+1,j+1,k)
        g0_psi_ijkp1=    g0_psi_adsijkp1     +psi_n(i,j,k+1)
        g0_psi_ip1jkp1=  g0_psi_adsip1jkp1   +psi_n(i+1,j,k+1)
        g0_psi_ijp1kp1=  g0_psi_adsijp1kp1   +psi_n(i,j+1,k+1)
        g0_psi_ip1jp1kp1=g0_psi_adsip1jp1kp1 +psi_n(i+1,j+1,k+1)

          ! interpolate metric components from 
          ! (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1),(i,j+1,k+1),(i+1,j+1,k+1)
          AH_g0_tt(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_tt_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_tt_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_tt_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_tt_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_tt_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_tt_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_tt_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_tt_ip1jp1kp1

!         write(*,*) "calc_ahmetric0:i0,j0,AH_g0_tt(i0,j0)="
!     &                             ,i0,j0,AH_g0_tt(i0,j0)

          AH_g0_tx(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_tx_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_tx_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_tx_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_tx_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_tx_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_tx_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_tx_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_tx_ip1jp1kp1

          AH_g0_ty(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_ty_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_ty_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_ty_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_ty_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_ty_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_ty_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_ty_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_ty_ip1jp1kp1

          AH_g0_tz(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_tz_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_tz_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_tz_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_tz_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_tz_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_tz_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_tz_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_tz_ip1jp1kp1

          AH_g0_xx(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_xx_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_xx_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_xx_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_xx_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_xx_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_xx_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_xx_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_xx_ip1jp1kp1

          AH_g0_xy(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_xy_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_xy_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_xy_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_xy_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_xy_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_xy_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_xy_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_xy_ip1jp1kp1

          AH_g0_xz(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_xz_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_xz_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_xz_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_xz_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_xz_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_xz_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_xz_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_xz_ip1jp1kp1

          AH_g0_yy(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_yy_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_yy_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_yy_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_yy_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_yy_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_yy_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_yy_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_yy_ip1jp1kp1

          AH_g0_yz(i0,j0)=(1-fx)*(1-fy)*(1-fz)  *g0_yz_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_yz_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_yz_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_yz_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_yz_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_yz_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_yz_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_yz_ip1jp1kp1

          AH_g0_psi(i0,j0)=(1-fx)*(1-fy)*(1-fz) *g0_psi_ijk+
     &                    (  fx)*(1-fy)*(1-fz)  *g0_psi_ip1jk+
     &                    (1-fx)*(  fy)*(1-fz)  *g0_psi_ijp1k+
     &                    (  fx)*(  fy)*(1-fz)  *g0_psi_ip1jp1k+
     &                    (1-fx)*(1-fy)*(  fz)  *g0_psi_ijkp1+
     &                    (  fx)*(1-fy)*(  fz)  *g0_psi_ip1jkp1+
     &                    (1-fx)*(  fy)*(  fz)  *g0_psi_ijp1kp1+
     &                    (  fx)*(  fy)*(  fz)  *g0_psi_ip1jp1kp1

        end if

        return
        end

!-----------------------------------------------------------------------------------------------

c-----------------------------------------------------------------------
c Kreiss-Oliger-like smoothing of R (or theta, or whatever)
c
c + removes *isolated* "spikes" that sometimes arise in late time
c distorted holes prior to merger with rem_spikes option
c-----------------------------------------------------------------------
        subroutine smooth_ah_r(AH_R,AH_w1,AH_eps,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_eps
        real*8 AH_w1(AH_Nchi,AH_Nphi)
        
        integer i,j,jp1,jp2,jm1,jm2,j0
        real*8 r0,r1,rip1,rip2,rim1,rim2,rjp1,rjp2,rjm1,rjm2,r
        real*8 avg_var,v1,v2,v3,v4

        logical rem_spikes
        parameter (rem_spikes=.true.)
        real*8 max_var
        parameter (max_var=5)

        ! initialize fixed-size variables 
        data i,j,jp1,jp2,jm1,jm2,j0/0,0,0,0,0,0,0/

        data r0,r1,rip1,rip2/0.0,0.0,0.0,0.0/
        data rim1,rim2,rjp1,rjp2/0.0,0.0,0.0,0.0/
        data rjm1,rjm2,r/0.0,0.0,0.0/
        data avg_var,v1,v2,v3,v4/0.0,0.0,0.0,0.0,0.0/
 
        !--------------------------------------------------------------

        if (AH_Nphi.gt.1.and.rem_spikes) then
           avg_var=0
           do i=1,AH_Nchi
              do j=1,AH_Nphi
                 AH_w1(i,j)=0
              end do
           end do
           do i=2,AH_Nchi-1
              do j=1,AH_Nphi-1
                 jm1=j-1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 jp1=j+1
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 r=AH_R(i,j)
                 v1=(r-AH_R(i+1,j))**2
                 v2=(r-AH_R(i-1,j))**2
                 v3=(r-AH_R(i,jp1))**2
                 v4=(r-AH_R(i,jm1))**2
                 !-----------------------------------------------------
                 ! by subtracting the max, we are only sensitive
                 ! to isolated spikes (i.e, the points adjacent to 
                 ! the spike will now have a small variance)
                 !-----------------------------------------------------
                 AH_w1(i,j)=sqrt(v1+v2+v3+v4-max(v1,v2,v3,v4))
                 avg_var=avg_var+AH_w1(i,j)
              end do
           end do
           avg_var=avg_var/(AH_Nphi-1)/(AH_Nchi-2)
           do i=2,AH_Nchi-1
              do j=1,AH_Nphi-1
                 jm1=j-1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 jp1=j+1
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 if (AH_w1(i,j)/avg_var.gt.max_var.and.
     &               AH_w1(i+1,j)/avg_var.lt.max_var.and.
     &               AH_w1(i-1,j)/avg_var.lt.max_var.and.
     &               AH_w1(i,jp1)/avg_var.lt.max_var.and.
     &               AH_w1(i,jm1)/avg_var.lt.max_var) then
                    AH_R(i,j)=0.25d0*(AH_R(i+1,j)+AH_R(i-1,j)+
     &                                AH_R(i,jp1)+AH_R(i,jm1))
                 end if
              end do
           end do
           do i=1,AH_Nchi
              AH_R(i,AH_Nphi)=AH_R(i,1)
           end do
        end if

        r0=0
        r1=0
        do i=1,AH_Nphi
           r0=r0+AH_R(1,i)
           r1=r1+AH_R(AH_Nchi,i)
        end do
        r0=r0/AH_Nphi
        r1=r1/AH_Nphi
        do i=1,AH_Nphi
           AH_R(1,i)=r0
           AH_R(AH_Nchi,i)=r1
        end do

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              AH_w1(i,j)=AH_R(i,j)
           end do
        end do

        do i=2,AH_Nchi-1
           do j=1,AH_Nphi
              r=AH_w1(i,j)
              rip1=AH_w1(i+1,j)
              rim1=AH_w1(i-1,j)
              if (i.eq.(AH_Nchi-1)) then
                 rip2=r
              else
                 rip2=AH_w1(i+2,j)
              end if
              if (i.eq.2) then
                 rim2=r
              else
                 rim2=AH_w1(i-2,j)
              end if
              if (AH_Nphi.lt.5) then
                 jp1=j
                 jm1=j
                 jp2=j
                 jm2=j
              else
                 jp1=j+1
                 jp2=j+2
                 jm1=j-1
                 jm2=j-2
                 if (jp1.gt.AH_Nphi) jp1=jp1-AH_Nphi+1
                 if (jp2.gt.AH_Nphi) jp2=jp2-AH_Nphi+1
                 if (jm1.lt.1) jm1=jm1+AH_Nphi-1
                 if (jm2.lt.1) jm2=jm2+AH_Nphi-1
              end if
              rjp1=AH_w1(i,jp1)
              rjp2=AH_w1(i,jp2)
              rjm1=AH_w1(i,jm1)
              rjm2=AH_w1(i,jm2)

              AH_R(i,j)=r-AH_eps/16*(
     &          rip2+rim2-4*(rip1+rim1)+6*r +
     &          rjp2+rjm2-4*(rjp1+rjm1)+6*r)
           end do
        end do
        
        do i=1,AH_Nchi
           AH_R(i,1)=AH_R(i,AH_Nphi)
        end do

        return 
        end

c-----------------------------------------------------------------------
c enforce on-axis regularity of R (or theta, or whatever)
c-----------------------------------------------------------------------
        subroutine reg_ah_r(AH_R,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi)
        
        integer i,j
        real*8 r0,r1

        ! initialize fixed-size variables 
        data i,j/0,0/
        data r0,r1/0.0,0.0/

        !--------------------------------------------------------------

        !HB!
        if (AH_Nphi.eq.1) then ! for the phi=0 point on the chi=0,PI poles, set by zero-derivative extrapolation along phi
           AH_r(1,1)=(AH_r(2,1)*4-AH_r(3,1))/3
           AH_r(AH_Nchi,1)=(AH_r(AH_Nchi-1,1)*4-AH_r(AH_Nchi-2,1))/3
        else ! for the rest of the phi points on the chi=0,PI poles, set by average value along phi
           r0=0
           r1=0
           do j=1,AH_Nphi
              r0=r0+AH_R(1,j)
              r1=r1+AH_R(AH_Nchi,j)
           end do
           r0=r0/AH_Nphi
           r1=r1/AH_Nphi
           do j=1,AH_Nphi
              AH_r(1,j)=r0
              AH_r(AH_Nchi,j)=r1
           end do
        end if

        return 
        end

c-----------------------------------------------------------------------
c the following changes AH_xc so that median(AH_R) is zero, and
c adjusts AH_R conversely. 
c
c NOTE: the manner in which AH_R is adjusted assumes that the change 
c in AH_xc is *small* !
c-----------------------------------------------------------------------
        subroutine adjust_ah_xc(AH_R,AH_xc,AH_Nchi,AH_Nphi,
     &                          dx,dy,dz,axisym)
        implicit none
        integer axisym
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3)
        
        integer i,j
        real*8 dahchi,dahphi,AH_chi,AH_phi
        real*8 xmin,xmax,ymin,ymax,zmin,zmax,x0,y0,z0,dx0,dy0,dz0
        real*8 dx,dy,dz
        real*8 csx,csy,csz

        real*8 PI
        parameter (PI=3.141592653589793d0)

!LOOKING!
        logical ltrace
        parameter (ltrace=.false.)
!        parameter (ltrace=.true.)

        ! initialize fixed-size variables
        data i,j/0,0/

        data dahchi,dahphi,AH_chi,AH_phi/0.0,0.0,0.0,0.0/
        data xmin,xmax,ymin,ymax/0.0,0.0,0.0,0.0/
        data zmin,zmax,x0,y0,z0/0.0,0.0,0.0,0.0,0.0/
        data dx0,dy0,dz0/0.0,0.0,0.0/
        data csx,csy,csz/0.0,0.0,0.0/

        !--------------------------------------------------------------

        xmin=0
        ymin=0
        zmin=0
        xmax=0
        ymax=0
        zmax=0

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        do i=1,AH_Nchi
           AH_chi=(i-1)*dahchi
           do j=1,AH_Nphi
              AH_phi=(j-1)*dahphi

              !x0,y0,z0 cartesian coordinates of point on AH, wrt center of AH
              !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
!              x0=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!              y0=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!              z0=AH_R(i,j)*cos(AH_chi)+AH_xc(3)

              x0=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
              y0=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
              z0=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

              xmin=min(xmin,x0)
              xmax=max(xmax,x0)
              ymin=min(ymin,y0)
              ymax=max(ymax,y0)
              zmin=min(zmin,z0)
              zmax=max(zmax,z0)
           end do
        end do

        ! for now, don't move xc (yc,zc) when xc=0 (yc=0,zc=0)
        if (AH_xc(1).eq.0) then
          dx0=0
        else
          dx0=(xmax+xmin)/2
        end if
        if (AH_xc(2).eq.0) then
          dy0=0
        else
          dy0=(ymax+ymin)/2
        end if
        if (AH_xc(3).eq.0) then
          dz0=0
        else
          dz0=(zmax+zmin)/2
        end if

        AH_xc(1)=AH_xc(1)+dx0
        AH_xc(2)=AH_xc(2)+dy0
        AH_xc(3)=AH_xc(3)+dz0

        if (ltrace) then
          write(*,*) '-----------------------------------'
          write(*,*) 'dx0,AH_xc(1)=',dx0,AH_xc(1)
          write(*,*) 'dy0,AH_xc(2)=',dy0,AH_xc(2)
          write(*,*) 'dz0,AH_xc(3)=',dz0,AH_xc(3)
        end if

        do i=1,AH_Nchi
           AH_chi=(i-1)*dahchi
           do j=1,AH_Nphi
              AH_phi=(j-1)*dahphi

              !x0,y0,z0 cartesian coordinates of point on AH, wrt center of AH
              !AH_R,AH_chi,AH_phi polar coordinates of point on AH, wrt center of AH
!              x0=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(1)
!              y0=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(2)
!              z0=AH_R(i,j)*cos(AH_chi)+AH_xc(3)

              x0=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
              y0=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
              z0=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

              csx=x0/AH_R(i,j)
              csy=y0/AH_R(i,j)
              csz=z0/AH_R(i,j)
              AH_R(i,j)=sqrt((abs(x0)-csx*dx0)**2+
     &                       (abs(y0)-csy*dy0)**2+
     &                       (abs(z0)-csz*dz0)**2)

           end do
        end do

        return
        end

c-----------------------------------------------------------------------
c More agressive averaging than smooth_ah_r() above
c
c eps not used yet
c 
c-----------------------------------------------------------------------
        subroutine smooth_ah_r_b(AH_R,AH_w1,AH_eps,AH_Nchi,AH_Nphi)
        implicit none
        integer AH_Nchi,AH_Nphi
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_eps
        real*8 AH_w1(AH_Nchi,AH_Nphi)
        
        integer i,j,i1,j1,i0,j0
        integer avg_rad
        real*8 maxr,sum,r0,r1,num,r
        parameter (avg_rad=2)

        ! initialize fixed-size variables
        data i,j,i1,j1,i0,j0/0,0,0,0,0,0/

        data maxr,sum,r0,r1,num,r/0.0,0.0,0.0,0.0,0.0,0.0/

        !--------------------------------------------------------------

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              AH_w1(i,j)=AH_R(i,j)
           end do
        end do

        do i=1,AH_Nchi
           do j=1,AH_Nphi
              maxr=0
              sum=0
              num=0
              do i1=i-avg_rad,i+avg_rad
                 do j1=j-avg_rad,j+avg_rad
                    r=((i1-i)**2+(j1-j)**2)**(0.5d0)
                    if (r.lt.(avg_rad+0.01)) then
                       i0=i1
                       j0=j1
                       if (i0.lt.1) i0=2-i0
                       if (i0.gt.AH_Nchi) i0=AH_Nchi-(i0-AH_Nchi)
                       if (j0.lt.1) j0=j0+AH_Nphi-1
                       if (j0.gt.AH_Nphi) j0=j0-AH_Nphi+1
                       sum=sum+(avg_rad+1-r)*AH_w1(i0,j0)
                       num=num+(avg_rad+1-r)
                       if (abs(AH_w1(i0,j0)).gt.abs(maxr)) 
     &                    maxr=AH_w1(i0,j0)
                    end if
                 end do
              end do
              sum=sum/num
              !sum=(sum-maxr)/(num-1)
              AH_R(i,j)=sum
           end do
        end do

        r0=0
        r1=0
        do i=1,AH_Nphi
           r0=r0+AH_R(1,i)
           r1=r1+AH_R(AH_Nchi,i)
        end do
        r0=r0/AH_Nphi
        r1=r1/AH_Nphi
        do i=1,AH_Nphi
           AH_R(1,i)=r0
           AH_R(AH_Nchi,i)=r1
        end do

        do i=1,AH_Nchi
           AH_R(i,1)=AH_R(i,AH_Nphi)
        end do

        return 
        end

