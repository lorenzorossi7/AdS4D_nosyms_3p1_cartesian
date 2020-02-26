//============================================================================
// in cartesian coordinates t,x,y,z for x,y,z in [-1,1]
// using r=2*rho/(1-rho^2) compactification for rho=sqrt(x^2+y^2)
//
// application interface functions for AdS4D
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <m_util_r8.h>
#include <bbhutil.h>
#include <mpi.h>
#include "AdS4D.h"
#include "apph.h"

//=============================================================================
// if axisym=1, then 2+1 simulation in (x,y) plane 
//=============================================================================
int axisym=1;

//=============================================================================
// set in fparam for now
//=============================================================================
real AdS_L;

//=============================================================================
// Carsten's "constraint-damping" parameters
//=============================================================================
real kappa_cd,rho_cd;

//=============================================================================
// id (and other) parameters
//=============================================================================

// gaussians
real phi1_amp_1,phi1_B_1,phi1_C_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[3],phi1_ecc_1[3];
real phi1_amp_2,phi1_B_2,phi1_C_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[3],phi1_ecc_2[3];

// if > 0, initialize with exact BH
real ief_bh_r0;

//gauge parameters:
int gauge_t;
int gauge_i;
real c1_t,c2_t,c3_t;
real c1_i,c2_i,c3_i;
real rho1_t,rho2_t,rho3_t,rho4_t,xi1_t,xi2_t,cbulk_t;
real rho1_i,rho2_i,rho3_i,rho4_i,xi1_i,xi2_i,cbulk_i;
real rhoa,rhob;
real rhoc,rhod;

int cp_version; 

// excision parameters x-ex_xc  ex_r
real ex_rbuf[MAX_BHS];
int ex_reset_rbuf;
int ex_max_repop,ex_repop_buf,ex_repop_io;

// "internal" excision parameters, set by AH finder (eventually)
real ex_r[MAX_BHS][3],ex_xc[MAX_BHS][3];

int background,skip_constraints;
int output_ires,output_bdyquantities,output_AdS_mass,output_relkretschcentregrid,output_kretsch,output_relkretsch;
int reduced_ascii,reduction_factor;
int alltimes_ascii,timestep_ascii;

// new parameters in rtfile
int interptype,i_shift,regtype,stype;

int harmonize;

//extra dissipation
real diss_eps_k,diss_eps_y_cutoff;
int diss_kmax,diss_eps_k_cutoff_n,diss_bdy_k,diss_all_past_k,diss_all;

// excision parameters (more)
real ex_rbuf_a[MAX_BHS];

// AH parameters
int AH_Nchi[MAX_BHS],AH_Nphi[MAX_BHS],AH_Lmin[MAX_BHS],AH_Lmax[MAX_BHS],AH_find_best_fit[MAX_BHS];
int AH_max_iter[MAX_BHS],AH_freq[MAX_BHS],AH_freq_aft[MAX_BHS],AH_rsteps[MAX_BHS],AH_maxinc[MAX_BHS];
real AH_tol[MAX_BHS],AH_tol_aft[MAX_BHS],AH_r0[MAX_BHS],AH_lambda[MAX_BHS],AH_lambda_min[MAX_BHS];
real AH_eps[MAX_BHS],AH_r1[MAX_BHS],AH_tol_scale[MAX_BHS],AH_reset_scale[MAX_BHS];
real AH_xc[MAX_BHS][3],AH_max_tol_inc[MAX_BHS],AH_tmin[MAX_BHS],AH_omt_scale[MAX_BHS];
int use_AH_new_smooth,use_AH_new;
int c_AH;

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *cl_res;

real *phi1,*phi1_n,*phi1_np1,*phi1_nm1; // MGH, AMRH n/np1/nm1

real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1; 
real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1; 
real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1; 
real *gb_tz,*gb_tz_n,*gb_tz_np1,*gb_tz_nm1;
real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1; 
real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1; 
real *gb_xz,*gb_xz_n,*gb_xz_np1,*gb_xz_nm1;
real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1; 
real *gb_yz,*gb_yz_n,*gb_yz_np1,*gb_yz_nm1;
real *psi,*psi_n,*psi_np1,*psi_nm1; 

real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
real *Hb_y,*Hb_y_n,*Hb_y_np1,*Hb_y_nm1;
real *Hb_z,*Hb_z_n,*Hb_z_np1,*Hb_z_nm1;

real *phi1_t,*phi1_t_n;
real *gb_tt_t,*gb_tt_t_n;
real *gb_tx_t,*gb_tx_t_n;
real *gb_ty_t,*gb_ty_t_n;
real *gb_tz_t,*gb_tz_t_n;
real *gb_xx_t,*gb_xx_t_n;
real *gb_xy_t,*gb_xy_t_n;
real *gb_xz_t,*gb_xz_t_n;
real *gb_yy_t,*gb_yy_t_n;
real *gb_yz_t,*gb_yz_t_n;
real *psi_t,*psi_t_n;
real *Hb_t_t,*Hb_t_t_n;
real *Hb_x_t,*Hb_x_t_n;
real *Hb_y_t,*Hb_y_t_n;
real *Hb_z_t,*Hb_z_t_n;


real *w1,*mg_w1;
real *w2,*mg_w2;
real *w3,*mg_w3;
real *w4,*mg_w4;

real *mask,*mask_mg,*chr,*chr_mg;
real *chrbdy;
real *kg_ires,*alpha,*ricci,*theta,*f,*K;

real *phi1_res,*gb_res;
real *efe_all_ires;
real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
real *efe_tz_ires;
real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires,*efe_psi_ires;
real *efe_xz_ires,*efe_yz_ires;
int numbdypoints;
int basenumbdypoints;
int basebdy_Nchi,basebdy_Nxi;
int *vecbdypoints, *dsplsbdypoints;
int uniSize;
real *quasiset_tt,*quasiset_tchi,*quasiset_txi;
real *quasiset_chichi,*quasiset_chixi,*quasiset_xixi;
real *quasiset_massdensity, *quasiset_trace, *AdS_mass;
real *locoeffphi1;
real *xextrap,*yextrap,*zextrap;
real *kretsch;
real *relkretsch;
real *relkretschcentregrid, *lrelkretschcentregrid0, *maxrelkretschcentregrid0,*minrelkretschcentregrid0,*relkretschcentregrid0;
real *relkretsch, *lrelkretsch0, *maxrelkretsch0,*minrelkretsch0,*relkretsch0;
real *kretsch, *lkretsch0, *maxkretsch0,*minkretsch0,*kretsch0;

real *tfunction,*test1,*test2,*test3,*test4;
real *iresall,*irestt,*irestx,*iresty,*iresxx,*iresxy,*iresyy,*irespsi;
real *irestz;
real *iresxz;
real *iresyz;
real *qstt,*qstx,*qsty,*qsxx,*qsxy,*qsyy,*qspsi;
real *qsmass;
real *qsone;

real *zeta,*zeta_res,*zeta_lop,*zeta_rhs;

real *hb_t_res,*hb_i_res;
real *Hb_t_0,*Hb_x_0,*Hb_y_0,*Hb_z_0;

real *g_norms;

real *x,*y,*z;
int shape[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size,g_rank;
real base_bbox[6],bbox[6],dx,dy,dz,dt,dx_Lc,dy_Lc,dz_Lc;
int g_L;

int cl_res_gfn;

int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 

int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn; 
int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn; 
int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn; 
int gb_tz_gfn,gb_tz_n_gfn,gb_tz_np1_gfn,gb_tz_nm1_gfn;
int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn; 
int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn; 
int gb_xz_gfn,gb_xz_n_gfn,gb_xz_np1_gfn,gb_xz_nm1_gfn;
int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn; 
int gb_yz_gfn,gb_yz_n_gfn,gb_yz_np1_gfn,gb_yz_nm1_gfn;
int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn; 

int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;
int Hb_z_gfn,Hb_z_n_gfn,Hb_z_np1_gfn,Hb_z_nm1_gfn;

int phi1_t_gfn,phi1_t_n_gfn;
int gb_tt_t_gfn,gb_tt_t_n_gfn;
int gb_tx_t_gfn,gb_tx_t_n_gfn;
int gb_ty_t_gfn,gb_ty_t_n_gfn;
int gb_tz_t_gfn,gb_tz_t_n_gfn;
int gb_xx_t_gfn,gb_xx_t_n_gfn;
int gb_xy_t_gfn,gb_xy_t_n_gfn;
int gb_xz_t_gfn,gb_xz_t_n_gfn;
int gb_yy_t_gfn,gb_yy_t_n_gfn;
int gb_yz_t_gfn,gb_yz_t_n_gfn;
int psi_t_gfn,psi_t_n_gfn;
int Hb_t_t_gfn,Hb_t_t_n_gfn;
int Hb_x_t_gfn,Hb_x_t_n_gfn;
int Hb_y_t_gfn,Hb_y_t_n_gfn;
int Hb_z_t_gfn,Hb_z_t_n_gfn;


int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn, chrbdy_gfn;
real *gu_tt,*gu_tx,*gu_ty,*gu_tz,*gu_xx,*gu_xy,*gu_xz,*gu_yy,*gu_yz,*gu_psi,*m_g_det;
int kg_ires_gfn,alpha_gfn,theta_gfn,f_gfn;

int phi1_res_gfn,gb_res_gfn;
int efe_all_ires_gfn;
int efe_tt_ires_gfn,efe_tx_ires_gfn,efe_ty_ires_gfn;
int efe_tz_ires_gfn;
int efe_xx_ires_gfn,efe_xy_ires_gfn,efe_yy_ires_gfn,efe_psi_ires_gfn;
int efe_xz_ires_gfn,efe_yz_ires_gfn;
int quasiset_tt_gfn,quasiset_tchi_gfn,quasiset_txi_gfn;
int quasiset_chichi_gfn,quasiset_chixi_gfn,quasiset_xixi_gfn;
int quasiset_massdensity_gfn,quasiset_trace_gfn,AdS_mass_gfn;
int locoeffphi1_gfn;
int xextrap_gfn,yextrap_gfn,zextrap_gfn;
int kretsch_gfn;
int relkretsch_gfn;
int relkretschcentregrid_gfn;

//for saving kretschmann related quantities as .sdf files
char name[256];

int kretsch_rank;
int kretsch_shape[3];
char kretsch_cnames[256];
double *kretsch_coords;

int baseNx,baseNy,baseNz;

int tfunction_gfn,test1_gfn,test2_gfn,test3_gfn,test4_gfn;
int iresall_gfn,irestt_gfn,irestx_gfn,iresty_gfn,irestz_gfn,iresxx_gfn,iresxy_gfn,iresxz_gfn,iresyy_gfn,iresyz_gfn,irespsi_gfn;
int qstt_gfn,qstx_gfn,qsty_gfn,qsxx_gfn,qsxy_gfn,qsyy_gfn,qspsi_gfn;
int qsmass_gfn;
int qsone_gfn;

int zeta_gfn,zeta_res_gfn,zeta_lop_gfn,zeta_rhs_gfn;

int hb_t_res_gfn,hb_i_res_gfn;
int Hb_t_0_gfn,Hb_x_0_gfn,Hb_y_0_gfn,Hb_z_0_gfn;

int w1_gfn,mg_w1_gfn;
int w2_gfn,mg_w2_gfn;
int w3_gfn,mg_w3_gfn;
int w4_gfn,mg_w4_gfn;

int skip_ires=0;
int skip_exp=0;

int gu_tt_gfn,gu_tx_gfn,gu_ty_gfn,gu_tz_gfn,gu_xx_gfn;
int gu_xy_gfn,gu_xz_gfn,gu_yy_gfn,gu_yz_gfn,gu_psi_gfn,m_g_det_gfn;

//=============================================================================
// arrays holding AH shape and other hypersurface info
//=============================================================================
real *AH_theta[MAX_BHS],*AH_R[MAX_BHS],*AH_w1[MAX_BHS],*AH_w2[MAX_BHS],*AH_w3[MAX_BHS];
real *AH_theta_ads[MAX_BHS];
int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

//=============================================================================
// arrays holding quasiset of CFT on ESU at outer bdy
//=============================================================================
real *lquasiset_tt0,*lquasiset_tchi0,*lquasiset_txi0;
real *lquasiset_chichi0,*lquasiset_chixi0,*lquasiset_xixi0;
real *lquasiset_massdensity0,*lquasiset_trace0,*lAdS_mass0;
real *llocoeffphi10;
real *maxquasiset_tt0,*maxquasiset_tchi0,*maxquasiset_txi0;
real *maxquasiset_chichi0,*maxquasiset_chixi0,*maxquasiset_xixi0;
real *maxquasiset_massdensity0,*maxquasiset_trace0;
real *maxlocoeffphi10;
real *minquasiset_tt0,*minquasiset_tchi0,*minquasiset_txi0;
real *minquasiset_chichi0,*minquasiset_chixi0,*minquasiset_xixi0;
real *minquasiset_massdensity0,*minquasiset_trace0;
real *minlocoeffphi10;
real *quasiset_tt0,*quasiset_tchi0,*quasiset_txi0;
real *quasiset_chichi0,*quasiset_chixi0,*quasiset_xixi0;
real *quasiset_massdensity0,*quasiset_trace0,*AdS_mass0;
real *locoeffphi10;

real *xextrap0,*yextrap0,*zextrap0;
real *rhoextrap0,*chiextrap0,*xiextrap0;
real *rhobdy0,*chibdy0,*xibdy0;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((cl_res_gfn   = PAMR_get_gfn("cl_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_gfn     = PAMR_get_gfn("phi1",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_nm1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_n_gfn   = PAMR_get_gfn("phi1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_np1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_gfn     = PAMR_get_gfn("gb_tt",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_nm1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_n_gfn   = PAMR_get_gfn("gb_tt",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_np1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tx_gfn     = PAMR_get_gfn("gb_tx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_nm1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_n_gfn   = PAMR_get_gfn("gb_tx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_np1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_ty_gfn     = PAMR_get_gfn("gb_ty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_nm1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_n_gfn   = PAMR_get_gfn("gb_ty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_np1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tz_gfn     = PAMR_get_gfn("gb_tz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tz_nm1_gfn = PAMR_get_gfn("gb_tz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tz_n_gfn   = PAMR_get_gfn("gb_tz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tz_np1_gfn = PAMR_get_gfn("gb_tz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xx_gfn     = PAMR_get_gfn("gb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_nm1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_n_gfn   = PAMR_get_gfn("gb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_np1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xy_gfn     = PAMR_get_gfn("gb_xy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_nm1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_n_gfn   = PAMR_get_gfn("gb_xy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_np1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xz_gfn     = PAMR_get_gfn("gb_xz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xz_nm1_gfn = PAMR_get_gfn("gb_xz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xz_n_gfn   = PAMR_get_gfn("gb_xz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xz_np1_gfn = PAMR_get_gfn("gb_xz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_yy_gfn     = PAMR_get_gfn("gb_yy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_nm1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_n_gfn   = PAMR_get_gfn("gb_yy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_np1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_yz_gfn     = PAMR_get_gfn("gb_yz",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yz_nm1_gfn = PAMR_get_gfn("gb_yz",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yz_n_gfn   = PAMR_get_gfn("gb_yz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yz_np1_gfn = PAMR_get_gfn("gb_yz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((psi_gfn     = PAMR_get_gfn("psi",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_nm1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_n_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_np1_gfn = PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_gfn      = PAMR_get_gfn("Hb_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_nm1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_n_gfn    = PAMR_get_gfn("Hb_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_np1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_t_gfn = PAMR_get_gfn("phi1_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_t_n_gfn = PAMR_get_gfn("phi1_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_t_gfn = PAMR_get_gfn("gb_tt_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_n_gfn = PAMR_get_gfn("gb_tt_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tx_t_gfn = PAMR_get_gfn("gb_tx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_n_gfn = PAMR_get_gfn("gb_tx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_ty_t_gfn = PAMR_get_gfn("gb_ty_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_t_n_gfn = PAMR_get_gfn("gb_ty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tz_t_gfn = PAMR_get_gfn("gb_tz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tz_t_n_gfn = PAMR_get_gfn("gb_tz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xx_t_gfn = PAMR_get_gfn("gb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_n_gfn = PAMR_get_gfn("gb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xy_t_gfn = PAMR_get_gfn("gb_xy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_t_n_gfn = PAMR_get_gfn("gb_xy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_xz_t_gfn = PAMR_get_gfn("gb_xz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xz_t_n_gfn = PAMR_get_gfn("gb_xz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_yy_t_gfn = PAMR_get_gfn("gb_yy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_t_n_gfn = PAMR_get_gfn("gb_yy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_yz_t_gfn = PAMR_get_gfn("gb_yz_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yz_t_n_gfn = PAMR_get_gfn("gb_yz_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((psi_t_gfn = PAMR_get_gfn("psi_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_n_gfn = PAMR_get_gfn("psi_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_t_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_n_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_n_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_n_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_t_gfn  = PAMR_get_gfn("Hb_z_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_t_n_gfn  = PAMR_get_gfn("Hb_z_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_x_gfn      = PAMR_get_gfn("Hb_x",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_nm1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_n_gfn    = PAMR_get_gfn("Hb_x",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_np1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_y_gfn      = PAMR_get_gfn("Hb_y",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_nm1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_n_gfn    = PAMR_get_gfn("Hb_y",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_np1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_z_gfn      = PAMR_get_gfn("Hb_z",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_nm1_gfn  = PAMR_get_gfn("Hb_z",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_n_gfn    = PAMR_get_gfn("Hb_z",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_np1_gfn  = PAMR_get_gfn("Hb_z",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);


    if ((zeta_gfn     = PAMR_get_gfn("zeta",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_res_gfn = PAMR_get_gfn("zeta_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_lop_gfn = PAMR_get_gfn("zeta_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_rhs_gfn = PAMR_get_gfn("zeta_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((chrbdy_gfn     = PAMR_get_gfn("chrbdy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi1_res_gfn  = PAMR_get_gfn("phi1_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_res_gfn    = PAMR_get_gfn("gb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_all_ires_gfn    = PAMR_get_gfn("efe_all_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_tt_ires_gfn    = PAMR_get_gfn("efe_tt_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_tx_ires_gfn    = PAMR_get_gfn("efe_tx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_ty_ires_gfn    = PAMR_get_gfn("efe_ty_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_tz_ires_gfn    = PAMR_get_gfn("efe_tz_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_xx_ires_gfn    = PAMR_get_gfn("efe_xx_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_xy_ires_gfn    = PAMR_get_gfn("efe_xy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_xz_ires_gfn    = PAMR_get_gfn("efe_xz_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_yy_ires_gfn    = PAMR_get_gfn("efe_yy_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_yz_ires_gfn    = PAMR_get_gfn("efe_yz_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_psi_ires_gfn    = PAMR_get_gfn("efe_psi_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tt_gfn    = PAMR_get_gfn("quasiset_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tchi_gfn    = PAMR_get_gfn("quasiset_tchi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_txi_gfn    = PAMR_get_gfn("quasiset_txi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_chichi_gfn    = PAMR_get_gfn("quasiset_chichi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_chixi_gfn    = PAMR_get_gfn("quasiset_chixi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_xixi_gfn    = PAMR_get_gfn("quasiset_xixi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_massdensity_gfn    = PAMR_get_gfn("quasiset_massdensity",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_trace_gfn    = PAMR_get_gfn("quasiset_trace",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((locoeffphi1_gfn    = PAMR_get_gfn("locoeffphi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((AdS_mass_gfn    = PAMR_get_gfn("AdS_mass",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((xextrap_gfn    = PAMR_get_gfn("xextrap",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((yextrap_gfn    = PAMR_get_gfn("yextrap",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((zextrap_gfn    = PAMR_get_gfn("zextrap",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((kretsch_gfn    = PAMR_get_gfn("kretsch",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((relkretsch_gfn    = PAMR_get_gfn("relkretsch",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((relkretschcentregrid_gfn    = PAMR_get_gfn("relkretschcentregrid",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_t_res_gfn  = PAMR_get_gfn("hb_t_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_i_res_gfn  = PAMR_get_gfn("hb_i_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_0_gfn  = PAMR_get_gfn("Hb_t_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_0_gfn  = PAMR_get_gfn("Hb_x_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_0_gfn  = PAMR_get_gfn("Hb_y_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_z_0_gfn  = PAMR_get_gfn("Hb_z_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w1_gfn   = PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w2_gfn   = PAMR_get_gfn("w2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w3_gfn   = PAMR_get_gfn("w3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w4_gfn   = PAMR_get_gfn("w4",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((kg_ires_gfn= PAMR_get_gfn("kg_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((alpha_gfn= PAMR_get_gfn("alpha",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((theta_gfn= PAMR_get_gfn("theta",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((f_gfn= PAMR_get_gfn("f",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((tfunction_gfn  = PAMR_get_gfn("tfunction",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((test1_gfn  = PAMR_get_gfn("test1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test2_gfn  = PAMR_get_gfn("test2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test3_gfn  = PAMR_get_gfn("test3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((test4_gfn  = PAMR_get_gfn("test4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((iresall_gfn  = PAMR_get_gfn("iresall",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irestt_gfn  = PAMR_get_gfn("irestt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irestx_gfn  = PAMR_get_gfn("irestx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresty_gfn  = PAMR_get_gfn("iresty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irestz_gfn  = PAMR_get_gfn("irestz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresxx_gfn  = PAMR_get_gfn("iresxx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresxy_gfn  = PAMR_get_gfn("iresxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresxz_gfn  = PAMR_get_gfn("iresxz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresyy_gfn  = PAMR_get_gfn("iresyy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresyz_gfn  = PAMR_get_gfn("iresyz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((irespsi_gfn  = PAMR_get_gfn("irespsi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qstt_gfn  = PAMR_get_gfn("qstt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qstx_gfn  = PAMR_get_gfn("qstx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsty_gfn  = PAMR_get_gfn("qsty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsxx_gfn  = PAMR_get_gfn("qsxx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsxy_gfn  = PAMR_get_gfn("qsxy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsyy_gfn  = PAMR_get_gfn("qsyy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qspsi_gfn  = PAMR_get_gfn("qspsi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsmass_gfn  = PAMR_get_gfn("qsmass",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((qsone_gfn  = PAMR_get_gfn("qsone",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w1_gfn   = PAMR_get_gfn("mg_w1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w2_gfn   = PAMR_get_gfn("mg_w2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w3_gfn   = PAMR_get_gfn("mg_w3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w4_gfn   = PAMR_get_gfn("mg_w4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((gu_tt_gfn   = PAMR_get_gfn("gu_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_tx_gfn   = PAMR_get_gfn("gu_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_ty_gfn   = PAMR_get_gfn("gu_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_tz_gfn   = PAMR_get_gfn("gu_tz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_xx_gfn   = PAMR_get_gfn("gu_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_xy_gfn   = PAMR_get_gfn("gu_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_xz_gfn   = PAMR_get_gfn("gu_xz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_yy_gfn   = PAMR_get_gfn("gu_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_yz_gfn   = PAMR_get_gfn("gu_yz",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gu_psi_gfn   = PAMR_get_gfn("gu_psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((m_g_det_gfn = PAMR_get_gfn("m_g_det",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    g_norms=AMRD_get_global_norms();
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr_bbox(void)
{
   real dx0[3];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
      if (PAMR_get_max_lev(PAMR_AMRH)>1) PAMR_get_dxdt(2,dx0,&dt); else PAMR_get_dxdt(1,dx0,&dt);
      dx_Lc=dx0[0];
      dy_Lc=dx0[1];
      dz_Lc=dx0[2];
   }

   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   dx=dx0[0];
   dy=dx0[1];
   dz=dx0[2];

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   if ((bbox[4]-base_bbox[4])<dz/2) phys_bdy[4]=1; else phys_bdy[4]=0;
   if ((base_bbox[5]-bbox[5])<dz/2) phys_bdy[5]=1; else phys_bdy[5]=0;

   Nx=shape[0];
   Ny=shape[1];
   Nz=shape[2];

   size=Nx*Ny*Nz;
}

void ldptr(void)
{
   real *x0[3],*gfs[PAMR_MAX_GFNS];

   ldptr_bbox();

   PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1]; 
   z=x0[2];

   PAMR_get_g_gfs(gfs);

   cl_res   = gfs[cl_res_gfn-1];

   phi1     = gfs[phi1_gfn-1];
   phi1_n   = gfs[phi1_n_gfn-1];
   phi1_np1 = gfs[phi1_np1_gfn-1];
   phi1_nm1 = gfs[phi1_nm1_gfn-1];

   gb_tt     = gfs[gb_tt_gfn-1];
   gb_tt_n   = gfs[gb_tt_n_gfn-1];
   gb_tt_np1 = gfs[gb_tt_np1_gfn-1];
   gb_tt_nm1 = gfs[gb_tt_nm1_gfn-1];

   gb_tx     = gfs[gb_tx_gfn-1];
   gb_tx_n   = gfs[gb_tx_n_gfn-1];
   gb_tx_np1 = gfs[gb_tx_np1_gfn-1];
   gb_tx_nm1 = gfs[gb_tx_nm1_gfn-1];

   gb_ty     = gfs[gb_ty_gfn-1];
   gb_ty_n   = gfs[gb_ty_n_gfn-1];
   gb_ty_np1 = gfs[gb_ty_np1_gfn-1];
   gb_ty_nm1 = gfs[gb_ty_nm1_gfn-1];

   gb_tz     = gfs[gb_tz_gfn-1];
   gb_tz_n   = gfs[gb_tz_n_gfn-1];
   gb_tz_np1 = gfs[gb_tz_np1_gfn-1];
   gb_tz_nm1 = gfs[gb_tz_nm1_gfn-1];

   gb_xx     = gfs[gb_xx_gfn-1];
   gb_xx_n   = gfs[gb_xx_n_gfn-1];
   gb_xx_np1 = gfs[gb_xx_np1_gfn-1];
   gb_xx_nm1 = gfs[gb_xx_nm1_gfn-1];

   gb_xy     = gfs[gb_xy_gfn-1];
   gb_xy_n   = gfs[gb_xy_n_gfn-1];
   gb_xy_np1 = gfs[gb_xy_np1_gfn-1];
   gb_xy_nm1 = gfs[gb_xy_nm1_gfn-1];

   gb_xz     = gfs[gb_xz_gfn-1];
   gb_xz_n   = gfs[gb_xz_n_gfn-1];
   gb_xz_np1 = gfs[gb_xz_np1_gfn-1];
   gb_xz_nm1 = gfs[gb_xz_nm1_gfn-1];

   gb_yy     = gfs[gb_yy_gfn-1];
   gb_yy_n   = gfs[gb_yy_n_gfn-1];
   gb_yy_np1 = gfs[gb_yy_np1_gfn-1];
   gb_yy_nm1 = gfs[gb_yy_nm1_gfn-1];

   gb_yz     = gfs[gb_yz_gfn-1];
   gb_yz_n   = gfs[gb_yz_n_gfn-1];
   gb_yz_np1 = gfs[gb_yz_np1_gfn-1];
   gb_yz_nm1 = gfs[gb_yz_nm1_gfn-1];

   psi     = gfs[psi_gfn-1];
   psi_n   = gfs[psi_n_gfn-1];
   psi_np1 = gfs[psi_np1_gfn-1];
   psi_nm1 = gfs[psi_nm1_gfn-1];

   Hb_t      = gfs[Hb_t_gfn-1];
   Hb_t_n    = gfs[Hb_t_n_gfn-1];
   Hb_t_nm1  = gfs[Hb_t_nm1_gfn-1];
   Hb_t_np1  = gfs[Hb_t_np1_gfn-1];

   Hb_x      = gfs[Hb_x_gfn-1];
   Hb_x_n    = gfs[Hb_x_n_gfn-1];
   Hb_x_nm1  = gfs[Hb_x_nm1_gfn-1];
   Hb_x_np1  = gfs[Hb_x_np1_gfn-1];

   Hb_y      = gfs[Hb_y_gfn-1];
   Hb_y_n    = gfs[Hb_y_n_gfn-1];
   Hb_y_nm1  = gfs[Hb_y_nm1_gfn-1];
   Hb_y_np1  = gfs[Hb_y_np1_gfn-1];

   Hb_z      = gfs[Hb_z_gfn-1];
   Hb_z_n    = gfs[Hb_z_n_gfn-1];
   Hb_z_nm1  = gfs[Hb_z_nm1_gfn-1];
   Hb_z_np1  = gfs[Hb_z_np1_gfn-1];

   phi1_t = gfs[phi1_t_gfn-1];
   phi1_t_n = gfs[phi1_t_n_gfn-1];

   gb_tt_t = gfs[gb_tt_t_gfn-1];
   gb_tt_t_n = gfs[gb_tt_t_n_gfn-1];

   gb_tx_t = gfs[gb_tx_t_gfn-1];
   gb_tx_t_n = gfs[gb_tx_t_n_gfn-1];

   gb_ty_t = gfs[gb_ty_t_gfn-1];
   gb_ty_t_n = gfs[gb_ty_t_n_gfn-1];

   gb_tz_t = gfs[gb_tz_t_gfn-1];
   gb_tz_t_n = gfs[gb_tz_t_n_gfn-1];

   gb_xx_t = gfs[gb_xx_t_gfn-1];
   gb_xx_t_n = gfs[gb_xx_t_n_gfn-1];

   gb_xy_t = gfs[gb_xy_t_gfn-1];
   gb_xy_t_n = gfs[gb_xy_t_n_gfn-1];

   gb_xz_t = gfs[gb_xz_t_gfn-1];
   gb_xz_t_n = gfs[gb_xz_t_n_gfn-1];

   gb_yy_t = gfs[gb_yy_t_gfn-1];
   gb_yy_t_n = gfs[gb_yy_t_n_gfn-1];

   gb_yz_t = gfs[gb_yz_t_gfn-1];
   gb_yz_t_n = gfs[gb_yz_t_n_gfn-1];

   psi_t = gfs[psi_t_gfn-1];
   psi_t_n = gfs[psi_t_n_gfn-1];

   Hb_t_t  = gfs[Hb_t_t_gfn-1];
   Hb_t_t_n  = gfs[Hb_t_t_n_gfn-1];
   Hb_x_t  = gfs[Hb_x_t_gfn-1];
   Hb_x_t_n  = gfs[Hb_x_t_n_gfn-1];
   Hb_y_t  = gfs[Hb_y_t_gfn-1];
   Hb_y_t_n  = gfs[Hb_y_t_n_gfn-1];
   Hb_z_t  = gfs[Hb_z_t_gfn-1];
   Hb_z_t_n  = gfs[Hb_z_t_n_gfn-1];

   zeta     = gfs[zeta_gfn-1];
   zeta_lop = gfs[zeta_lop_gfn-1];
   zeta_res = gfs[zeta_res_gfn-1];
   zeta_rhs = gfs[zeta_rhs_gfn-1];

   mask    = gfs[mask_gfn-1];
   mask_mg = gfs[mask_mg_gfn-1];
   chr = gfs[chr_gfn-1]; 
   chr_mg = gfs[chr_mg_gfn-1]; 
   chrbdy = gfs[chrbdy_gfn-1];

   phi1_res  = gfs[phi1_res_gfn-1];
   gb_res    = gfs[gb_res_gfn-1];
   efe_all_ires  = gfs[efe_all_ires_gfn-1];
   efe_tt_ires  = gfs[efe_tt_ires_gfn-1];
   efe_tx_ires  = gfs[efe_tx_ires_gfn-1];
   efe_ty_ires  = gfs[efe_ty_ires_gfn-1];
   efe_tz_ires  = gfs[efe_tz_ires_gfn-1];
   efe_xx_ires  = gfs[efe_xx_ires_gfn-1];
   efe_xy_ires  = gfs[efe_xy_ires_gfn-1];
   efe_xz_ires  = gfs[efe_xz_ires_gfn-1];
   efe_yy_ires  = gfs[efe_yy_ires_gfn-1];
   efe_yz_ires  = gfs[efe_yz_ires_gfn-1];
   efe_psi_ires  = gfs[efe_psi_ires_gfn-1];
   quasiset_tt  = gfs[quasiset_tt_gfn-1];
   quasiset_tchi  = gfs[quasiset_tchi_gfn-1];
   quasiset_txi  = gfs[quasiset_txi_gfn-1];
   quasiset_chichi  = gfs[quasiset_chichi_gfn-1];
   quasiset_chixi  = gfs[quasiset_chixi_gfn-1];
   quasiset_xixi  = gfs[quasiset_xixi_gfn-1];
   quasiset_massdensity  = gfs[quasiset_massdensity_gfn-1];
   quasiset_trace  = gfs[quasiset_trace_gfn-1];
   locoeffphi1  = gfs[locoeffphi1_gfn-1];
   AdS_mass  = gfs[AdS_mass_gfn-1];
   xextrap  = gfs[xextrap_gfn-1];
   yextrap  = gfs[yextrap_gfn-1];
   zextrap  = gfs[zextrap_gfn-1];
   kretsch  = gfs[kretsch_gfn-1];
   relkretsch  = gfs[relkretsch_gfn-1];
   relkretschcentregrid = gfs[relkretschcentregrid_gfn-1];
   hb_t_res  = gfs[hb_t_res_gfn-1];
   hb_i_res  = gfs[hb_i_res_gfn-1];
   Hb_t_0  = gfs[Hb_t_0_gfn-1];
   Hb_x_0  = gfs[Hb_x_0_gfn-1];
   Hb_y_0  = gfs[Hb_y_0_gfn-1];
   Hb_z_0  = gfs[Hb_z_0_gfn-1];
   w1 = gfs[w1_gfn-1];
   w2 = gfs[w2_gfn-1];
   w3 = gfs[w3_gfn-1];
   w4 = gfs[w4_gfn-1];
   kg_ires = gfs[kg_ires_gfn-1];
   alpha = gfs[alpha_gfn-1];
   theta = gfs[theta_gfn-1];
   f = gfs[f_gfn-1];

   tfunction = gfs[tfunction_gfn-1];
   test1 = gfs[test1_gfn-1];
   test2 = gfs[test2_gfn-1];
   test3 = gfs[test3_gfn-1];
   test4 = gfs[test4_gfn-1];
   iresall = gfs[iresall_gfn-1];
   irestt = gfs[irestt_gfn-1];
   irestx = gfs[irestx_gfn-1];
   iresty = gfs[iresty_gfn-1];
   irestz = gfs[irestz_gfn-1];
   iresxx = gfs[iresxx_gfn-1];
   iresxy = gfs[iresxy_gfn-1];
   iresxz = gfs[iresxz_gfn-1];
   iresyy = gfs[iresyy_gfn-1];
   iresyz = gfs[iresyz_gfn-1];
   irespsi = gfs[irespsi_gfn-1];
   qstt = gfs[qstt_gfn-1];
   qstx = gfs[qstx_gfn-1];
   qsty = gfs[qsty_gfn-1];
   qsxx = gfs[qsxx_gfn-1];
   qsxy = gfs[qsxy_gfn-1];
   qsyy = gfs[qsyy_gfn-1];
   qspsi = gfs[qspsi_gfn-1];
   qsmass = gfs[qsmass_gfn-1];
   qsone = gfs[qsone_gfn-1];

   mg_w1    =gfs[mg_w1_gfn-1]; 
   mg_w2    =gfs[mg_w2_gfn-1]; 
   mg_w3    =gfs[mg_w3_gfn-1]; 
   mg_w4    =gfs[mg_w4_gfn-1]; 

   gu_tt = gfs[gu_tt_gfn-1];
   gu_tx = gfs[gu_tx_gfn-1];
   gu_ty = gfs[gu_ty_gfn-1];
   gu_tz = gfs[gu_tz_gfn-1];
   gu_xx = gfs[gu_xx_gfn-1];
   gu_xy = gfs[gu_xy_gfn-1];
   gu_xz = gfs[gu_xz_gfn-1];
   gu_yy = gfs[gu_yy_gfn-1];
   gu_yz = gfs[gu_yz_gfn-1];
   gu_psi = gfs[gu_psi_gfn-1];
   m_g_det = gfs[m_g_det_gfn-1];

}

//=============================================================================
// PAMR_get_dxdt() only works with AMR hierarchy levels ... here we use
// lambda for dt, but this only works if rhosp=rhotm
//=============================================================================
void ldptr_mg(void)
{
   real lambda;

   ldptr();

   dx=x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0];
   PAMR_get_lambda(&lambda);
   dt=lambda*dx;
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny*Nz; i++) f[i]=c;
}

void zero_f(real *f)
{
   const_f(f,0);
}

void zero_f_ex(real *f, real *chr)
{
   int i;

   for (i=0; i<Nx*Ny*Nz; i++) if (chr[i]==AMRD_ex) f[i]=0;
}

real norm_l2(real *f, real *cmask, real *chr)
{
   int i;
   real norm=0;
   int sum=0;

   for (i=0; i<Nx*Ny*Nz; i++) 
      if (cmask[i]==AMRD_CMASK_ON && (chr[i]!=AMRD_ex)) { sum++; norm+=f[i]*f[i]; }

   if (!sum) sum=1;
   return (sqrt(norm/sum));
}

//=============================================================================
// the following zeros the AH_max_iter and ex_r (semiaxes of best-fit ellipses) 
// of engulfed AHs
//=============================================================================
void remove_redundant_AH()
{
   int i1,i2,is_1in2,is_2in1,intr;
   real d,min_bhr,bhr;
   int rem[MAX_BHS];

   for (i1=0; i1<MAX_BHS; i1++) rem[i1]=0;

   for (i1=0; i1<MAX_BHS; i1++)
   {
      for (i2=i1+1; i2<MAX_BHS; i2++)
      {
         if (ex_r[i1][0]>0 && ex_r[i2][0]>0)
         {
            is_inside_(&is_1in2,&intr,ex_r[i1],ex_xc[i1],ex_r[i2],ex_xc[i2],&AMRD_dim);
            is_inside_(&is_2in1,&intr,ex_r[i2],ex_xc[i2],ex_r[i1],ex_xc[i1],&AMRD_dim);

//(unless this is disabled, AH[4] is removed at first time step)
            if (is_1in2 || is_2in1)
            {
               if (my_rank==0) printf("\n\nremove_redundant_AH: the following two excised regions "
                  "overlap ... removing the smallest one (is_1in2=%i, is_2in1=%i).\n"
                  "1: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf],  "
                  "2: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n\n",is_1in2,is_2in1,
                  ex_r[i1][0],ex_r[i1][1],ex_r[i1][2],ex_xc[i1][0],ex_xc[i1][1],ex_xc[i1][2],
                  ex_r[i2][0],ex_r[i2][1],ex_r[i2][2],ex_xc[i2][0],ex_xc[i2][1],ex_xc[i2][2]);
               if (is_2in1) rem[i2]=1; else rem[i1]=1;
            }

            // special case
            if (i2==MAX_BHS-1 && rem[i2]==1)
            { 
               if (my_rank==0) printf("\nremove_redundant_AH: 'unremoving' the tentative encapsulating BH\n");
               rem[i2]=0; 
            }
            if (i2==MAX_BHS-1 && i1!=MAX_BHS-2 && intr) 
            {
               if (my_rank==0) printf("\nremove_redundant_AH: encapsulating BH found, removing intersecting BHs.\n"
                  "1: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf],  "
                  "2: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n\n",
                  ex_r[i1][0],ex_r[i1][1],ex_r[i1][2],ex_xc[i1][0],ex_xc[i1][1],ex_xc[i1][2],
                  ex_r[i2][0],ex_r[i2][1],ex_r[i2][2],ex_xc[i2][0],ex_xc[i2][1],ex_xc[i2][2]);
               rem[i1]=1;
            }

            if (is_1in2 && is_2in1 && my_rank==0) printf("WARNING ... both is_1in2 && is_2in1\n");
         }
      }
   }

   for (i1=0; i1<MAX_BHS; i1++) if (rem[i1])
   {
       ex_r[i1][0]=0;
       ex_r[i1][1]=0;
       ex_r[i1][2]=0;
       AH_max_iter[i1]=0;
   }
}

//=============================================================================
// the following returns 1 if the test ellipsoid does not intersect any
// existing horizons, *and* is entirely contained within the finest overlapping
// level
//=============================================================================
#define MAX_GRIDS 20
#define EX_MIN_PTS_IN 4
int no_AH_intersect(real ex_r0[3],real ex_xc0[3],int inew)
{
   int i,i1;
   real d;
   real sgh[6*MAX_GRIDS],dx0[3],dt0,bound[6],max_xy;
   int num,L,Lf,contained,ex0_in_exi,ex0_int_exi,ltrace=1;
   int exi_in_ex0,exi_int_ex0;

   for (i=0; i<MAX_BHS; i++)
   {
      if (ex_r[i][0]>0 && i!=inew)
      {
         is_inside_(&ex0_in_exi,&ex0_int_exi,ex_r0,ex_xc0,ex_r[i],ex_xc[i],&AMRD_dim);
         is_inside_(&exi_in_ex0,&exi_int_ex0,ex_r[i],ex_xc[i],ex_r0,ex_xc0,&AMRD_dim);

         if (ex0_int_exi)
         {
            // special case
            if (inew==MAX_BHS-1)
            {
              if (my_rank==0) printf("\nno_AH_intersect: encapsulating BH found, ignoring intersections\n");
              return 1;
            }
            
            if (my_rank==0) printf("\nno_AH_intersect: new BH intersects existing BH\n"
                   "new: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf.%lf],"
                   "current(%i): ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n",
                   ex_r0[0],ex_r0[1],ex_r0[2],ex_xc0[0],ex_xc0[1],ex_xc0[2],i,
                   ex_r[i][0],ex_r[i][1],ex_r[i][2],ex_xc[i][0],ex_xc[i][1],ex_xc[i][2]);
            if (my_rank==0) printf("ex0_in_exi=%i,ex0_int_exi=%i,exi_in_ex0=%i,exi_int_ex0=%i\n\n",ex0_in_exi,ex0_int_exi,exi_in_ex0,exi_int_ex0);
         }

         if (ex0_in_exi)
         {
            if (my_rank==0) printf("no_AH_intersect: WARNING ... new BH is *entirely* contained in existing BH\n"
                   "new: ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf.%lf],"
                   "current(%i): ex_r=[%lf,%lf,%lf], ex_xc=[%lf,%lf,%lf]\n",
                   ex_r0[0],ex_r0[1],ex_r0[2],ex_xc0[0],ex_xc0[1],ex_xc0[2],i,
                   ex_r[i][0],ex_r[i][1],ex_r[i][2],ex_xc[i][0],ex_xc[i][1],ex_xc[i][2]);
            if (my_rank==0) printf("ex0_in_exi=%i,ex0_int_exi=%i,exi_in_ex0=%i,exi_int_ex0=%i\n\n",ex0_in_exi,ex0_int_exi,exi_in_ex0,exi_int_ex0);
            return 0;
         }

         if (ex0_int_exi && !exi_in_ex0) return 0;
      }
   }

   return 1;
}


//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int AdS4D_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void AdS4D_var_pre_init(char *pfile)
{
   AMRD_echo_params=1;
   AMRD_int_param(pfile,"echo_params",&AMRD_echo_params,1);

   cp_version=ADS5D_CP_VERSION;
   AMRD_int_param(pfile,"cp_version",&cp_version,1);
   
   AdS_L=1.0; AMRD_real_param(pfile,"AdS_L",&AdS_L,1);

   regtype=1; AMRD_int_param(pfile,"regtype",&regtype,1);
   stype=1; AMRD_int_param(pfile,"stype",&stype,1);
   interptype=2; AMRD_int_param(pfile,"interptype",&interptype,1);
   i_shift=0; AMRD_int_param(pfile,"i_shift",&i_shift,1);

   return;
}

void AdS4D_var_post_init(char *pfile)
{
   int i,j,k,ind,l;
   char buf[64];
   real rmin,deltar;

   int valid;


   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading AdS4D parameters:\n\n");
      fflush(stdout);
   }

   phi1_amp_1=phi1_B_1=phi1_C_1=phi1_r0_1=phi1_x0_1[0]=phi1_x0_1[1]=phi1_x0_1[2]=phi1_ecc_1[0]=phi1_ecc_1[1]=phi1_ecc_1[2]=0;
   phi1_amp_2=phi1_B_2=phi1_C_2=phi1_r0_1=phi1_x0_2[0]=phi1_x0_2[1]=phi1_x0_2[2]=phi1_ecc_2[0]=phi1_ecc_2[1]=phi1_ecc_2[2]=0;

   AMRD_real_param(pfile,"phi1_amp_1",&phi1_amp_1,1);
   AMRD_real_param(pfile,"phi1_B_1",&phi1_B_1,1);
   AMRD_real_param(pfile,"phi1_C_1",&phi1_C_1,1);
   AMRD_real_param(pfile,"phi1_r0_1",&phi1_r0_1,1);
   AMRD_real_param(pfile,"phi1_delta_1",&phi1_delta_1,1);
   AMRD_real_param(pfile,"phi1_x0_1",phi1_x0_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_ecc_1",phi1_ecc_1,AMRD_dim);

   AMRD_real_param(pfile,"phi1_amp_2",&phi1_amp_2,1);
   AMRD_real_param(pfile,"phi1_B_2",&phi1_B_2,1);
   AMRD_real_param(pfile,"phi1_C_2",&phi1_C_2,1);
   AMRD_real_param(pfile,"phi1_r0_2",&phi1_r0_2,1);
   AMRD_real_param(pfile,"phi1_delta_2",&phi1_delta_2,1);
   AMRD_real_param(pfile,"phi1_x0_2",phi1_x0_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_ecc_2",phi1_ecc_2,AMRD_dim);

   kappa_cd=0; AMRD_real_param(pfile,"kappa_cd",&kappa_cd,1);
   rho_cd=0; AMRD_real_param(pfile,"rho_cd",&rho_cd,1);

   diss_eps_k=0; AMRD_real_param(pfile,"diss_eps_k",&diss_eps_k,1);
   diss_eps_y_cutoff=1; AMRD_real_param(pfile,"diss_eps_y_cutoff",&diss_eps_y_cutoff,1);
   diss_kmax=0; AMRD_int_param(pfile,"diss_kmax",&diss_kmax,1);
   diss_bdy_k=0; AMRD_int_param(pfile,"diss_bdy_k",&diss_bdy_k,1);
   diss_all_past_k=0; AMRD_int_param(pfile,"diss_all_past_k",&diss_all_past_k,1);
   diss_eps_k_cutoff_n=0; AMRD_int_param(pfile,"diss_eps_k_cutoff_n",&diss_eps_k_cutoff_n,1);
   diss_all=1; AMRD_int_param(pfile,"diss_all",&diss_all,1);

   background=0; AMRD_int_param(pfile,"background",&background,1);
   skip_constraints=0; AMRD_int_param(pfile,"skip_constraints",&skip_constraints,1);
   output_ires=0; AMRD_int_param(pfile,"output_ires",&output_ires,1);
   output_bdyquantities=0; AMRD_int_param(pfile,"output_bdyquantities",&output_bdyquantities,1);
   output_AdS_mass=0; AMRD_int_param(pfile,"output_AdS_mass",&output_AdS_mass,1);
   reduced_ascii=0; AMRD_int_param(pfile,"reduced_ascii",&reduced_ascii,1);
   reduction_factor=1; AMRD_int_param(pfile,"reduction_factor",&reduction_factor,1);
   alltimes_ascii=0; AMRD_int_param(pfile,"alltimes_ascii",&alltimes_ascii,1);
   timestep_ascii=0; AMRD_int_param(pfile,"timestep_ascii",&timestep_ascii,1);

   harmonize=0; AMRD_int_param(pfile,"harmonize",&harmonize,1);

   gauge_t=0; AMRD_int_param(pfile,"gauge_t",&gauge_t,1);
   c1_t=0; AMRD_real_param(pfile,"c1_t",&c1_t,1);
   c2_t=0; AMRD_real_param(pfile,"c2_t",&c2_t,1);
   c3_t=0; AMRD_real_param(pfile,"c3_t",&c3_t,1);
   rho1_t=0; AMRD_real_param(pfile,"rho1_t",&rho1_t,1);
   rho2_t=0; AMRD_real_param(pfile,"rho2_t",&rho2_t,1);
   rho3_t=1; AMRD_real_param(pfile,"rho3_t",&rho3_t,1);
   rho4_t=1; AMRD_real_param(pfile,"rho4_t",&rho4_t,1);
   xi1_t=0; AMRD_real_param(pfile,"xi1_t",&xi1_t,1);
   xi2_t=0; AMRD_real_param(pfile,"xi2_t",&xi2_t,1);
   cbulk_t=2; AMRD_real_param(pfile,"cbulk_t",&cbulk_t,1);

   gauge_i=0; AMRD_int_param(pfile,"gauge_i",&gauge_i,1);
   c1_i=0; AMRD_real_param(pfile,"c1_i",&c1_i,1);
   c2_i=0; AMRD_real_param(pfile,"c2_i",&c2_i,1);
   c3_i=0; AMRD_real_param(pfile,"c3_i",&c3_i,1);
   rho1_i=0; AMRD_real_param(pfile,"rho1_i",&rho1_i,1);
   rho2_i=0; AMRD_real_param(pfile,"rho2_i",&rho2_i,1);
   rho3_i=1; AMRD_real_param(pfile,"rho3_i",&rho3_i,1);
   rho4_i=1; AMRD_real_param(pfile,"rho4_i",&rho4_i,1);
   xi1_i=0; AMRD_real_param(pfile,"xi1_i",&xi1_i,1);
   xi2_i=0; AMRD_real_param(pfile,"xi2_i",&xi2_i,1);
   cbulk_i=2; AMRD_real_param(pfile,"cbulk_i",&cbulk_i,1);

   rhoa=1; AMRD_real_param(pfile,"rhoa",&rhoa,1);
   rhob=1; AMRD_real_param(pfile,"rhob",&rhob,1);

   rhoc=1; AMRD_real_param(pfile,"rhoc",&rhoc,1);
   rhod=1; AMRD_real_param(pfile,"rhod",&rhod,1);

   ex_reset_rbuf=0; AMRD_int_param(pfile,"ex_reset_rbuf",&ex_reset_rbuf,1);
   output_relkretschcentregrid=0; AMRD_int_param(pfile,"output_relkretschcentregrid",&output_relkretschcentregrid,1);
//   output_relkretsch=0; AMRD_int_param(pfile,"output_relkretsch",&output_relkretsch,1);
   output_kretsch=0; AMRD_int_param(pfile,"output_kretsch",&output_kretsch,1);

  //allocate memory for relative Kretschmann scalar at the centre of the grid
  if (output_relkretschcentregrid)
  {
      lrelkretschcentregrid0= malloc(sizeof(real));
      maxrelkretschcentregrid0= malloc(sizeof(real));
      minrelkretschcentregrid0= malloc(sizeof(real));
      relkretschcentregrid0= malloc(sizeof(real));
  }
  if (output_kretsch)
  {  
      lrelkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      maxrelkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      minrelkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      relkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));


      lkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      maxkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      minkretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));
      kretsch0= malloc(AMRD_base_shape[0]*AMRD_base_shape[1]*AMRD_base_shape[2]*sizeof(real));

  }

   // set fraction, 1-ex_rbuf, of AH radius to be excised
   for (j=0; j<MAX_BHS; j++)
   {
      if (!AMRD_cp_restart || ex_reset_rbuf || ex_r[j][0]<=0)
      {
         if (j==0) { if (!AMRD_cp_restart) ex_rbuf[j]=0; sprintf(buf,"ex_rbuf"); }
         else { if (!AMRD_cp_restart) ex_rbuf[j]=ex_rbuf[0]; sprintf(buf,"ex_rbuf_%i",j+1); }
         AMRD_real_param(pfile,buf,&ex_rbuf[j],1);
         if (ex_rbuf[j]<0 || ex_rbuf[j]>1 ) printf("WARNING ... ex_rbuf[%i]=%lf is outside of standard range\n",j,ex_rbuf[j]);
      }

      if (!AMRD_cp_restart)
      {
         ex_xc[j][0]=ex_xc[j][1]=ex_xc[j][2]=0;
         ex_r[j][0]=ex_r[j][1]=ex_r[j][2]=0;
      }
   }

   // set AH parameters
   int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
   int found_count_AH[MAX_BHS];

   for (l=0; l<MAX_BHS; l++)
   {
      // if restoring from check-point files, then this just searches for AH from scratch
      if (AMRD_cp_restart) found_AH[l]=0;

      // because the AH shape is saved, we can't currently change that upon a restart (unless we haven't yet found one)
      if (!AMRD_cp_restart || !found_AH[l])
      {
         if (AMRD_cp_restart) free(AH_R[l]);
         if (l==0) { AH_Nchi[l]=17; sprintf(buf,"AH_Nchi"); }
         else { AH_Nchi[l]=AH_Nchi[0]; sprintf(buf,"AH_Nchi_%i",l+1); }
         AMRD_int_param(pfile,buf,&AH_Nchi[l],1);

         if (l==0) { AH_Nphi[l]=17; sprintf(buf,"AH_Nphi"); }
         else { AH_Nphi[l]=AH_Nphi[0]; sprintf(buf,"AH_Nphi_%i",l+1); }
         AMRD_int_param(pfile,buf,&AH_Nphi[l],1);

         if (!( ((AH_Nchi[l]-1)/2)==(AH_Nchi[l]/2) ) ||
             !( ((AH_Nphi[l]-1)/2)==(AH_Nphi[l]/2) ))
         {
            printf("\n\n\n\n\n WARNING: SMOOTHING AND REGULARIZATION IN AH ROUTINES ASSUME\n"
                   " AN ODD NUMBER OF POINTS IN AH_Nchi,AH_Nphi\n\n\n\n\n");
         }
      }

      if (l==0) { AH_Lmin[l]=2; sprintf(buf,"AH_Lmin"); }
      else { AH_Lmin[l]=AH_Lmin[0]; sprintf(buf,"AH_Lmin_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_Lmin[l],1);

      if (l==0) { AH_Lmax[l]=100; sprintf(buf,"AH_Lmax"); }
      else { AH_Lmax[l]=AH_Lmax[0]; sprintf(buf,"AH_Lmax_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_Lmax[l],1);

      if (l==0) { AH_max_iter[l]=0; sprintf(buf,"AH_max_iter"); }
      else { AH_max_iter[l]=AH_max_iter[0]; sprintf(buf,"AH_max_iter_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_max_iter[l],1);

      if (l==0) { AH_freq[l]=1; sprintf(buf,"AH_freq"); }
      else { AH_freq[l]=AH_freq[0]; sprintf(buf,"AH_freq_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_freq[l],1);

      if (l==0) { AH_freq_aft[l]=1; sprintf(buf,"AH_freq_aft"); }
      else { AH_freq_aft[l]=AH_freq_aft[0]; sprintf(buf,"AH_freq_aft_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_freq_aft[l],1);

      if (AMRD_cp_restart)
      {
         if (found_AH[l]) freq0[l]=AH_freq_aft[l];
         else freq0[l]=AH_freq[l];
      }

      if (l==0) { AH_rsteps[l]=1; sprintf(buf,"AH_rsteps"); }
      else { AH_rsteps[l]=AH_rsteps[0]; sprintf(buf,"AH_rsteps_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_rsteps[l],1);

      if (l==0) { AH_maxinc[l]=10; sprintf(buf,"AH_maxinc"); }
      else { AH_maxinc[l]=AH_maxinc[0]; sprintf(buf,"AH_maxinc_%i",l+1); }
      AMRD_int_param(pfile,buf,&AH_maxinc[l],1);

      if (l==0) { AH_tol[l]=1e-2; sprintf(buf,"AH_tol"); }
      else { AH_tol[l]=AH_tol[0]; sprintf(buf,"AH_tol_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_tol[l],1);

      AH_tol_aft[l]=AH_tol[l];
      if (l==0) sprintf(buf,"AH_tol_aft"); else sprintf(buf,"AH_tol_aft_%i",l+1);
      AMRD_real_param(pfile,buf,&AH_tol_aft[l],1);

      if (l==0) { AH_r0[l]=0.1; sprintf(buf,"AH_r0"); }
      else { AH_r0[l]=AH_r0[0]; sprintf(buf,"AH_r0_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_r0[l],1);

      if (l==0) { AH_r1[l]=0.2; sprintf(buf,"AH_r1"); }
      else { AH_r1[l]=AH_r1[0]; sprintf(buf,"AH_r1_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_r1[l],1);

      if (l==0) { AH_lambda[l]=0.1; sprintf(buf,"AH_lambda"); }
      else { AH_lambda[l]=AH_lambda[0]; sprintf(buf,"AH_lambda_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_lambda[l],1);

      if (l==0) { AH_lambda_min[l]=0.1; sprintf(buf,"AH_lambda_min"); }
      else { AH_lambda_min[l]=AH_lambda_min[0]; sprintf(buf,"AH_lambda_min_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_lambda_min[l],1);

      if (l==0) { AH_eps[l]=0.0; sprintf(buf,"AH_eps"); }
      else { AH_eps[l]=AH_eps[0]; sprintf(buf,"AH_eps_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_eps[l],1);

      if (l==0) { AH_tmin[l]=0.0; sprintf(buf,"AH_tmin"); }
      else { AH_tmin[l]=AH_tmin[0]; sprintf(buf,"AH_tmin_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_tmin[l],1);

      if (l==0) { AH_reset_scale[l]=0; sprintf(buf,"AH_reset_scale"); }
      else { AH_reset_scale[l]=AH_reset_scale[0]; sprintf(buf,"AH_reset_scale_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_reset_scale[l],1);

      if (l==0) { AH_omt_scale[l]=1.1; sprintf(buf,"AH_omt_scale"); }
      else { AH_omt_scale[l]=AH_omt_scale[0]; sprintf(buf,"AH_omt_scale_%i",l+1); }
      AMRD_real_param(pfile,buf,&AH_omt_scale[l],1);

      if (!AMRD_cp_restart || !found_AH[l])
      {
         AH_xc[l][0]=AH_xc[l][1]=0;
         if (l==0) sprintf(buf,"AH_xc");
         else sprintf(buf,"AH_xc_%i",l+1);
         AMRD_real_param(pfile,buf,AH_xc[l],AMRD_dim);
      }

      if (AH_rsteps[l]<1) AMRD_stop("error ... AH_rsteps<1\n","");

      AH_theta[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      if (!AMRD_cp_restart || !found_AH[l]) AH_R[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      AH_w1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      AH_w2[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      AH_w3[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      AH_theta_ads[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
      AH_own[l]=(int *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(int));
      AH_lev[l]=(int *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(int));

   }

   ex_max_repop=0; AMRD_int_param(pfile,"ex_max_repop",&ex_max_repop,1);
   ex_repop_buf=1; AMRD_int_param(pfile,"ex_repop_buf",&ex_repop_buf,1);
   ex_repop_io=1; AMRD_int_param(pfile,"ex_repop_io",&ex_repop_io,1);
   if (abs(ex_repop_io)<1 || abs(ex_repop_io)>4) AMRD_stop("invalid |ex_repop_io| ... must be 1,2,3 or 4","");

   ief_bh_r0=0; AMRD_real_param(pfile,"ief_bh_r0",&ief_bh_r0,1);


   // later have a proper routine to convert from AH to excision parameters.
   // rh is radius in uncompcatified coordicates, rhoh in compactified,
   // ief_bh_r0 is BH radius parameter, ex_r is excision radius
   int ah_finder_is_off=1; 
   for (l=0; l<MAX_BHS; l++) {if (AH_max_iter[l]!=0) ah_finder_is_off=0;}
    if (ief_bh_r0)
    { 
     real rh,rhoh,mh;
     rh=-pow(AdS_L,2)
         /(pow(3,(1.0/3.0)))
         /(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
         +(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
         /(pow(3,(2.0/3.0)));
     mh=ief_bh_r0/2;
     rhoh=(-1 + sqrt(1 + pow(rh,2)))/rh;
//     ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh;  //*(1-ex_rbuf[0]); //HB
     if (my_rank==0) 
     {
       printf("\nBH initial data\n"
                            "r0/L=%lf, rh/L=%lf, mass M = r0/2 = rh*(1+rh^2/L^2)/2 = %lf\n"
                            "Initial BH Schwarzschild radius=%lf, (%lf in compactified (code) coords)\n" 
                            "Excision buffer (i.e. size of the evolved region within the AH) ex_rbuf[0]=%lf\n\n"
                            ,ief_bh_r0/AdS_L,rh/AdS_L,mh,rh,rhoh,ex_rbuf[0]);
     }
  
       ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh;   
     if (ah_finder_is_off) 
     {
//       ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh;
       if (my_rank==0)   printf("\n ... AH finder is off so we excise, AT ALL TIME STEPS, points with compactified radius smaller than rhoh*(1-ex_rbuf[0])=%lf ... \n",rhoh*(1-ex_rbuf[0]));
     }
     else //if we start from bh initial data and AH finder is not off
     {
      //we change AH_r0[l] to avoid looking for AH in the excised region
      for (l=0; l<MAX_BHS; l++)
      { 
       if (rhoh*(1-ex_rbuf[0])>AH_r0[l]) AH_r0[l]=rhoh;  
       if (rhoh*(1-ex_rbuf[0])>AH_r1[l]) AH_r1[l]=rhoh;
      }
     }

    }
    else
    {
//     ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=(1-ex_rbuf[0]); //HB
//     if (my_rank==0) printf("\nscalar field initial data from Hamiltonian constraint solver, no BH\n"
//                            "We excise inside fixed excision radius=%lf\n\n",ex_r[0][0]);
     if (my_rank==0) printf("\nscalar field initial data from Hamiltonian constraint solver\n");
     if (ah_finder_is_off)
     {
       ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=1;
       if (my_rank==0)   printf("\n ... AH finder is off so we excise, AT ALL TIME STEPS, points with compactified radius smaller than (1-ex_rbuf[0])=%lf ... \n",(1-ex_rbuf[0]));
     }
    }

   if (AMRD_do_ex==0) AMRD_stop("require excision to be on","");

   PAMR_excision_on("chr",&AdS4D_fill_ex_mask,AMRD_ex,1);


   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void AdS4D_AMRH_var_clear(void)
{
   int i,j,k,ind;

   ldptr();

   zero_f(phi1_n); 

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//         ind=i+Nx*(j+Ny*k);
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         printf("AdS4D_AMRH_var_clear-PRE init_ghb_ads\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

   init_ghb_ads_(gb_tt_n,gb_tx_n,gb_ty_n,
                 gb_tz_n,
                 gb_xx_n,gb_xy_n,
                 gb_xz_n,
                 gb_yy_n,
                 gb_yz_n,
                 psi_n,
                 Hb_t_n,Hb_x_n,Hb_y_n,
                 Hb_z_n,
                 &AdS_L,x,y,z,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        { 
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_ghb_ads\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2) ... following vars also used in 
// t0_cnst_data
//=============================================================================
void AdS4D_free_data(void)
{
   int i;

   ldptr();

   AdS4D_AMRH_var_clear(); // constrained variables are set post-MG

   zero_f(phi1_t_n); // holds initial time derivatives for ID

   gauss3d_(phi1_n,&phi1_amp_1,&phi1_B_1,&phi1_C_1,&phi1_r0_1,&phi1_delta_1,&phi1_x0_1[0],&phi1_x0_1[1],&phi1_x0_1[2],
            &phi1_ecc_1[0],&phi1_ecc_1[1],&phi1_ecc_1[2],&AdS_L,x,y,z,&Nx,&Ny,&Nz,&rhoc,&rhod,&stype);

   gauss3d_(w1,&phi1_amp_2,&phi1_B_2,&phi1_C_2,&phi1_r0_2,&phi1_delta_2,&phi1_x0_2[0],&phi1_x0_2[1],&phi1_x0_2[2],
            &phi1_ecc_2[0],&phi1_ecc_2[1],&phi1_ecc_2[2],&AdS_L,x,y,z,&Nx,&Ny,&Nz,&rhoc,&rhod,&stype);

   for (i=0; i<size; i++) phi1_n[i]+=w1[i]; 

   return;
}  

//=============================================================================
// Initialize any "elliptic_vars_t0" post construction of MGH, but before
// the start of vcycling.
//=============================================================================
void AdS4D_elliptic_vars_t0_init(void)
{
   // initializes dt, dx, dy
   ldptr_mg();

   // initializes zeta conformal factor to background AdS value zeta=1
   const_f(zeta,1);

}

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//
// Here we also initialize past time level information if 
// AMRD_id_pl_method==3
//
// NOTE: not cleaning up memory allocation after reading in square black hole 
//       data
//
// NOTE: np1,n,nm1 variables are allocated only at the top level of the MG hierarchy,
//       so do an if(f_nm1){...}, for example, to make sure we are at the top level
//=============================================================================
void AdS4D_t0_cnst_data(void)
{
   int i,j,k,ind;
   real ct,rho;

   ldptr_mg();

//   printf("AdS4D_t0_cnst_data is called");

   // initialize time derivatives of gbars,hbars
   if (gb_xx_nm1)
   {

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-PRE init_ghbdot_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//         printf("gb_xx_nm1[ind]=%lf,gb_xx_n[ind]=%lf,gb_xx_np1[ind]=%lf\n",gb_xx_nm1[ind],gb_xx_n[ind],gb_xx_np1[ind]);
//        }
//       }
//      }
//   }


     init_ghbdot_(gb_tt_n,gb_tx_n,gb_ty_n,
                  gb_tz_n,
                  gb_xx_n,gb_xy_n,
                  gb_xz_n,
                  gb_yy_n,
                  gb_yz_n,
                  psi_n,gb_tt_t_n,gb_tx_t_n,gb_ty_t_n,
                  gb_tz_t_n,
                  gb_xx_t_n,gb_xy_t_n,
                  gb_xz_t_n,
                  gb_yy_t_n,
                  gb_yz_t_n,
                  psi_t_n,Hb_t_n,Hb_x_n,Hb_y_n,
                  Hb_z_n,
                  Hb_t_t_n,Hb_x_t_n,Hb_y_t_n,
                  Hb_z_t_n,
                  &AdS_L,phys_bdy,x,y,z,&dt,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_ghbdot_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }


   }


   // initialize gbars and nm1, np1 time levels
   if ((background || skip_constraints) && ief_bh_r0==0)
   {

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-PRE init_ghb_ads_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

     init_ghb_ads_(gb_tt,gb_tx,gb_ty,
                   gb_tz,
                   gb_xx,gb_xy,
                   gb_xz,
                   gb_yy,
                   gb_yz,
                   psi,Hb_t,Hb_x,Hb_y,
                   Hb_z,
                   &AdS_L,x,y,z,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_ghb_ads_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }


     for (i=0; i<size; i++)
     {
       gb_tt_np1[i]=gb_tt_nm1[i]=gb_tt[i];
       gb_tx_np1[i]=gb_tx_nm1[i]=gb_tx[i];
       gb_ty_np1[i]=gb_ty_nm1[i]=gb_ty[i];
       gb_tz_np1[i]=gb_tz_nm1[i]=gb_tz[i];
       gb_xx_np1[i]=gb_xx_nm1[i]=gb_xx[i];
       gb_xy_np1[i]=gb_xy_nm1[i]=gb_xy[i];
       gb_xz_np1[i]=gb_xz_nm1[i]=gb_xz[i];
       gb_yy_np1[i]=gb_yy_nm1[i]=gb_yy[i];
       gb_yz_np1[i]=gb_yz_nm1[i]=gb_yz[i];
       psi_np1[i]=psi_nm1[i]=psi[i];
       Hb_t_np1[i]=Hb_t_nm1[i]=Hb_t_n[i];
       Hb_x_np1[i]=Hb_x_nm1[i]=Hb_x_n[i];
       Hb_y_np1[i]=Hb_y_nm1[i]=Hb_y_n[i];
       Hb_z_np1[i]=Hb_z_nm1[i]=Hb_z_n[i];
       phi1_np1[i]=phi1_nm1[i]=phi1[i];
     }
   }
   else if (background || skip_constraints)
   {
     if (gb_xx_nm1) //"np1,n,nm1" variables only allocated on finest MG level
     {

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-PRE init_ads4d_bh_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

       init_ads4d_bh_(&ief_bh_r0,&AdS_L,gb_tt,gb_tx,gb_ty,
                      gb_tz,
                      gb_xx,gb_xy,
                      gb_xz,
                      gb_yy,
                      gb_yz,
                      psi,gb_tt_t_n,gb_tx_t_n,gb_ty_t_n,
                      gb_tz_t_n,
                      gb_xx_t_n,gb_xy_t_n,
                      gb_xz_t_n,
                      gb_yy_t_n,
                      gb_yz_t_n,
                      psi_t_n,Hb_t,Hb_x,Hb_y,
                      Hb_z,
                      Hb_t_t_n,Hb_x_t_n,Hb_y_t_n,
                      Hb_z_t_n,
                      phys_bdy,x,y,z,&dt,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_ads4d_bh_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

       for (i=0; i<size; i++)
       {
         gb_tt_np1[i]=gb_tt_nm1[i]=gb_tt[i];
         gb_tx_np1[i]=gb_tx_nm1[i]=gb_tx[i];
         gb_ty_np1[i]=gb_ty_nm1[i]=gb_ty[i];
         gb_tz_np1[i]=gb_tz_nm1[i]=gb_tz[i];
         gb_xx_np1[i]=gb_xx_nm1[i]=gb_xx[i];
         gb_xy_np1[i]=gb_xy_nm1[i]=gb_xy[i];
         gb_xz_np1[i]=gb_xz_nm1[i]=gb_xz[i];
         gb_yy_np1[i]=gb_yy_nm1[i]=gb_yy[i];
         gb_yz_np1[i]=gb_yz_nm1[i]=gb_yz[i];
         psi_np1[i]=psi_nm1[i]=psi[i];
         Hb_t_np1[i]=Hb_t_nm1[i]=Hb_t_n[i];
         Hb_x_np1[i]=Hb_x_nm1[i]=Hb_x_n[i];
         Hb_y_np1[i]=Hb_y_nm1[i]=Hb_y_n[i];
         Hb_z_np1[i]=Hb_z_nm1[i]=Hb_z_n[i];
         phi1_np1[i]=phi1_nm1[i]=phi1[i];
       }
     }
   }
   else
   {

    if (gb_xx_nm1) //"np1,n,nm1" variables only allocated on finest MG level //IT WAS NOT HERE IN PREVIOUS VERSIONS, ASK HANS
    {

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-PRE init_ghb_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//         printf("gb_xx_nm1[ind]=%lf,gb_xx_n[ind]=%lf,gb_xx_np1[ind]=%lf\n",gb_xx_nm1[ind],gb_xx_n[ind],gb_xx_np1[ind]);
//        }
//       }
//      }
//   }

     init_ghb_(zeta,
               gb_tt,gb_tx,gb_ty,
               gb_tz,
               gb_xx,gb_xy,
               gb_xz,
               gb_yy,
               gb_yz,
               psi,Hb_t,Hb_x,Hb_y,
               Hb_z,
               &AdS_L,mask_mg,phys_bdy,x,y,z,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype,&rhoa,&rhob);

//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_ghb_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//         printf("gb_xx_nm1[ind]=%lf,gb_xx_n[ind]=%lf,gb_xx_np1[ind]=%lf\n",gb_xx_nm1[ind],gb_xx_n[ind],gb_xx_np1[ind]);
//        }
//       }
//      }
//   }




    }

   }

   // initialize hbars 
   if (AMRD_id_pl_method==3 && gb_xx_nm1)
   {

//   for (i=0; i<Nx; i++)
//   {  
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-PRE init_hb_ AND init_nm1_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

//      for (i=0; i<Nx; i++)
//      {
//         for (j=0; j<Ny; j++)
//         {
//          for (k=0; k<Nz; k++)
//          {
//            ind=i+Nx*(j+Ny*k);
//            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
//            if ((chr[ind]==AMRD_ex)&&(rho<0.7))
//            printf("PT EXCISED: i=%i,j=%i,k=%i,x=%lf,y=%lf,z=%lf,rho=%lf\n",i,j,k,x[i],y[j],z[k],rho);
//          }
//         }
//       }

     init_hb_(gb_tt_np1,gb_tt_n,gb_tt_nm1,
              gb_tx_np1,gb_tx_n,gb_tx_nm1,
              gb_ty_np1,gb_ty_n,gb_ty_nm1,
              gb_tz_np1,gb_tz_n,gb_tz_nm1,
              gb_xx_np1,gb_xx_n,gb_xx_nm1,
              gb_xy_np1,gb_xy_n,gb_xy_nm1,
              gb_xz_np1,gb_xz_n,gb_xz_nm1,
              gb_yy_np1,gb_yy_n,gb_yy_nm1,
              gb_yz_np1,gb_yz_n,gb_yz_nm1,
              psi_np1,psi_n,psi_nm1,
              Hb_t_n,Hb_x_n,Hb_y_n,
              Hb_z_n,
              &AdS_L,phys_bdy,x,y,z,&dt,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

       init_nm1_(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
                 gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
                 gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
                 gb_tz_np1,gb_tz_n,gb_tz_nm1,gb_tz_t_n,
                 gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
                 gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
                 gb_xz_np1,gb_xz_n,gb_xz_nm1,gb_xz_t_n,
                 gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
                 gb_yz_np1,gb_yz_n,gb_yz_nm1,gb_yz_t_n,
                 psi_np1,psi_n,psi_nm1,psi_t_n,
                 Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
                 Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
                 Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t_n,
                 Hb_z_np1,Hb_z_n,Hb_z_nm1,Hb_z_t_n,
                 phi1_np1,phi1_n,phi1_nm1,phi1_t_n,tfunction,
                 &AdS_L,phys_bdy,x,y,z,&dt,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);


//   for (i=0; i<Nx; i++)
//   {  
//      for (j=0; j<Ny; j++)
//      {
//       for (k=0; k<Nz; k++)
//       {
//        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
//        {
//         ind=i+Nx*(j+Ny*k);
//         printf("AdS4D_AMRH_var_clear-POST init_hb_ AND init_nm1_\n");
//         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"
//                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));
//         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);
//        }
//       }
//      }
//   }

   }

   // store initial source functions 
   if (gb_xx_nm1)
   {
     for (i=0; i<size; i++)
     {
       Hb_t_0[i]=Hb_t[i];
       Hb_x_0[i]=Hb_x[i];
       Hb_y_0[i]=Hb_y[i];
       Hb_z_0[i]=Hb_z[i];
     }
   }

   // fill test functions with conformal factor diagnostic profiles
   for (i=0; i<Nx; i++)
   {
      for (j=0; j<Ny; j++)
      {
       for (k=0; k<Nz; k++)
       {
         ind=i+Nx*(j+Ny*k);
         test1[ind]=zeta[ind];
         test2[ind]=zeta_lop[ind];
         test3[ind]=zeta_rhs[ind];
         test4[ind]=zeta_res[ind];
       }
      } 
   }
  
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void AdS4D_pre_io_calc(void)
{
   ldptr();

   int i,j,k,ind;
   int is,ie,js,je,ks,ke;
   int i0,j0,k0;
   int j_shift;
   int n;
   int Lf,Lc;
   real ct,rho;

   dx=x[1]-x[0]; dy=y[1]-y[0]; dz=z[1]-z[0];
   ct=PAMR_get_time(g_L);

   if (regtype==6) j_shift=3;
   if (regtype==5 || regtype==4 || regtype==3) j_shift=2;
   if (regtype==2 || regtype==1) j_shift=1;

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   int ivecNt=AMRD_steps/AMRD_save_ivec0[3]+1; //+1 to include t=0
   int lsteps=AMRD_lsteps[Lc-1];

//   printf("AdS4D_pre_io_calc is called,lsteps=%i\n",lsteps);
//           fflush(stdout);


         MPI_Comm_size(MPI_COMM_WORLD,&uniSize);

//   // output independent residual
//   if (output_ires)
//   {
//    if (gb_xx_nm1)
//    {
      // compute independent residuals of the AdS4D system
    if (ct!=0)
    {
      //(NOTE: for t>t0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
      // so here, time level n is the most advanced time level)

     // output independent residual
     if ((output_ires)||(output_kretsch)||(output_relkretschcentregrid))
     { 

      ires_(efe_all_ires,
         efe_tt_ires,efe_tx_ires,efe_ty_ires,
         efe_tz_ires,
         efe_xx_ires,efe_xy_ires,
         efe_xz_ires,
         efe_yy_ires,
         efe_yz_ires,
         efe_psi_ires,
         kretsch,
         relkretsch,
         relkretschcentregrid,
         gb_tt_n,gb_tt_nm1,gb_tt_np1,
         gb_tx_n,gb_tx_nm1,gb_tx_np1,
         gb_ty_n,gb_ty_nm1,gb_ty_np1,
         gb_tz_n,gb_tz_nm1,gb_tz_np1,
         gb_xx_n,gb_xx_nm1,gb_xx_np1,
         gb_xy_n,gb_xy_nm1,gb_xy_np1,
         gb_xz_n,gb_xz_nm1,gb_xz_np1,
         gb_yy_n,gb_yy_nm1,gb_yy_np1,
         gb_yz_n,gb_yz_nm1,gb_yz_np1,
         psi_n,psi_nm1,psi_np1,
         Hb_t_n,Hb_t_nm1,Hb_t_np1,
         Hb_x_n,Hb_x_nm1,Hb_x_np1,
         Hb_y_n,Hb_y_nm1,Hb_y_np1,
         Hb_z_n,Hb_z_nm1,Hb_z_np1,
         phi1_n,phi1_nm1,phi1_np1,
         x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

     } 
    }
    else
    {
      //(NOTE: for t=t0, have *not* cycled time sequence, so still np1,n,nm1,
      // so here, time level np1 is the most advanced time level)
     if ((output_ires)||(output_kretsch)||(output_relkretschcentregrid))
     { 

      ires_(efe_all_ires,
         efe_tt_ires,efe_tx_ires,efe_ty_ires,
         efe_tz_ires,
         efe_xx_ires,efe_xy_ires,
         efe_xz_ires,
         efe_yy_ires,
         efe_yz_ires,
         efe_psi_ires,
         kretsch,
         relkretsch,
         relkretschcentregrid,
         gb_tt_np1,gb_tt_n,gb_tt_nm1,
         gb_tx_np1,gb_tx_n,gb_tx_nm1,
         gb_ty_np1,gb_ty_n,gb_ty_nm1,
         gb_tz_np1,gb_tz_n,gb_tz_nm1,
         gb_xx_np1,gb_xx_n,gb_xx_nm1,
         gb_xy_np1,gb_xy_n,gb_xy_nm1,
         gb_xz_np1,gb_xz_n,gb_xz_nm1,
         gb_yy_np1,gb_yy_n,gb_yy_nm1,
         gb_yz_np1,gb_yz_n,gb_yz_nm1,
         psi_np1,psi_n,psi_nm1,
         Hb_t_np1,Hb_t_n,Hb_t_nm1,
         Hb_x_np1,Hb_x_n,Hb_x_nm1,
         Hb_y_np1,Hb_y_n,Hb_y_nm1,
         Hb_z_np1,Hb_z_n,Hb_z_nm1,
         phi1_np1,phi1_n,phi1_nm1,
         x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

//      for (i=0; i<Nx; i++)
//      {
//         for (j=0; j<Ny; j++)
//         {
//          for (k=0; k<Nz; k++)
//          {
//            ind=i+Nx*(j+Ny*k);
//            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
//            if (rho<1) printf("PRE_IO_CALC:my_rank=%i,ct=%lf,x=%lf,y=%lf,z=%lf,rho=%lf,relkretsch[ind]=%lf\n",my_rank,ct,x[i],y[j],z[k],rho,relkretsch[ind]);
//          }
//         }
//      }


      }

         if (output_relkretschcentregrid)
         {
          *lrelkretschcentregrid0= *relkretschcentregrid;
         }

         if (output_kretsch)
         {

          is=(bbox[0]-base_bbox[0])/dx;    // for "left-most"  processors, includes i=0
          ie=(bbox[1]-base_bbox[0])/dx+1;  // for "right-most" processors, includes i=Nx
          js=(bbox[2]-base_bbox[2])/dy;    // for "left-most"  processors, includes j=0
          je=(bbox[3]-base_bbox[2])/dy+1;  // for "right-most" processors, includes j=Ny
          ks=(bbox[4]-base_bbox[4])/dz;    // for "left-most"  processors, includes k=0
          ke=(bbox[5]-base_bbox[4])/dz+1;  // for "right-most" processors, includes k=Ny

           baseNx=AMRD_base_shape[0];
           baseNy=AMRD_base_shape[1];
           baseNz=AMRD_base_shape[2];

           for (i=is; i<ie; i++)
           {
            for (j=js; j<je; j++)
            {
             for (k=ks; k<ke; k++)
             {
              lrelkretsch0[i+baseNx*(j+baseNy*k)]=relkretsch[(i-is)+Nx*((j-js)+Ny*(k-ks))];
              lkretsch0[i+baseNx*(j+baseNy*k)]=kretsch[(i-is)+Nx*((j-js)+Ny*(k-ks))];
             }
            }
           }
         }


        if (lsteps==0) 
        {
         if (output_relkretschcentregrid)
         {
          MPI_Allreduce((lrelkretschcentregrid0),(maxrelkretschcentregrid0),1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
          MPI_Allreduce((lrelkretschcentregrid0),(minrelkretschcentregrid0),1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

          if (uniSize>1)
          {
           *relkretschcentregrid0=*maxrelkretschcentregrid0+*minrelkretschcentregrid0;
          }
          else
          {
           *relkretschcentregrid0=*maxrelkretschcentregrid0;
          }
           if (my_rank==0)
           {
            // save relkretschcentregrid as ascii
            FILE * fp;
            fp = fopen ("ascii_t_relkretschcentregrid.txt", "a+");
                fprintf(fp,"%24.16e %24.16e \n",ct,*relkretschcentregrid0);
            fclose(fp);
           }
          }

       if (output_kretsch)
       {


           kretsch_rank=3;
           kretsch_shape[0]=baseNx;
           kretsch_shape[1]=baseNy;
           kretsch_shape[2]=baseNz;
           sprintf(kretsch_cnames,"x|y|z");
           kretsch_coords = malloc(sizeof(double)*(baseNx+baseNy+baseNz));
           for (i=0;i<baseNx;i++) {kretsch_coords[i] = base_bbox[0]+i*dx;}
           for (j=0;j<baseNy;j++) {kretsch_coords[j+baseNx] = base_bbox[2]+j*dy;}
           for (k=0;k<baseNz;k++) {kretsch_coords[k+baseNx+baseNy] = base_bbox[4]+k*dz;}


             MPI_Allreduce((lrelkretsch0),(maxrelkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
             MPI_Allreduce((lrelkretsch0),(minrelkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
             MPI_Allreduce((lkretsch0),(maxkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
             MPI_Allreduce((lkretsch0),(minkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
   
             for (i=0; i<baseNx*baseNy*baseNz; i++)
             {
              if (uniSize>1)
              {
               relkretsch0[i]=maxrelkretsch0[i]+minrelkretsch0[i];
               kretsch0[i]=maxkretsch0[i]+minkretsch0[i];
              }
              else
              {
               relkretsch0[i]=maxrelkretsch0[i];
               kretsch0[i]=maxkretsch0[i];
              }
             }
   

           if (my_rank==0)
           {
            sprintf(name,"%srelkretsch",AMRD_save_tag);
            gft_out_full(name,ct,kretsch_shape,kretsch_cnames,kretsch_rank,kretsch_coords,relkretsch0);
//            gft_out_bbox(name,ct,kretsch_shape,kretsch_rank,base_bbox,relkretsch0);
            sprintf(name,"%skretsch",AMRD_save_tag);
            gft_out_full(name,ct,kretsch_shape,kretsch_cnames,kretsch_rank,kretsch_coords,kretsch0);
//            gft_out_bbox(name,ct,kretsch_shape,kretsch_rank,base_bbox,kretsch0);
           }
   
          } //closes condition on output_kretsch


         } //closes condition on lsteps

    } //closes condition on ct==0


      // fill in independent residual evaluator test functions
      for (i=0; i<Nx; i++)
      {
         for (j=0; j<Ny; j++)
         {
          for (k=0; k<Nz; k++)
          {
            ind=i+Nx*(j+Ny*k);
            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);

            // excise rho>1-1.5*dx pts (pure AdS diverges at rho=1, so cannot use these pts in difference stencils) 
            if (chr[ind]==AMRD_ex || 1-rho<1.5*dx_Lc)
            {
              iresall[ind]=0;  
              irestt[ind]=0;  
              irestx[ind]=0;  
              iresty[ind]=0;  
              irestz[ind]=0;
              iresxx[ind]=0;  
              iresxy[ind]=0;  
              iresxz[ind]=0;
              iresyy[ind]=0;  
              iresyz[ind]=0;
              irespsi[ind]=0;  
            }
            else
            {
              iresall[ind]=efe_all_ires[ind];
              irestt[ind]=efe_tt_ires[ind];
              irestx[ind]=efe_tx_ires[ind];
              iresty[ind]=efe_ty_ires[ind];
              irestz[ind]=efe_tz_ires[ind];
              iresxx[ind]=efe_xx_ires[ind];
              iresxy[ind]=efe_xy_ires[ind];
              iresxz[ind]=efe_xz_ires[ind];
              iresyy[ind]=efe_yy_ires[ind];
              iresyz[ind]=efe_yz_ires[ind];
              irespsi[ind]=efe_psi_ires[ind];
            }
          }
         } 
      }

//    }

//   }



   return;
}

// to be callable from fortran
void check_nan_(real *x, int *is_nan)
{
   if (isnan(*x)) *is_nan=1; else *is_nan=0;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables ... just
// use the value from the most recent evolution step
//=============================================================================
#define LIN_ZERO_BND 1
int lin_zero_bnd_all=1;
real AdS4D_evo_residual(void)
{
   real l2norm=0,l2norm_phi1,l2norm_gb,l2norm_hb_t,l2norm_hb_i;
   int is_nan;

   ldptr();

   if (LIN_ZERO_BND) 
   {
      lin_zero_bnd_res_(phi1_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&Nz);
      lin_zero_bnd_res_(gb_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&Nz);
      lin_zero_bnd_res_(hb_t_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&Nz);
      lin_zero_bnd_res_(hb_i_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny,&Nz);
   }

   l2norm_phi1=norm_l2(phi1_res,mask,chr);
   l2norm_gb=norm_l2(gb_res,mask,chr);
   l2norm_hb_t=norm_l2(hb_t_res,mask,chr);
   l2norm_hb_i=norm_l2(hb_i_res,mask,chr);

   l2norm=l2norm_phi1;

   if (!background) l2norm+=(l2norm_gb+l2norm_hb_t+l2norm_hb_i);

   check_nan_(&l2norm,&is_nan);

   if (is_nan)
   {
      printf("\nl2norm_phi1=%lf, l2norm_gb=%lf, l2norm_hb_t=%lf, l2norm_hb_i=%lf, g_norms[phi1_n_gfn-1]=%lf\n",
              l2norm_phi1,l2norm_gb,l2norm_hb_t,l2norm_hb_i,g_norms[phi1_n_gfn-1]);
      printf("[Nx,Ny,Nz]=[%i,%i,%i],L=%i\n",Nx,Ny,Nz,g_L);
      AMRD_stop("l2norm is nan ... stopping","");
      l2norm=0;
   }

   return l2norm;
}

////=============================================================================
////---lower frequency KO dissipation filter, called from AdS4D_evolve()
////   after pointers set up
////=============================================================================
//void apply_diss_eps_k(void)
//{
//   int pbt_even[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_EVEN,PAMR_UNKNOWN};
//   int pbt_unknown[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_UNKNOWN};
//   int pbt_odd[4]={PAMR_UNKNOWN,PAMR_UNKNOWN,PAMR_ODD,PAMR_UNKNOWN};
//   int even=PAMR_EVEN,odd=PAMR_ODD,i,j,ind,Nz;
//
//   // define eps as an array, and dissipate each dimension one at a time 
//   int do_ex=-1,ind_sweeps=1;
//
//   real rho;
//   dx=x[1]-x[0]; dy=y[1]-y[0];
//
//   // using w2 for eps array, with new diss_eps_y_cutoff flag (off for diss_eps_y_cutoff=1) 
//   for (i=0; i<Nx; i++)
//      for (j=0; j<Ny; j++)
//      {
//         ind=i+j*Nx;
//         rho=sqrt(x[i]*x[i]+y[j]*y[j]);
//         if (rho>=1 || y[j]>diss_eps_y_cutoff) 
//         {
//           w2[ind]=0; 
//         }
//         else 
//         {
//           w2[ind]=pow(rho,diss_eps_k_cutoff_n)*diss_eps_k;
//         }
//      } 
//
//   // only 2D Kreiss-Oliger dissipation
//   Nz=1;
//
//   // define effective excision mask w3 for dmdiss3d
//   for (i=0; i<Nx; i++)
//   {
//      for (j=0; j<Ny; j++)
//      {
//         ind=i+j*Nx;
//         rho=sqrt(x[i]*x[i]+y[j]*y[j]);
//         w3[ind]=0;
//         if (1-rho<i_shift*dx+dx/2) w3[ind]=AMRD_ex;
//      }
//   }
//
//   // use effective excision mask w3 instead of excision mask chr with new diss_all flag
//   if (diss_all==1)
//   {
//   dmdiss3d_ex_gen_(gb_tt_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(gb_tx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(gb_ty_n,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(gb_xx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(gb_xy_n,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(gb_yy_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(psi_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(Hb_t_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(Hb_x_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(Hb_y_n,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   dmdiss3d_ex_gen_(phi1_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   }
//   else {
//   dmdiss3d_ex_gen_(gb_xx_n,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   }
//
//   if (diss_all_past_k)
//   {
//      dmdiss3d_ex_gen_(gb_tt_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(gb_tx_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(gb_ty_nm1,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(gb_xx_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(gb_xy_nm1,w1,w2,&diss_bdy_k,pbt_odd,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(gb_yy_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(psi_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(Hb_t_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(Hb_x_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(Hb_y_nm1,w1,w2,&diss_bdy_k,pbt_unknown,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//      dmdiss3d_ex_gen_(phi1_nm1,w1,w2,&diss_bdy_k,pbt_even,&even,&odd,&Nx,&Ny,&Nz,w3,&AMRD_ex,&do_ex,&ind_sweeps,&diss_kmax);
//   }
//
//   return;
//}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//
// NOTE: at this point, the time sequence is: np1,n,nm1
//=============================================================================
void AdS4D_evolve(int iter)
{
   int i,zero_i=0;
   int ltrace=0;
   real ct,zero=0;

   ldptr();

   ct=PAMR_get_time(g_L);

   if (ltrace) printf("AdS4D_evolve: iter=%i , time=%lf, lev=%i, rank=%i\n",iter,ct,g_L,my_rank);

   zero_f(hb_t_res);
   zero_f(hb_i_res);

   // when kmax nonzero, apply extra dissipation
//   if (diss_kmax>=1 && diss_eps_k>0 && iter==1) apply_diss_eps_k();

   if (!background)
   {
      hb_t_evo_(hb_t_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_tz_np1,gb_tz_n,gb_tz_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_xz_np1,gb_xz_n,gb_xz_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                gb_yz_np1,gb_yz_n,gb_yz_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                Hb_z_np1,Hb_z_n,Hb_z_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,z,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,&Nz,
                Hb_t_0,Hb_x_0,Hb_y_0,
                Hb_z_0,
                &gauge_t,&ct,&rho1_t,&rho2_t,&rho3_t,&rho4_t,&xi1_t,&xi2_t,
                &c1_t,&c2_t,&c3_t,&cbulk_t);

      hb_i_evo_(hb_i_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_tz_np1,gb_tz_n,gb_tz_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,
                gb_xz_np1,gb_xz_n,gb_xz_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                gb_yz_np1,gb_yz_n,gb_yz_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                Hb_z_np1,Hb_z_n,Hb_z_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,z,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,&Nz,
                Hb_t_0,Hb_x_0,Hb_y_0,
                Hb_z_0,
                &gauge_i,&ct,&rho1_i,&rho2_i,&rho3_i,&rho4_i,&xi1_i,&xi2_i,
                &c1_i,&c2_i,&c3_i,&cbulk_i);

     g_evo_opt_(gb_res,phi1_res,cl_res,
                gb_tt_np1,gb_tt_n,gb_tt_nm1,
                gb_tx_np1,gb_tx_n,gb_tx_nm1,
                gb_ty_np1,gb_ty_n,gb_ty_nm1,
                gb_tz_np1,gb_tz_n,gb_tz_nm1,
                gb_xx_np1,gb_xx_n,gb_xx_nm1,
                gb_xy_np1,gb_xy_n,gb_xy_nm1,     
                gb_xz_np1,gb_xz_n,gb_xz_nm1,
                gb_yy_np1,gb_yy_n,gb_yy_nm1,
                gb_yz_np1,gb_yz_n,gb_yz_nm1,
                psi_np1,psi_n,psi_nm1,
                Hb_t_np1,Hb_t_n,Hb_t_nm1,
                Hb_x_np1,Hb_x_n,Hb_x_nm1,
                Hb_y_np1,Hb_y_n,Hb_y_nm1,
                Hb_z_np1,Hb_z_n,Hb_z_nm1,
                phi1_np1,phi1_n,phi1_nm1,
                &AdS_L,x,y,z,&dt,chr,&AMRD_ex,
                phys_bdy,ghost_width,&Nx,&Ny,&Nz,
                &background,&kappa_cd,&rho_cd,
                &interptype,&i_shift,&regtype,
                &diss_kmax,tfunction);
   }

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
// 
// outside rho=1 grid is also excised
//=============================================================================
void AdS4D_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
   int i,j,k,ind,l;
   real x,y,z,dx,dy,dz,rho,xp,yp,zp,ex_r_xp,ex_r_yp,ex_r_zp,r;

   dx=(bbox[1]-bbox[0])/(shape[0]-1);
   dy=(bbox[3]-bbox[2])/(shape[1]-1);
   dz=(bbox[5]-bbox[4])/(shape[2]-1);

//   printf("AdS4D_fill_ex_mask is called");

   for (i=0; i<shape[0]; i++)
   {
      x=bbox[0]+i*dx;
      for (j=0; j<shape[1]; j++)
      {
       y=bbox[2]+j*dy;
       for (k=0; k<shape[2]; k++)
       {
         z=bbox[4]+k*dz;
         rho=sqrt(x*x+y*y+z*z);
         ind=i+shape[0]*(j+shape[1]*k);

         if (rho>=(1-dx_Lc/2)) 
         {
            //excise outside rho>=1-dx
            mask[ind]=excised; 
         }
         else 
         {
            mask[ind]=excised-1;
//       printf("point not excised");
//       printf("EX_MASK: i=%i,j=%i,k=%i,ind=%i\n",i,j,k,ind);
//       printf("EX_MASK: x=%lf,y=%lf,z=%lf,rho=%lf\n",x,y,z,rho);
//       printf("EX_MASK: mask[ind]=%lf\n",mask[ind]);
            //when MAX_BHS.ne.0, excise inside a fraction ex_r of the horizon
            for (l=0; l<MAX_BHS; l++)
            {
               if (ex_r[l][0]>0)
               {
                 xp=(x-ex_xc[l][0]);
                 yp=(y-ex_xc[l][1]);
                 zp=(z-ex_xc[l][2]);
                 ex_r_xp=(ex_r[l][0]*(1-ex_rbuf[l]));
                 ex_r_yp=(ex_r[l][1]*(1-ex_rbuf[l]));
                 ex_r_zp=(ex_r[l][2]*(1-ex_rbuf[l]));

//                 printf("xp=%lf,yp=%lf,zp=%lf\n",xp,yp,zp);
//                 printf("ex_r_xp=%lf,ex_r_yp=%lf,ex_r_zp=%lf\n",ex_r_xp,ex_r_yp,ex_r_zp);

                 if ((r=sqrt(xp*xp/ex_r_xp/ex_r_xp+yp*yp/ex_r_yp/ex_r_yp+zp*zp/ex_r_zp/ex_r_zp))<1) 
                 {
                  mask[ind]=excised;
                 }
               }
            }
         }
       }
      }
   }
}


//=============================================================================
//=============================================================================
void AdS4D_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
   if (max_num<MAX_BHS) AMRD_stop("AdS4D_fill_bh_bboxes: error max_num too small\n","");

   *num=0;
}

//=============================================================================
// The following routine searches for AH's, manages excision, 
// and if t==0, determines past time level information using evolve
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
#define AH_RESET_AFTER_FAIL 0
int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
int found_count_AH[MAX_BHS];
int pre_tstep_global_first=1,valid;

int mem_alloc_bdyquantities_first=1;
int is_bdy,ie_bdy;

void AdS4D_pre_tstep(int L)
{
   char name[256];
   int AH[MAX_BHS];
   int AH_shape[2],got_an_AH,do_reinit_ex,do_repop;
   real M,J,c_equat,c_polar;
   real AH_bbox[4],AH_min_resid0,AH_min_resid1,AH_min_resid2;
   real ex_r0[3],ex_xc0[3],dt;

   real ct;
   real new_rbuf;
   real tol_save;     

   int omt;

   int n,i,j,k,ind,j_red,l,e,Lf,Lc;
   real rho;
   real rh,mh,rhoh;

   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

//   printf("AdS4D_pre_tstep is called");
//           fflush(stdout);


   if (AMRD_state!=AMRD_STATE_EVOLVE) return; // if disable, enable(?) reset_AH_shapes below

   if (pre_tstep_global_first)
   {
      for (l=0; l<MAX_BHS; l++) { AH_count[l]=found_AH[l]=found_count_AH[l]=0; freq0[l]=AH_freq[l]; }
      pre_tstep_global_first=0;
   }


     //allocate memory for quasi-local boundary stress-energy tensor before the first time step of every run. Reallocate if restarting from checkpoint
     if (L==Lc)
     {
      if (output_bdyquantities)
      {

       valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
       while(valid)
       {

        if (mem_alloc_bdyquantities_first)
        {
         mem_alloc_bdyquantities_first=0;

         ldptr();

//routine that sets a mask for near bdy points. We will call these "nexttobdypoints". The number of nexttobdypoints is also the number of points at the boundary where we will extrapolate the stress-energy tensor in AdS4D_pre_tstep and AdS4D_post_tstep. We call this number numbdypoints.
       nexttobdypoints_(chrbdy,&numbdypoints,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

         MPI_Comm_size(MPI_COMM_WORLD,&uniSize); 
//         chrbdy = malloc(size*sizeof(real)); //recall: size=Nx*Ny*Nz, where Nx/y/z are the number of grid points along x/y/z of the current process
         vecbdypoints = malloc(uniSize*sizeof(int));
         dsplsbdypoints = malloc(uniSize*sizeof(int));
 
//         numbdypoints=0; //initialize
//         for (i=0; i<size; i++)

//         }
         //routine that identifies points next to the boundary AND next to excised points. We will call these "nexttobdypoints". The number of nexttobdypoints is also the number of points at the boundary where we will extrapolate the stress-energy tensor. We call this number numbdypoints.
//         nexttobdypoints_(chrbdy,&numbdypoints,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

         //reduce boundary points if needed
//         int reduced_numbdypoints=10000;
//         if (numbdypoints*uniSize>reduced_numbdypoints) numbdypoints=roundl(reduced_numbdypoints/uniSize)+1;
 
          //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
          MPI_Allgather(&numbdypoints,1,MPI_INT,vecbdypoints,1,MPI_INT,MPI_COMM_WORLD);

          //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
          MPI_Allreduce(&numbdypoints,&basenumbdypoints,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

//          if (my_rank==0) {printf("basenumbdypoints=%i",basenumbdypoints); fflush(stdout);}


          lquasiset_tt0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_tchi0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_txi0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_chichi0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_chixi0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_xixi0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_massdensity0   = malloc((basenumbdypoints)*sizeof(real));
          lquasiset_trace0   = malloc((basenumbdypoints)*sizeof(real));
          lAdS_mass0   = malloc(sizeof(real));
          llocoeffphi10   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_tt0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_tchi0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_txi0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_chichi0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_chixi0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_xixi0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_massdensity0   = malloc((basenumbdypoints)*sizeof(real));
          maxquasiset_trace0   = malloc((basenumbdypoints)*sizeof(real));
          maxlocoeffphi10   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_tt0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_tchi0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_txi0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_chichi0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_chixi0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_xixi0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_massdensity0   = malloc((basenumbdypoints)*sizeof(real));
          minquasiset_trace0   = malloc((basenumbdypoints)*sizeof(real));
          minlocoeffphi10   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_tt0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_tchi0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_txi0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_chichi0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_chixi0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_xixi0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_massdensity0   = malloc((basenumbdypoints)*sizeof(real));
          quasiset_trace0   = malloc((basenumbdypoints)*sizeof(real));
          AdS_mass0   = malloc(sizeof(real));
          locoeffphi10   = malloc((basenumbdypoints)*sizeof(real));
 
          xextrap0   = malloc((basenumbdypoints)*sizeof(real));
          yextrap0   = malloc((basenumbdypoints)*sizeof(real));
          zextrap0   = malloc((basenumbdypoints)*sizeof(real));


          //initialize
          for (i=0;i<basenumbdypoints;i++)
          {
           lquasiset_tt0[i]           = 0;
           lquasiset_tchi0[i]         = 0;
           lquasiset_txi0[i]          = 0;
           lquasiset_chichi0[i]       = 0;
           lquasiset_chixi0[i]        = 0;
           lquasiset_xixi0[i]         = 0;
           lquasiset_massdensity0[i]  = 0;
           lquasiset_trace0[i]        = 0;
           *lAdS_mass0                = 0;
           llocoeffphi10[i]           = 0;
           maxquasiset_tt0[i]         = 0;
           maxquasiset_tchi0[i]       = 0;
           maxquasiset_txi0[i]        = 0;
           maxquasiset_chichi0[i]     = 0;
           maxquasiset_chixi0[i]      = 0;
           maxquasiset_xixi0[i]       = 0;
           maxquasiset_massdensity0[i]= 0;
           maxquasiset_trace0[i]      = 0;
           maxlocoeffphi10[i]         = 0;
           minquasiset_tt0[i]         = 0;
           minquasiset_tchi0[i]       = 0;
           minquasiset_txi0[i]        = 0;
           minquasiset_chichi0[i]     = 0;
           minquasiset_chixi0[i]      = 0;
           minquasiset_xixi0[i]       = 0;
           minquasiset_massdensity0[i]= 0;
           minquasiset_trace0[i]      = 0;
           minlocoeffphi10[i]         = 0;
           quasiset_tt0[i]            = 0;
           quasiset_tchi0[i]          = 0;
           quasiset_txi0[i]           = 0;
           quasiset_chichi0[i]        = 0;
           quasiset_chixi0[i]         = 0;
           quasiset_xixi0[i]          = 0;
           quasiset_massdensity0[i]   = 0;
           quasiset_trace0[i]         = 0;
           *AdS_mass0                 = 0;
           locoeffphi10[i]            = 0;
           xextrap0[i]                = 0;
           yextrap0[i]                = 0;
           zextrap0[i]                = 0;
          }


          //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
          is_bdy=0;
          if (my_rank==0)
          {
           ie_bdy=vecbdypoints[0];
          }
          else
          {
           for (j=0; j<my_rank; j++)
           {
            is_bdy=is_bdy+vecbdypoints[j];
           }
            ie_bdy=is_bdy+vecbdypoints[my_rank];
          }
   
          //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
          for (i=0; i<uniSize; i++)
          {
           dsplsbdypoints[i]=0;
          }
      
          for (i=0; i<uniSize; i++)
          {
           if (i!=0)
           {
            for (j=0; j<i; j++)
            {
             dsplsbdypoints[i]=dsplsbdypoints[i]+vecbdypoints[j];
            }
           }
          }

          xyzextrap_(xextrap,yextrap,zextrap,chrbdy,&numbdypoints,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);


          //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
          MPI_Allgatherv(xextrap,numbdypoints,MPI_DOUBLE,xextrap0,vecbdypoints,dsplsbdypoints,MPI_DOUBLE,MPI_COMM_WORLD);
          MPI_Allgatherv(yextrap,numbdypoints,MPI_DOUBLE,yextrap0,vecbdypoints,dsplsbdypoints,MPI_DOUBLE,MPI_COMM_WORLD);
          MPI_Allgatherv(zextrap,numbdypoints,MPI_DOUBLE,zextrap0,vecbdypoints,dsplsbdypoints,MPI_DOUBLE,MPI_COMM_WORLD);


//the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
         if (output_AdS_mass)
         {
          if (uniSize>1)
          {
           if (my_rank==0) printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
          }
          else //i.e. we're running on only 1 process
          {
           printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS ON ONLY 1 PROCESS\n");
           rhoextrap0 = malloc(sizeof(real));
           chiextrap0 = malloc((basenumbdypoints)*sizeof(real));
           xiextrap0  = malloc((basenumbdypoints)*sizeof(real));

            chixiextrap_(rhoextrap0,chiextrap0,xiextrap0,xextrap0,yextrap0,zextrap0,&basenumbdypoints);

            basebdy_Nchi=0;//initialize
            basebdy_Nxi=0; //initialize
            bdyn_(&basebdy_Nchi,&basebdy_Nxi,&basenumbdypoints,chiextrap0,xiextrap0);

           rhobdy0 = malloc(sizeof(real));
           chibdy0 = malloc(basebdy_Nchi*sizeof(real));
           xibdy0  = malloc(basebdy_Nxi*sizeof(real));

          }
         }
        
        }


         valid=PAMR_next_g();
       }
      }
     }


   // initialize qs objects at t=0, when at coarsest level L=Lc
   if (L==Lc && ct==0)
   {
     int lsteps=AMRD_lsteps[Lc-1];
     int ivecNt=AMRD_steps/AMRD_save_ivec0[3]+1; //+1 to include t=0
     real lmass,mass;


     valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
     while(valid)
     {
       ldptr();

     if ((output_ires)||(output_kretsch)||(output_relkretschcentregrid))
     {
       ires_(efe_all_ires,
          efe_tt_ires,efe_tx_ires,efe_ty_ires,
          efe_tz_ires,
          efe_xx_ires,efe_xy_ires,
          efe_xz_ires,
          efe_yy_ires,
          efe_yz_ires,
          efe_psi_ires,
          kretsch,
          relkretsch,
          relkretschcentregrid,
          gb_tt_np1,gb_tt_n,gb_tt_nm1,
          gb_tx_np1,gb_tx_n,gb_tx_nm1,
          gb_ty_np1,gb_ty_n,gb_ty_nm1,
          gb_tz_np1,gb_tz_n,gb_tz_nm1,
          gb_xx_np1,gb_xx_n,gb_xx_nm1,
          gb_xy_np1,gb_xy_n,gb_xy_nm1,
          gb_xz_np1,gb_xz_n,gb_xz_nm1,
          gb_yy_np1,gb_yy_n,gb_yy_nm1,
          gb_yz_np1,gb_yz_n,gb_yz_nm1,
          psi_np1,psi_n,psi_nm1,
          Hb_t_np1,Hb_t_n,Hb_t_nm1,
          Hb_x_np1,Hb_x_n,Hb_x_nm1,
          Hb_y_np1,Hb_y_n,Hb_y_nm1,
          Hb_z_np1,Hb_z_n,Hb_z_nm1,
          phi1_np1,phi1_n,phi1_nm1,    
          x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);


//      for (i=0; i<Nx; i++)
//      {
//         for (j=0; j<Ny; j++)
//         {
//          for (k=0; k<Nz; k++)
//          {
//            ind=i+Nx*(j+Ny*k);
//            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);
//            if (rho<1) printf("PRE_TSTEP:my_rank=%i,ct=%lf,x=%lf,y=%lf,z=%lf,rho=%lf,relkretsch[ind]=%lf\n",my_rank,ct,x[i],y[j],z[k],rho,relkretsch[ind]);
//          }
//         }
//      }

     }


     if (output_bdyquantities)
     {

          calc_locoeffphi1_(locoeffphi1,
                           phi1_np1,phi1_n,phi1_nm1,
                           xextrap,yextrap,zextrap,
                           chrbdy,&numbdypoints,
                           x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
  
         //routine that extrapolates the values of the component of the stress energy tensor at points at the boundary and the coordinates of the points at the boundary (i.e. xextrap[i]*xextrap[i]+yextrap[i]*yextrap[i]+zextrap[i]*zextrap[i]=1)
         quasiset_(quasiset_tt,quasiset_tchi,quasiset_txi,
                   quasiset_chichi,quasiset_chixi,
                   quasiset_xixi,
                   quasiset_trace,
                   quasiset_massdensity,
                   gb_tt_np1,gb_tt_n,gb_tt_nm1,
                   gb_tx_np1,gb_tx_n,gb_tx_nm1,
                   gb_ty_np1,gb_ty_n,gb_ty_nm1,
                   gb_tz_np1,gb_tz_n,gb_tz_nm1,
                   gb_xx_np1,gb_xx_n,gb_xx_nm1,
                   gb_xy_np1,gb_xy_n,gb_xy_nm1,
                   gb_xz_np1,gb_xz_n,gb_xz_nm1,
                   gb_yy_np1,gb_yy_n,gb_yy_nm1,
                   gb_yz_np1,gb_yz_n,gb_yz_nm1,
                   psi_np1,psi_n,psi_nm1,
                   xextrap,yextrap,zextrap,
                   chrbdy,&numbdypoints,
                   x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

    //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
        for (i=is_bdy; i<ie_bdy; i++)
        {
             lquasiset_tt0[i]=quasiset_tt[i-is_bdy];
             lquasiset_tchi0[i]=quasiset_tchi[i-is_bdy];
             lquasiset_txi0[i]=quasiset_txi[i-is_bdy];
             lquasiset_chichi0[i]=quasiset_chichi[i-is_bdy];
             lquasiset_chixi0[i]=quasiset_chixi[i-is_bdy];
             lquasiset_xixi0[i]=quasiset_xixi[i-is_bdy];
             lquasiset_massdensity0[i]=quasiset_massdensity[i-is_bdy];
             lquasiset_trace0[i]=quasiset_trace[i-is_bdy];
             llocoeffphi10[i]=locoeffphi1[i-is_bdy];
//             *lAdS_mass0=*AdS_mass;
         }
       }

       valid=PAMR_next_g();
     }

         if ((lsteps==0)&& output_bdyquantities)   //paste in post_tstep
         {
           // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
           //basenumbdypoints is set in AdS4D_post_init
           MPI_Allreduce(lquasiset_tt0,maxquasiset_tt0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_tchi0,maxquasiset_tchi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_txi0,maxquasiset_txi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chichi0,maxquasiset_chichi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chixi0,maxquasiset_chixi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_xixi0,maxquasiset_xixi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_massdensity0,maxquasiset_massdensity0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_trace0,maxquasiset_trace0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(llocoeffphi10,maxlocoeffphi10,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 
           MPI_Allreduce(lquasiset_tt0,minquasiset_tt0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_tchi0,minquasiset_tchi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_txi0,minquasiset_txi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chichi0,minquasiset_chichi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chixi0,minquasiset_chixi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_xixi0,minquasiset_xixi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_massdensity0,minquasiset_massdensity0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_trace0,minquasiset_trace0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(llocoeffphi10,minlocoeffphi10,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 
          for (i=0; i<basenumbdypoints; i++)
          { 
           if (uniSize>1)
           {
            quasiset_tt0[i]=maxquasiset_tt0[i]+minquasiset_tt0[i];
            quasiset_tchi0[i]=maxquasiset_tchi0[i]+minquasiset_tchi0[i];
            quasiset_txi0[i]=maxquasiset_txi0[i]+minquasiset_txi0[i];
            quasiset_chichi0[i]=maxquasiset_chichi0[i]+minquasiset_chichi0[i];
            quasiset_chixi0[i]=maxquasiset_chixi0[i]+minquasiset_chixi0[i];
            quasiset_xixi0[i]=maxquasiset_xixi0[i]+minquasiset_xixi0[i];
            quasiset_massdensity0[i]=maxquasiset_massdensity0[i]+minquasiset_massdensity0[i];
            quasiset_trace0[i]=maxquasiset_trace0[i]+minquasiset_trace0[i];
            locoeffphi10[i]=maxlocoeffphi10[i]+minlocoeffphi10[i];
           }
           else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
           {
            quasiset_tt0[i]=maxquasiset_tt0[i];
            quasiset_tchi0[i]=maxquasiset_tchi0[i];
            quasiset_txi0[i]=maxquasiset_txi0[i];
            quasiset_chichi0[i]=maxquasiset_chichi0[i];
            quasiset_chixi0[i]=maxquasiset_chixi0[i];
            quasiset_xixi0[i]=maxquasiset_xixi0[i];
            quasiset_massdensity0[i]=maxquasiset_massdensity0[i];
            quasiset_trace0[i]=maxquasiset_trace0[i];
            locoeffphi10[i]=maxlocoeffphi10[i];
           }
          }
//            MPI_Allreduce((lAdS_mass0),(AdS_mass0),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

          if (my_rank==0)
          {

           FILE * fp;
           fp = fopen ("ascii_xext_yext_zext_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %24.16e %i \n",xextrap0[j],yextrap0[j],zextrap0[j],j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_xext_yext_zext_indbdypoint.txt", "a+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %24.16e %i \n",xextrap0[j],yextrap0[j],zextrap0[j],j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


        if (alltimes_ascii)
        {


           fp = fopen ("ascii_t_bdyphi1_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_bdyphi1_indbdypoint.txt", "a+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           // save quasiset_ll as ascii
           fp = fopen ("ascii_t_quasisetll_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisetll_indbdypoint.txt", "a+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           fp = fopen ("ascii_t_quasisettrace_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisettrace_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           fp = fopen ("ascii_t_quasisetmassdensity_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisetmassdensity_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

        } //closes if(alltimes_ascii) condition


        if (timestep_ascii)
        {
           char filename[64];
           sprintf(filename, "ascii_t_bdyphi1_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_bdyphi1_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           // save quasiset_ll as ascii
           sprintf(filename, "ascii_t_quasisetll_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisetll_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           sprintf(filename, "ascii_t_quasisettrace_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisettrace_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           sprintf(filename, "ascii_t_quasisetmassdensity_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisetmassdensity_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

        } //closes if (timestep_ascii) condition


      } //closes if (my_rank==0) condition

//the following bit computes and prints AdS_mass0 (see below) if we're running on only 1 process
          if (output_AdS_mass)
          { 
           if (uniSize>1)
           {
            if (my_rank==0) printf("\nTHE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT RELIABLE...NOT COMPUTING AdS MASS\n");
           }
           else //i.e. we're running on only 1 process
           {
            printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");

            *rhobdy0=1;
            chibdy_xibdy_(chibdy0,xibdy0,xextrap0,yextrap0,zextrap0,&basenumbdypoints,chiextrap0,xiextrap0,&basebdy_Nchi,&basebdy_Nxi);    
            doubleintegralonsphere_(AdS_mass0,quasiset_massdensity0,xextrap0,yextrap0,zextrap0,&basenumbdypoints,rhobdy0,chibdy0,xibdy0,&basebdy_Nchi,&basebdy_Nxi); 
   
            FILE * fp;
            fp = fopen ("ascii_t_AdSmass.txt", "a+");
                 fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0);
            fclose(fp);
           }
          }

        }
    }


   // search for AHs at t>0
   do_repop=do_reinit_ex=got_an_AH=0;

   for (l=0; l<MAX_BHS; l++)
   {
      real prev_AH_R[AH_Nchi[l]*AH_Nphi[l]];
      real prev_AH_xc[3];
      for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) {prev_AH_R[i]=AH_R[l][i];} 
      for (i=0; i<3; i++) {prev_AH_xc[i]=AH_xc[l][i];}
      c_AH=l;
      if (AH_max_iter[l]>0 && L==AH_Lmin[l] && ct>=AH_tmin[l])
      {
         if (AH_count[l]<0) { AH[l]=1; if (AMRD_state==AMRD_STATE_EVOLVE) M=J=0;}
         else if (!(AH_count[l] % freq0[l]) && !(c_AH==3 && ct==0)) // for fourth BH, do not search at t=0
         {
            omt=0; // over-max-tolerance
            AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0); // AH finder

            if (AH[l]) { freq0[l]=AH_freq_aft[l]; found_AH[l]=1; got_an_AH=1; AH_tol[l]=AH_tol_aft[l]; found_count_AH[l]++; } // if this time found AH

            // if previously found but failed now
            if (found_AH[l] && !AH[l]) 
            {
               if (AH_reset_scale[l]>0)
               {
                  // expand old initial-guess surface
                  if (my_rank==0 && AMRD_evo_trace>=1)
                     printf("t=%lf ... lost AH[%i] ...\n" 
                            "expanding old initial-guess surface by a factor %lf\n",ct,l+1,AH_reset_scale[l]);
                  for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]*AH_reset_scale[l];

                  if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1)))
                  { 
                     // shrink old initial-guess surface
                     if (my_rank==0 && AMRD_evo_trace>=1) printf("... still can't find one (min_resid=%lf)\n"
                         "... shrinking old initial-guess surface by a factor %lf\n",AH_min_resid1,1/AH_reset_scale[l]);
                     for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]/AH_reset_scale[l]; 

                     if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2)))
                     {

                        // increase AH_tol to force an AH to be found, starting with shrunken initial-guess surface
                        if (AH_min_resid2<AH_min_resid1 && AH_min_resid2<AH_min_resid0) 
                        {
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from shrunken initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]/AH_reset_scale[l]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; will then use the resulting surface\n"
                              ,AH_min_resid2*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid2*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                        // increase AH_tol to force an AH to be found, starting with expanded initial-guess surface
                        else if (AH_min_resid1<AH_min_resid0) 
                        {  
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from expanded initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]*AH_reset_scale[l]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; use result surface\n",AH_min_resid1*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid1*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                        // increase AH_tol to force an AH to be found, starting with old initial-guess surface
                        else 
                        {
                           if (my_rank==0 && AMRD_evo_trace>=1)
                              printf("starting from old initial-guess surface\n");
                           for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]; 
                           if (my_rank==0 && AMRD_evo_trace>=1) 
                              printf("and temporarily increasing tolerance to %lf; use result surface\n",AH_min_resid0*AH_omt_scale[l]);
                           omt=1;
                           tol_save=AH_tol[l];
                           AH_tol[l]=AH_min_resid0*AH_omt_scale[l];
                           if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0)))
                           {
                              if (my_rank==0 && AMRD_evo_trace>=1) printf("BUG: couldn't find *same* AH\n");
                              if (AH_RESET_AFTER_FAIL) found_AH[l]=0;
                           }
                           if (my_rank==0 && AMRD_evo_trace>=1) printf("setting tolerance back to %lf\n",tol_save);
                           AH_tol[l]=tol_save;
                        }

                     }                 
                  }
               }
            }

            // if still not found
            if (!(found_AH[l])) for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) {AH_R[l][i]=prev_AH_R[i];}

            // if never found AH
            if (!(found_AH[l])) for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=AH_r0[l];

            // save AH grid functions if this time found AH
            if (AH[l]) 
            {
               AH_shape[0]=AH_Nchi[l];
               AH_shape[1]=AH_Nphi[l];
               AH_bbox[0]=0;
               AH_bbox[1]=M_PI;
               AH_bbox[2]=0;
               if (AH_xc[l][1]<dy) {AH_bbox[3]=M_PI;} else {AH_bbox[3]=2*M_PI;}
               int rank=2;
               sprintf(name,"%sAH_R_%i",AMRD_save_tag,l);
               gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_R[l]);
               sprintf(name,"%sAH_theta_%i",AMRD_save_tag,l);
               gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_theta[l]);
            }

            // fill in excision parameters 
            // ( ex_xc0[0],ex_xc0[1],ex_xc0[2] are filled with coordinate center for excision
            //   and ex_r0[0],ex_r0[1],ex_r0[2] are filled with principle axis radii for excision )
            if (found_AH[l] && AH[l] && AMRD_do_ex)
            {
               fill_ex_params_(AH_R[l],AH_xc[l],ex_r0,ex_xc0,&AH_Nchi[l],&AH_Nphi[l],&dx,&dy,&dz,&axisym);

               if (no_AH_intersect(ex_r0,ex_xc0,l))
               {
                 do_reinit_ex=1;
                 do_repop=1;
                 // saves local ex_r0,ex_xc0 to global ex_r, ex_xc
                 ex_r[l][0]=ex_r0[0]; ////AH ellipse x-semiaxis
                 ex_r[l][1]=ex_r0[1]; //AH ellipse y-semiaxis
                 ex_r[l][2]=ex_r0[2]; //AH ellipse z-semiaxis
                 ex_xc[l][0]=ex_xc0[0]; //excision ellipse x-coordinate center
                 ex_xc[l][1]=ex_xc0[1]; //excision ellipse y-coordinate center 
                 ex_xc[l][2]=ex_xc0[2]; //excision ellipse z-coordinate center
               }

               if (my_rank==0) 
               {
                 printf("AH ellipse x/y/z-semiaxis: (ex_r[0],ex_r[1],ex_r[2])=(%lf,%lf,%lf)\n",ex_r[l][0],ex_r[l][1],ex_r[l][2]);
                 printf("excision ellipse semiaxes: (ex_r[0]*(1-ex_rbuf[l]),ex_r[1]*(1-ex_rbuf[l]),ex_r[2]*(1-ex_rbuf[l]))=(%lf,%lf,%lf)\nexcision ellipse center:   (ex_xc[0],ex_xc[1],ex_xc[2])=(%lf,%lf,%lf)\n",ex_r[l][0]*(1-ex_rbuf[l]),ex_r[l][1]*(1-ex_rbuf[l]),ex_r[l][2]*(1-ex_rbuf[l]),ex_xc[l][0],ex_xc[l][1],ex_xc[l][2]);
                printf("Excision buffer (i.e. size of the evolved region within the AH) ex_rbuf[0]=%lf\n",ex_rbuf[0]);
               }
            }

         }
         AH_count[l]++;
      }
   }

   // repopulate if needed
   int repop_n=1;
   if (do_repop) AMRD_repopulate(repop_n,ex_repop_io); 

   do_reinit_ex=1; //REMOVE THIS LATER (but: not doing PAMR_excision_on at all causes problems)

   // re-initialize mask function
   if (do_reinit_ex)
   {
//     remove_redundant_AH();
     PAMR_excision_on("chr",&AdS4D_fill_ex_mask,AMRD_ex,1);
   }

     // displays initial black hole radius and excision radius
     if (my_rank==0)
     {
        if (ief_bh_r0!=0) 
        {
         rh=-pow(AdS_L,2)
            /(pow(3,(1.0/3.0)))
            /(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            +(pow((9*pow(AdS_L,2)*(ief_bh_r0/2)+sqrt(3.0)*sqrt(pow(AdS_L,6)+27*pow(AdS_L,4)*pow((ief_bh_r0/2),2))),(1.0/3.0)))
            /(pow(3,(2.0/3.0)));
         mh=ief_bh_r0/2;
         rhoh=(-1 + sqrt(1 + pow(rh,2)))/rh;
//         ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh*(1-ex_rbuf[0]);
         printf("\n ... we started with a BH of mass mh=%lf, Schwarzschild radius rh=%lf and compactified radius rhoh=%lf. \n",mh,rh,rhoh);
         printf("Excision buffer (i.e. size of the evolved region within the AH) ex_rbuf[0]=%lf\n",ex_rbuf[0]);
        }
     }



   return;
}

//=============================================================================
// The following routine prints diagnostic quantities
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void AdS4D_post_tstep(int L)
{
   int itrace=1,valid;
   static int local_first = 1;

   real ct;
   int n,i,j,k,j_red,Lf,Lc;
   int is,ie,js,je,ks,ke;

//   printf("AdS4D_post_tstep is called");
//           fflush(stdout);

   ct = PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; elise Lc=1;


   // qs objects at t>0, when at coarsest level L=Lc
   if (L==Lc)
   {
     int lsteps=AMRD_lsteps[Lc-1];
     int ivecNt=AMRD_steps/AMRD_save_ivec0[3]+1; //+1 to include t=0

     valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
     while(valid)
     {
       ldptr();

     if ((output_ires)||(output_kretsch)||(output_relkretschcentregrid))
     {

       ires_(efe_all_ires,
          efe_tt_ires,efe_tx_ires,efe_ty_ires,
          efe_tz_ires,
          efe_xx_ires,efe_xy_ires,
          efe_xz_ires,
          efe_yy_ires,
          efe_yz_ires,
          efe_psi_ires,
          kretsch,
          relkretsch, 
          relkretschcentregrid,
          gb_tt_n,gb_tt_nm1,gb_tt_np1,
          gb_tx_n,gb_tx_nm1,gb_tx_np1,
          gb_ty_n,gb_ty_nm1,gb_ty_np1,
          gb_tz_n,gb_tz_nm1,gb_tz_np1,
          gb_xx_n,gb_xx_nm1,gb_xx_np1,
          gb_xy_n,gb_xy_nm1,gb_xy_np1,
          gb_xz_n,gb_xz_nm1,gb_xz_np1,
          gb_yy_n,gb_yy_nm1,gb_yy_np1,
          gb_yz_n,gb_yz_nm1,gb_yz_np1,
          psi_n,psi_nm1,psi_np1,
          Hb_t_n,Hb_t_nm1,Hb_t_np1,
          Hb_x_n,Hb_x_nm1,Hb_x_np1,
          Hb_y_n,Hb_y_nm1,Hb_y_np1,
          Hb_z_n,Hb_z_nm1,Hb_z_np1,
          phi1_n,phi1_nm1,phi1_np1,     
          x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
      }

//we save here the values of the Kretschmann scalar at the centre of the grid at times greater than 0.
//The value of the Kretschmann scalar at the centre of the grid at t=0 is saved in pre_io_calc.

         if (output_relkretschcentregrid)
         {
          *lrelkretschcentregrid0= *relkretschcentregrid;
         }

         if (output_kretsch)
         {

          is=((bbox[0]-base_bbox[0])/dx);    // for "left-most"  processors, includes i=0
          ie=((bbox[1]-base_bbox[0])/dx)+1;  // for "right-most" processors, includes i=Nx
          js=((bbox[2]-base_bbox[2])/dy);    // for "left-most"  processors, includes j=0
          je=((bbox[3]-base_bbox[2])/dy)+1;  // for "right-most" processors, includes j=Ny
          ks=((bbox[4]-base_bbox[4])/dz);    // for "left-most"  processors, includes k=0
          ke=((bbox[5]-base_bbox[4])/dz)+1;  // for "right-most" processors, includes k=Ny

           baseNx=AMRD_base_shape[0];
           baseNy=AMRD_base_shape[1];
           baseNz=AMRD_base_shape[2];

           for (i=is; i<ie; i++)
           {
            for (j=js; j<je; j++)
            {
             for (k=ks; k<ke; k++)
             {
              lrelkretsch0[i+baseNx*(j+baseNy*k)]=relkretsch[(i-is)+Nx*((j-js)+Ny*(k-ks))];
              lkretsch0[i+baseNx*(j+baseNy*k)]=kretsch[(i-is)+Nx*((j-js)+Ny*(k-ks))];
             }
            }
           }
         }



      if (output_bdyquantities)
      {

       calc_locoeffphi1_(locoeffphi1,
                           phi1_n,phi1_nm1,phi1_np1,
                           xextrap,yextrap,zextrap,
                           chrbdy,&numbdypoints,
                           x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

       quasiset_(quasiset_tt,quasiset_tchi,quasiset_txi,
                 quasiset_chichi,quasiset_chixi,
                 quasiset_xixi,
                 quasiset_trace,
                 quasiset_massdensity,
                 gb_tt_n,gb_tt_nm1,gb_tt_np1,
                 gb_tx_n,gb_tx_nm1,gb_tx_np1,
                 gb_ty_n,gb_ty_nm1,gb_ty_np1,
                 gb_tz_n,gb_tz_nm1,gb_tz_np1,
                 gb_xx_n,gb_xx_nm1,gb_xx_np1,
                 gb_xy_n,gb_xy_nm1,gb_xy_np1,
                 gb_xz_n,gb_xz_nm1,gb_xz_np1,
                 gb_yy_n,gb_yy_nm1,gb_yy_np1,
                 gb_yz_n,gb_yz_nm1,gb_yz_np1,
                 psi_n,psi_nm1,psi_np1,
                 xextrap,yextrap,zextrap,
                 chrbdy,&numbdypoints,
                 x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);

     //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywher except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution.
         for (i=is_bdy; i<ie_bdy; i++)
         {
             lquasiset_tt0[i]=quasiset_tt[i-is_bdy];
             lquasiset_tchi0[i]=quasiset_tchi[i-is_bdy];
             lquasiset_txi0[i]=quasiset_txi[i-is_bdy];
             lquasiset_chichi0[i]=quasiset_chichi[i-is_bdy];
             lquasiset_chixi0[i]=quasiset_chixi[i-is_bdy];
             lquasiset_xixi0[i]=quasiset_xixi[i-is_bdy];
             lquasiset_massdensity0[i]=quasiset_massdensity[i-is_bdy];
             lquasiset_trace0[i]=quasiset_trace[i-is_bdy];
             llocoeffphi10[i]=locoeffphi1[i-is_bdy];
//             *lAdS_mass0=*AdS_mass;
         }

      } //closes condition on output_bdyquantities

       valid=PAMR_next_g();
     }

  
         if ((lsteps%AMRD_save_ivec0[3]==0)&&(lsteps!=0))   //paste in post_tstep
         {

          if (output_relkretschcentregrid)
          {
            MPI_Allreduce((lrelkretschcentregrid0),(maxrelkretschcentregrid0),1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce((lrelkretschcentregrid0),(minrelkretschcentregrid0),1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
            if (uniSize>1)
            {
             *relkretschcentregrid0=*maxrelkretschcentregrid0+*minrelkretschcentregrid0;
            }
            else
            {
             *relkretschcentregrid0=*maxrelkretschcentregrid0;
            }
             if (my_rank==0)
             {
             // save relkretschcentregrid as ascii
             FILE * fp;
             fp = fopen ("ascii_t_relkretschcentregrid.txt", "a+");
                   fprintf(fp,"%24.16e %24.16e \n",ct,*relkretschcentregrid0);
             fclose(fp);
             }
          }

          if (output_kretsch)
          {


           kretsch_rank=3;
           kretsch_shape[0]=baseNx;
           kretsch_shape[1]=baseNy;
           kretsch_shape[2]=baseNz;
           sprintf(kretsch_cnames,"x|y|z");
           kretsch_coords = malloc(sizeof(double)*(baseNx+baseNy+baseNz));
           for (i=0;i<baseNx;i++) {kretsch_coords[i] = base_bbox[0]+i*dx;}
           for (j=0;j<baseNy;j++) {kretsch_coords[j+baseNx] = base_bbox[2]+j*dy;}
           for (k=0;k<baseNz;k++) {kretsch_coords[k+baseNx+baseNy] = base_bbox[4]+k*dz;}

            MPI_Allreduce((lrelkretsch0),(maxrelkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce((lrelkretsch0),(minrelkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
            MPI_Allreduce((lkretsch0),(maxkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
            MPI_Allreduce((lkretsch0),(minkretsch0),baseNx*baseNy*baseNz,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

  
            for (i=0; i<baseNx*baseNy*baseNz; i++)
            {
             if (uniSize>1)
             {
              relkretsch0[i]=maxrelkretsch0[i]+minrelkretsch0[i];
              kretsch0[i]=maxkretsch0[i]+minkretsch0[i];
             }
             else
             {
              relkretsch0[i]=maxrelkretsch0[i];
              kretsch0[i]=maxkretsch0[i];
             }
            }

            if (my_rank==0)
            {
             sprintf(name,"%srelkretsch",AMRD_save_tag);
             gft_out_full(name,ct,kretsch_shape,kretsch_cnames,kretsch_rank,kretsch_coords,relkretsch0);
//             gft_out_bbox(name,ct,kretsch_shape,kretsch_rank,base_bbox,relkretsch0);
             sprintf(name,"%skretsch",AMRD_save_tag);
             gft_out_full(name,ct,kretsch_shape,kretsch_cnames,kretsch_rank,kretsch_coords,kretsch0);
//             gft_out_bbox(name,ct,kretsch_shape,kretsch_rank,base_bbox,kretsch0);

            }


          }



          if (output_bdyquantities)
          {

           // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
           //basenumbdypoints is set in AdS4D_post_init
           MPI_Allreduce(lquasiset_tt0,maxquasiset_tt0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_tchi0,maxquasiset_tchi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_txi0,maxquasiset_txi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chichi0,maxquasiset_chichi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chixi0,maxquasiset_chixi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_xixi0,maxquasiset_xixi0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_massdensity0,maxquasiset_massdensity0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_trace0,maxquasiset_trace0,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
           MPI_Allreduce(llocoeffphi10,maxlocoeffphi10,basenumbdypoints,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

           MPI_Allreduce(lquasiset_tt0,minquasiset_tt0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_tchi0,minquasiset_tchi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_txi0,minquasiset_txi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chichi0,minquasiset_chichi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_chixi0,minquasiset_chixi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_xixi0,minquasiset_xixi0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_massdensity0,minquasiset_massdensity0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(lquasiset_trace0,minquasiset_trace0,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
           MPI_Allreduce(llocoeffphi10,minlocoeffphi10,basenumbdypoints,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);

          for (i=0; i<basenumbdypoints; i++)
          {
           if (uniSize>1)
           {
            quasiset_tt0[i]=maxquasiset_tt0[i]+minquasiset_tt0[i];
            quasiset_tchi0[i]=maxquasiset_tchi0[i]+minquasiset_tchi0[i];
            quasiset_txi0[i]=maxquasiset_txi0[i]+minquasiset_txi0[i];
            quasiset_chichi0[i]=maxquasiset_chichi0[i]+minquasiset_chichi0[i];
            quasiset_chixi0[i]=maxquasiset_chixi0[i]+minquasiset_chixi0[i];
            quasiset_xixi0[i]=maxquasiset_xixi0[i]+minquasiset_xixi0[i];
            quasiset_massdensity0[i]=maxquasiset_massdensity0[i]+minquasiset_massdensity0[i];
            quasiset_trace0[i]=maxquasiset_trace0[i]+minquasiset_trace0[i];
            locoeffphi10[i]=maxlocoeffphi10[i]+minlocoeffphi10[i];
           }
           else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
           {
            quasiset_tt0[i]=maxquasiset_tt0[i];
            quasiset_tchi0[i]=maxquasiset_tchi0[i];
            quasiset_txi0[i]=maxquasiset_txi0[i];
            quasiset_chichi0[i]=maxquasiset_chichi0[i];
            quasiset_chixi0[i]=maxquasiset_chixi0[i];
            quasiset_xixi0[i]=maxquasiset_xixi0[i];
            quasiset_massdensity0[i]=maxquasiset_massdensity0[i];
            quasiset_trace0[i]=maxquasiset_trace0[i];
            locoeffphi10[i]=maxlocoeffphi10[i];
           }
          }
//            MPI_Allreduce((lAdS_mass0),(AdS_mass0),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


          if (my_rank==0)
          {

           FILE * fp;
           // save quasiset_ll as ascii

        if (alltimes_ascii)
        {




           fp = fopen ("ascii_t_bdyphi1_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {

                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,locoeffphi10[j],j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_bdyphi1_indbdypoint.txt", "a+");
            j_red=0; 
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if (j%reduction_factor==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           fp = fopen ("ascii_t_quasisetll_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisetll_indbdypoint.txt", "a+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           fp = fopen ("ascii_t_quasisettrace_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisettrace_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


           fp = fopen ("ascii_t_quasisetmassdensity_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           fp = fopen ("ascii_reduced_t_quasisetmassdensity_indbdypoint.txt", "a+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }


        }  //closes if (alltimes_ascii) condition

        if (timestep_ascii)
        {
           char filename[64];
           sprintf(filename, "ascii_t_bdyphi1_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_bdyphi1_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            j_red=0; 
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",ct,locoeffphi10[j],j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           // save quasiset_ll as ascii
           sprintf(filename, "ascii_t_quasisetll_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisetll_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            j_red=0;
            for( j = 0; j < basenumbdypoints; j++ )
              {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                            ct,
                            quasiset_tt0[j],quasiset_tchi0[j],quasiset_txi0[j],quasiset_chichi0[j],quasiset_chixi0[j],quasiset_xixi0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           sprintf(filename, "ascii_t_quasisettrace_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisettrace_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_trace0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

           sprintf(filename, "ascii_t_quasisetmassdensity_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
              {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j);
              }
           fclose(fp);

          if (reduced_ascii)
          {
           sprintf(filename, "ascii_reduced_t_quasisetmassdensity_indbdypoint_tstep%d.txt", lsteps);
           fp = fopen (filename, "w+");
            for( j = 0; j < basenumbdypoints; j++ )
            {
               if ((j%reduction_factor)==0)
               {
                fprintf(fp,"%24.16e %24.16e %i \n",
                            ct,
                            quasiset_massdensity0[j],
                            j_red);
                j_red=j_red+1;
               }
              }
           fclose(fp);
          }

        }  //closes if (timestep_ascii) condition 

       }  //closes condition my_rank==0

//the following bit computes and prints AdS_mass0 (see below) if we're running on only 1 process
          if (output_AdS_mass)
          {
           if (uniSize>1)
           {
            if (my_rank==0) printf("\nrunning on more than 1 process...not computing AdS mass\n");
           }
           else //i.e. we're running on only 1 process
           {
            printf("\nrunning on 1 process...computing AdS mass");

            *rhobdy0=1;
            chibdy_xibdy_(chibdy0,xibdy0,xextrap0,yextrap0,zextrap0,&basenumbdypoints,chiextrap0,xiextrap0,&basebdy_Nchi,&basebdy_Nxi);
            doubleintegralonsphere_(AdS_mass0,quasiset_massdensity0,xextrap0,yextrap0,zextrap0,&basenumbdypoints,rhobdy0,chibdy0,xibdy0,&basebdy_Nchi,&basebdy_Nxi);

            FILE * fp;
            fp = fopen ("ascii_t_AdSmass.txt", "a+");
                 fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0);
            fclose(fp);
           }
          }


         } //closes if condition on output_bdyquantities
        } //closes if condition on lsteps

   } //closes if condition on L==Lc 

   if (my_rank==0) printf("\n===================================================================\n");


   if (AMRD_state!=AMRD_STATE_EVOLVE) return; // if disable, enable(?) reset_AH_shapes below


   return;
}


#define ACTION_RELAX 1
#define ACTION_LOP 3
#define ACTION_RESIDUAL 2
//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real AdS4D_MG_residual(void)
{
   int action=ACTION_RESIDUAL;
   real norm;

   ldptr_mg();

   if (background || skip_constraints)
   {
      zero_f(zeta_res);
      return 0;
   }

   // solves for zeta conformal factor at t=0; residual 
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex, 
           x,y,z,&norm,&Nx,&Ny,&Nz);

   return norm;
}


//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real AdS4D_MG_relax(void)
{
   int action=ACTION_RELAX;
   real norm;

   ldptr_mg();

   if (background || skip_constraints)
   {
      const_f(zeta,1);
      return 0;
    }

   // solves for zeta conformal factor at t=0; relaxation 
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,
           x,y,z,&norm,&Nx,&Ny,&Nz);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void AdS4D_L_op(void)
{
   int action=ACTION_LOP;
   real norm;

   ldptr_mg();

   if (background || skip_constraints) 
   {
      zero_f(zeta_lop);
      return;
   }

   // solves for zeta conformal factor at t=0; elliptic operator
   mg_sup_(&action,zeta,zeta_rhs,zeta_lop,zeta_res,phi1,
           &AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,
           x,y,z,&norm,&Nx,&Ny,&Nz);

   return;
}

//=============================================================================
//=============================================================================
void AdS4D_scale_tre(void)
{
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void AdS4D_post_regrid(void)
{
   int i;

   if (!background || skip_constraints) return;

   ldptr();

   init_ghb_ads_(gb_tt_n,gb_tx_n,gb_ty_n,
                 gb_tz_n,
                 gb_xx_n,gb_xy_n,
                 gb_xz_n,
                 gb_yy_n,
                 gb_yz_n,
                 psi_n,Hb_t_n,Hb_x_n,Hb_y_n,
                 Hb_z_n,
                 &AdS_L,x,y,z,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);

   for (i=0; i<size; i++)
   {
      gb_tt_nm1[i]=gb_tt_np1[i]=gb_tt_n[i];
      gb_tx_nm1[i]=gb_tx_np1[i]=gb_tx_n[i];
      gb_ty_nm1[i]=gb_ty_np1[i]=gb_ty_n[i];
      gb_tz_nm1[i]=gb_tz_np1[i]=gb_tz_n[i];
      gb_xx_nm1[i]=gb_xx_np1[i]=gb_xx_n[i];
      gb_xy_nm1[i]=gb_xy_np1[i]=gb_xy_n[i];
      gb_xz_nm1[i]=gb_xz_np1[i]=gb_xz_n[i];
      gb_yy_nm1[i]=gb_yy_np1[i]=gb_yy_n[i];
      gb_yz_nm1[i]=gb_yz_np1[i]=gb_yz_n[i];
      psi_nm1[i]=psi_np1[i]=psi_n[i];
   }
}

//=============================================================================
//check-pointing
//=============================================================================
#define CP_DATA_SIZE 50000
void AdS4D_copy_block(char *p, char **q, int n, int dir, int *tot_size)
{
   char *p0;
   int n0;

   if (n==0) return;

   if ((*tot_size+n) > (CP_DATA_SIZE))
      AMRD_stop("AdS4D_copy_block: error ... CP_DATA_SIZE too small\n","");
   *tot_size+=n;

   n0=n;
   p0=p;

   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*p0++;
   else while(n0--) *p0++=*(*q)++;

   return;
}

void AdS4D_cp(int dir, char *data)
{
   int size=0;
   
   if (dir==AMRD_CP_SAVE)
   {
       cp_version=ADS5D_CP_VERSION;
       AdS4D_copy_block((char *)&cp_version,&data,sizeof(int),dir,&size);
   }
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd_set_app_user_cp_hook(AdS4D_cp,CP_DATA_SIZE);
   amrd_set_app_pre_tstep_hook(AdS4D_pre_tstep);
   amrd_set_elliptic_vars_t0_init(AdS4D_elliptic_vars_t0_init);
   amrd(argc,argv,&AdS4D_id,&AdS4D_var_pre_init,
        &AdS4D_var_post_init, &AdS4D_AMRH_var_clear,
        &AdS4D_free_data, &AdS4D_t0_cnst_data,
        &AdS4D_evo_residual, &AdS4D_MG_residual,
        &AdS4D_evolve, &AdS4D_MG_relax, &AdS4D_L_op, 
        &AdS4D_pre_io_calc, &AdS4D_scale_tre, 
        &AdS4D_post_regrid, &AdS4D_post_tstep,
        &AdS4D_fill_ex_mask, &AdS4D_fill_bh_bboxes);
}

