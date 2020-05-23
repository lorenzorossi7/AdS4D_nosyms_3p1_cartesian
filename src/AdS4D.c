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
int output_ires,output_relkretschcentregrid,output_kretsch,output_relkretsch;
int output_metricatAH;
int output_bdyquantities,output_AdS_mass;
int reduced_ascii,reduction_factor_ascii;
int alltimes_ascii,timestep_ascii;
int bdy_extrap_freepts;
int bdy_extrap_fixedpts;
int output_bdy_extraporder1;
int output_bdy_extraporder2;
int output_bdy_extraporder3;
int output_bdy_extraporder4;
int bdy_extrap_order;
real ratio_Lhighres_Llowres;
int resolution_degree;
int max_resolution_degree;
int reduction_factor;

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

//variables that are defined in hyperbolic_vars but just because we want them to be synchronized (synchronization only affects the current time level)
real *relkretsch,*relkretsch_n,*relkretsch_np1,*relkretsch_nm1;

real *w1,*mg_w1;
real *w2,*mg_w2;
real *w3,*mg_w3;
real *w4,*mg_w4;

real *mask,*mask_mg,*chr,*chr_mg;

real *chrbdy_freepts_extraporder1;
real *chrbdy_freepts_extraporder2;
real *chrbdy_fixedpts_extraporder1;
real *chrbdy_fixedpts_extraporder2;

real *kg_ires,*alpha,*ricci,*theta,*f,*K;

real *phi1_res,*gb_res;
real *efe_all_ires;
real *efe_tt_ires,*efe_tx_ires,*efe_ty_ires;
real *efe_tz_ires;
real *efe_xx_ires,*efe_xy_ires,*efe_yy_ires,*efe_psi_ires;
real *efe_xz_ires,*efe_yz_ires;
int currentres_reduction_factor;
real currentres_ratio_Lhighres_Llowres;
int ind_distance_fixedpts;
int num_fixed_coords;
real *fixed_coords;

int numbdypoints_freepts_extraporder1;
int basenumbdypoints_freepts_extraporder1;
int basebdy_Nchi_freepts_extraporder1,basebdy_Nxi_freepts_extraporder1;

int numbdypoints_freepts_extraporder2;
int basenumbdypoints_freepts_extraporder2;
int basebdy_Nchi_freepts_extraporder2,basebdy_Nxi_freepts_extraporder2;

int numbdypoints_fixedpts_extraporder1;
int basenumbdypoints_fixedpts_extraporder1;
int basebdy_Nchi_fixedpts_extraporder1,basebdy_Nxi_fixedpts_extraporder1;

int numbdypoints_fixedpts_extraporder2;
int basenumbdypoints_fixedpts_extraporder2;
int basebdy_Nchi_fixedpts_extraporder2,basebdy_Nxi_fixedpts_extraporder2;

int *vecbdypoints_freepts_extraporder1, *dsplsbdypoints_freepts_extraporder1;
int *vecbdypoints_freepts_extraporder2, *dsplsbdypoints_freepts_extraporder2;
int *vecbdypoints_fixedpts_extraporder1, *dsplsbdypoints_fixedpts_extraporder1;
int *vecbdypoints_fixedpts_extraporder2, *dsplsbdypoints_fixedpts_extraporder2;

int uniSize;
real *leadordcoeff_phi1;
real *quasiset_tt_ll,*quasiset_tchi_ll,*quasiset_txi_ll;
real *quasiset_chichi_ll,*quasiset_chixi_ll,*quasiset_xixi_ll;
real *quasiset_tracell, *quasiset_massdensityll;

real *quasiset_tt_freepts_extraporder1;
real *quasiset_tchi_freepts_extraporder1;
real *quasiset_txi_freepts_extraporder1;
real *quasiset_chichi_freepts_extraporder1;
real *quasiset_chixi_freepts_extraporder1;
real *quasiset_xixi_freepts_extraporder1;
real *quasiset_trace_freepts_extraporder1;
real *quasiset_massdensity_freepts_extraporder1;
real *AdS_mass_freepts_extraporder1;
real *bdyphi_freepts_extraporder1;
real *xextrap_freepts_extraporder1;
real *yextrap_freepts_extraporder1;
real *zextrap_freepts_extraporder1;

real *quasiset_tt_freepts_extraporder2;
real *quasiset_tchi_freepts_extraporder2;
real *quasiset_txi_freepts_extraporder2;
real *quasiset_chichi_freepts_extraporder2;
real *quasiset_chixi_freepts_extraporder2;
real *quasiset_xixi_freepts_extraporder2;
real *quasiset_trace_freepts_extraporder2;
real *quasiset_massdensity_freepts_extraporder2;
real *AdS_mass_freepts_extraporder2;
real *bdyphi_freepts_extraporder2;
real *xextrap_freepts_extraporder2;
real *yextrap_freepts_extraporder2;
real *zextrap_freepts_extraporder2;

real *quasiset_tt_fixedpts_extraporder1;
real *quasiset_tchi_fixedpts_extraporder1;
real *quasiset_txi_fixedpts_extraporder1;
real *quasiset_chichi_fixedpts_extraporder1;
real *quasiset_chixi_fixedpts_extraporder1;
real *quasiset_xixi_fixedpts_extraporder1;
real *quasiset_trace_fixedpts_extraporder1;
real *quasiset_massdensity_fixedpts_extraporder1;
real *AdS_mass_fixedpts_extraporder1;
real *bdyphi_fixedpts_extraporder1;
real *xextrap_fixedpts_extraporder1;
real *yextrap_fixedpts_extraporder1;
real *zextrap_fixedpts_extraporder1;

real *quasiset_tt_fixedpts_extraporder2;
real *quasiset_tchi_fixedpts_extraporder2;
real *quasiset_txi_fixedpts_extraporder2;
real *quasiset_chichi_fixedpts_extraporder2;
real *quasiset_chixi_fixedpts_extraporder2;
real *quasiset_xixi_fixedpts_extraporder2;
real *quasiset_trace_fixedpts_extraporder2;
real *quasiset_massdensity_fixedpts_extraporder2;
real *AdS_mass_fixedpts_extraporder2;
real *bdyphi_fixedpts_extraporder2;
real *xextrap_fixedpts_extraporder2;
real *yextrap_fixedpts_extraporder2;
real *zextrap_fixedpts_extraporder2;

real *relkretschcentregrid, *lrelkretschcentregrid0, *maxrelkretschcentregrid0,*minrelkretschcentregrid0,*relkretschcentregrid0;

real *tfunction,*test1,*test2,*test3,*test4;
real *iresall,*irestt,*irestx,*iresty,*iresxx,*iresxy,*iresyy,*irespsi;
real *irestz;
real *iresxz;
real *iresyz;
real *ireskg;

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

//variables that are defined in hyperbolic_vars but just because we want them to be synchronized (synchronization only affects the current time level)
int relkretsch_gfn,relkretsch_n_gfn,relkretsch_np1_gfn,relkretsch_nm1_gfn;

int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;

int chrbdy_freepts_extraporder1_gfn;
int chrbdy_freepts_extraporder2_gfn;
int chrbdy_fixedpts_extraporder1_gfn;
int chrbdy_fixedpts_extraporder2_gfn;

real *gu_tt,*gu_tx,*gu_ty,*gu_tz,*gu_xx,*gu_xy,*gu_xz,*gu_yy,*gu_yz,*gu_psi,*m_g_det;
int kg_ires_gfn,alpha_gfn,theta_gfn,f_gfn;

int phi1_res_gfn,gb_res_gfn;
int efe_all_ires_gfn;
int efe_tt_ires_gfn,efe_tx_ires_gfn,efe_ty_ires_gfn;
int efe_tz_ires_gfn;
int efe_xx_ires_gfn,efe_xy_ires_gfn,efe_yy_ires_gfn,efe_psi_ires_gfn;
int efe_xz_ires_gfn,efe_yz_ires_gfn;
int leadordcoeff_phi1_gfn;
int quasiset_tt_ll_gfn;
int quasiset_tchi_ll_gfn,quasiset_txi_ll_gfn;
int quasiset_chichi_ll_gfn,quasiset_chixi_ll_gfn,quasiset_xixi_ll_gfn;
int quasiset_tracell_gfn,quasiset_massdensityll_gfn;
int relkretschcentregrid_gfn;

int tfunction_gfn,test1_gfn,test2_gfn,test3_gfn,test4_gfn;
int iresall_gfn,irestt_gfn,irestx_gfn,iresty_gfn,irestz_gfn,iresxx_gfn,iresxy_gfn,iresxz_gfn,iresyy_gfn,iresyz_gfn,irespsi_gfn;
int ireskg_gfn;

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
real *AH_g0_tt[MAX_BHS];
real *AH_g0_tx[MAX_BHS];
real *AH_g0_ty[MAX_BHS];
real *AH_g0_tz[MAX_BHS];
real *AH_g0_xx[MAX_BHS];
real *AH_g0_xy[MAX_BHS];
real *AH_g0_xz[MAX_BHS];
real *AH_g0_yy[MAX_BHS];
real *AH_g0_yz[MAX_BHS];
real *AH_g0_psi[MAX_BHS];

real *AH_wtt1[MAX_BHS];
real *AH_wtx1[MAX_BHS];
real *AH_wty1[MAX_BHS];
real *AH_wtz1[MAX_BHS];
real *AH_wxx1[MAX_BHS];
real *AH_wxy1[MAX_BHS];
real *AH_wxz1[MAX_BHS];
real *AH_wyy1[MAX_BHS];
real *AH_wyz1[MAX_BHS];
real *AH_wpsi1[MAX_BHS];



real *AH_theta_ads[MAX_BHS];
int *AH_lev[MAX_BHS],*AH_own[MAX_BHS];

//=============================================================================
// arrays holding quasiset of CFT on ESU at outer bdy
//=============================================================================
real *lquasiset_tt0_freepts_extraporder1;
real *lquasiset_tchi0_freepts_extraporder1;
real *lquasiset_txi0_freepts_extraporder1;
real *lquasiset_chichi0_freepts_extraporder1;
real *lquasiset_chixi0_freepts_extraporder1;
real *lquasiset_xixi0_freepts_extraporder1;
real *lquasiset_trace0_freepts_extraporder1;
real *lquasiset_massdensity0_freepts_extraporder1;
real *lAdS_mass0_freepts_extraporder1;
real *lbdyphi0_freepts_extraporder1;

real *lquasiset_tt0_freepts_extraporder2;
real *lquasiset_tchi0_freepts_extraporder2;
real *lquasiset_txi0_freepts_extraporder2;
real *lquasiset_chichi0_freepts_extraporder2;
real *lquasiset_chixi0_freepts_extraporder2;
real *lquasiset_xixi0_freepts_extraporder2;
real *lquasiset_trace0_freepts_extraporder2;
real *lquasiset_massdensity0_freepts_extraporder2;
real *lAdS_mass0_freepts_extraporder2;
real *lbdyphi0_freepts_extraporder2;

real *lquasiset_tt0_fixedpts_extraporder1;
real *lquasiset_tchi0_fixedpts_extraporder1;
real *lquasiset_txi0_fixedpts_extraporder1;
real *lquasiset_chichi0_fixedpts_extraporder1;
real *lquasiset_chixi0_fixedpts_extraporder1;
real *lquasiset_xixi0_fixedpts_extraporder1;
real *lquasiset_trace0_fixedpts_extraporder1;
real *lquasiset_massdensity0_fixedpts_extraporder1;
real *lAdS_mass0_fixedpts_extraporder1;
real *lbdyphi0_fixedpts_extraporder1;

real *lquasiset_tt0_fixedpts_extraporder2;
real *lquasiset_tchi0_fixedpts_extraporder2;
real *lquasiset_txi0_fixedpts_extraporder2;
real *lquasiset_chichi0_fixedpts_extraporder2;
real *lquasiset_chixi0_fixedpts_extraporder2;
real *lquasiset_xixi0_fixedpts_extraporder2;
real *lquasiset_trace0_fixedpts_extraporder2;
real *lquasiset_massdensity0_fixedpts_extraporder2;
real *lAdS_mass0_fixedpts_extraporder2;
real *lbdyphi0_fixedpts_extraporder2;

real *maxquasiset_tt0_freepts_extraporder1;
real *maxquasiset_tchi0_freepts_extraporder1;
real *maxquasiset_txi0_freepts_extraporder1;
real *maxquasiset_chichi0_freepts_extraporder1;
real *maxquasiset_chixi0_freepts_extraporder1;
real *maxquasiset_xixi0_freepts_extraporder1;
real *maxquasiset_trace0_freepts_extraporder1;
real *maxquasiset_massdensity0_freepts_extraporder1;
real *maxbdyphi0_freepts_extraporder1;

real *maxquasiset_tt0_freepts_extraporder2;
real *maxquasiset_tchi0_freepts_extraporder2;
real *maxquasiset_txi0_freepts_extraporder2;
real *maxquasiset_chichi0_freepts_extraporder2;
real *maxquasiset_chixi0_freepts_extraporder2;
real *maxquasiset_xixi0_freepts_extraporder2;
real *maxquasiset_trace0_freepts_extraporder2;
real *maxquasiset_massdensity0_freepts_extraporder2;
real *maxbdyphi0_freepts_extraporder2;

real *maxquasiset_tt0_fixedpts_extraporder1;
real *maxquasiset_tchi0_fixedpts_extraporder1;
real *maxquasiset_txi0_fixedpts_extraporder1;
real *maxquasiset_chichi0_fixedpts_extraporder1;
real *maxquasiset_chixi0_fixedpts_extraporder1;
real *maxquasiset_xixi0_fixedpts_extraporder1;
real *maxquasiset_trace0_fixedpts_extraporder1;
real *maxquasiset_massdensity0_fixedpts_extraporder1;
real *maxbdyphi0_fixedpts_extraporder1;

real *maxquasiset_tt0_fixedpts_extraporder2;
real *maxquasiset_tchi0_fixedpts_extraporder2;
real *maxquasiset_txi0_fixedpts_extraporder2;
real *maxquasiset_chichi0_fixedpts_extraporder2;
real *maxquasiset_chixi0_fixedpts_extraporder2;
real *maxquasiset_xixi0_fixedpts_extraporder2;
real *maxquasiset_trace0_fixedpts_extraporder2;
real *maxquasiset_massdensity0_fixedpts_extraporder2;
real *maxbdyphi0_fixedpts_extraporder2;

real *minquasiset_tt0_freepts_extraporder1;
real *minquasiset_tchi0_freepts_extraporder1;
real *minquasiset_txi0_freepts_extraporder1;
real *minquasiset_chichi0_freepts_extraporder1;
real *minquasiset_chixi0_freepts_extraporder1;
real *minquasiset_xixi0_freepts_extraporder1;
real *minquasiset_trace0_freepts_extraporder1;
real *minquasiset_massdensity0_freepts_extraporder1;
real *minbdyphi0_freepts_extraporder1;

real *minquasiset_tt0_freepts_extraporder2;
real *minquasiset_tchi0_freepts_extraporder2;
real *minquasiset_txi0_freepts_extraporder2;
real *minquasiset_chichi0_freepts_extraporder2;
real *minquasiset_chixi0_freepts_extraporder2;
real *minquasiset_xixi0_freepts_extraporder2;
real *minquasiset_trace0_freepts_extraporder2;
real *minquasiset_massdensity0_freepts_extraporder2;
real *minbdyphi0_freepts_extraporder2;

real *minquasiset_tt0_fixedpts_extraporder1;
real *minquasiset_tchi0_fixedpts_extraporder1;
real *minquasiset_txi0_fixedpts_extraporder1;
real *minquasiset_chichi0_fixedpts_extraporder1;
real *minquasiset_chixi0_fixedpts_extraporder1;
real *minquasiset_xixi0_fixedpts_extraporder1;
real *minquasiset_trace0_fixedpts_extraporder1;
real *minquasiset_massdensity0_fixedpts_extraporder1;
real *minbdyphi0_fixedpts_extraporder1;

real *minquasiset_tt0_fixedpts_extraporder2;
real *minquasiset_tchi0_fixedpts_extraporder2;
real *minquasiset_txi0_fixedpts_extraporder2;
real *minquasiset_chichi0_fixedpts_extraporder2;
real *minquasiset_chixi0_fixedpts_extraporder2;
real *minquasiset_xixi0_fixedpts_extraporder2;
real *minquasiset_trace0_fixedpts_extraporder2;
real *minquasiset_massdensity0_fixedpts_extraporder2;
real *minbdyphi0_fixedpts_extraporder2;

real *quasiset_tt0_freepts_extraporder1;
real *quasiset_tchi0_freepts_extraporder1;
real *quasiset_txi0_freepts_extraporder1;
real *quasiset_chichi0_freepts_extraporder1;
real *quasiset_chixi0_freepts_extraporder1;
real *quasiset_xixi0_freepts_extraporder1;
real *quasiset_trace0_freepts_extraporder1;
real *quasiset_massdensity0_freepts_extraporder1;
real *AdS_mass0_freepts_extraporder1;
real *bdyphi0_freepts_extraporder1;

real *quasiset_tt0_freepts_extraporder2;
real *quasiset_tchi0_freepts_extraporder2;
real *quasiset_txi0_freepts_extraporder2;
real *quasiset_chichi0_freepts_extraporder2;
real *quasiset_chixi0_freepts_extraporder2;
real *quasiset_xixi0_freepts_extraporder2;
real *quasiset_trace0_freepts_extraporder2;
real *quasiset_massdensity0_freepts_extraporder2;
real *AdS_mass0_freepts_extraporder2;
real *bdyphi0_freepts_extraporder2;

real *quasiset_tt0_fixedpts_extraporder1;
real *quasiset_tchi0_fixedpts_extraporder1;
real *quasiset_txi0_fixedpts_extraporder1;
real *quasiset_chichi0_fixedpts_extraporder1;
real *quasiset_chixi0_fixedpts_extraporder1;
real *quasiset_xixi0_fixedpts_extraporder1;
real *quasiset_trace0_fixedpts_extraporder1;
real *quasiset_massdensity0_fixedpts_extraporder1;
real *AdS_mass0_fixedpts_extraporder1;
real *bdyphi0_fixedpts_extraporder1;

real *quasiset_tt0_fixedpts_extraporder2;
real *quasiset_tchi0_fixedpts_extraporder2;
real *quasiset_txi0_fixedpts_extraporder2;
real *quasiset_chichi0_fixedpts_extraporder2;
real *quasiset_chixi0_fixedpts_extraporder2;
real *quasiset_xixi0_fixedpts_extraporder2;
real *quasiset_trace0_fixedpts_extraporder2;
real *quasiset_massdensity0_fixedpts_extraporder2;
real *AdS_mass0_fixedpts_extraporder2;
real *bdyphi0_fixedpts_extraporder2;

real *xextrap0_freepts_extraporder1;
real *yextrap0_freepts_extraporder1;
real *zextrap0_freepts_extraporder1;

real *xextrap0_freepts_extraporder2;
real *yextrap0_freepts_extraporder2;
real *zextrap0_freepts_extraporder2;

real *xextrap0_fixedpts_extraporder1;
real *yextrap0_fixedpts_extraporder1;
real *zextrap0_fixedpts_extraporder1;

real *xextrap0_fixedpts_extraporder2;
real *yextrap0_fixedpts_extraporder2;
real *zextrap0_fixedpts_extraporder2;

real *rhoextrap0_freepts_extraporder1;
real *chiextrap0_freepts_extraporder1;
real *xiextrap0_freepts_extraporder1;
real *rhobdy0_freepts_extraporder1;
real *chibdy0_freepts_extraporder1;
real *xibdy0_freepts_extraporder1;

real *rhoextrap0_freepts_extraporder2;
real *chiextrap0_freepts_extraporder2;
real *xiextrap0_freepts_extraporder2;
real *rhobdy0_freepts_extraporder2;
real *chibdy0_freepts_extraporder2;
real *xibdy0_freepts_extraporder2;

real *rhoextrap0_fixedpts_extraporder1;
real *chiextrap0_fixedpts_extraporder1;
real *xiextrap0_fixedpts_extraporder1;
real *rhobdy0_fixedpts_extraporder1;
real *chibdy0_fixedpts_extraporder1;
real *xibdy0_fixedpts_extraporder1;

real *rhoextrap0_fixedpts_extraporder2;
real *chiextrap0_fixedpts_extraporder2;
real *xiextrap0_fixedpts_extraporder2;
real *rhobdy0_fixedpts_extraporder2;
real *chibdy0_fixedpts_extraporder2;
real *xibdy0_fixedpts_extraporder2;

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

//variables that are defined in hyperbolic_vars but just because we want them to be synchronized (synchronization only affects the current time level)
    if ((relkretsch_gfn      = PAMR_get_gfn("relkretsch",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((relkretsch_nm1_gfn   = PAMR_get_gfn("relkretsch",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((relkretsch_n_gfn   = PAMR_get_gfn("relkretsch",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((relkretsch_np1_gfn   = PAMR_get_gfn("relkretsch",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);


    if ((zeta_gfn     = PAMR_get_gfn("zeta",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_res_gfn = PAMR_get_gfn("zeta_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_lop_gfn = PAMR_get_gfn("zeta_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zeta_rhs_gfn = PAMR_get_gfn("zeta_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((chrbdy_freepts_extraporder1_gfn     = PAMR_get_gfn("chrbdy_freepts_extraporder1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chrbdy_freepts_extraporder2_gfn     = PAMR_get_gfn("chrbdy_freepts_extraporder2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chrbdy_fixedpts_extraporder1_gfn     = PAMR_get_gfn("chrbdy_fixedpts_extraporder1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chrbdy_fixedpts_extraporder2_gfn     = PAMR_get_gfn("chrbdy_fixedpts_extraporder2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

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
    if ((leadordcoeff_phi1_gfn    = PAMR_get_gfn("leadordcoeff_phi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tt_ll_gfn    = PAMR_get_gfn("quasiset_tt_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tchi_ll_gfn    = PAMR_get_gfn("quasiset_tchi_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_txi_ll_gfn    = PAMR_get_gfn("quasiset_txi_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_chichi_ll_gfn    = PAMR_get_gfn("quasiset_chichi_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_chixi_ll_gfn    = PAMR_get_gfn("quasiset_chixi_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_xixi_ll_gfn    = PAMR_get_gfn("quasiset_xixi_ll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_tracell_gfn    = PAMR_get_gfn("quasiset_tracell",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((quasiset_massdensityll_gfn    = PAMR_get_gfn("quasiset_massdensityll",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

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
    if ((ireskg_gfn  = PAMR_get_gfn("ireskg",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

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
    Hb_z_t_n  = gfs[Hb_z_t_n_gfn-1];    //variables that are defined in hyperbolic_vars but just because we want them to be synchronized (synchronization only affects the current time level)
    relkretsch     = gfs[relkretsch_gfn-1];
    relkretsch_n   = gfs[relkretsch_n_gfn-1];
    relkretsch_nm1   = gfs[relkretsch_nm1_gfn-1];
    relkretsch_np1   = gfs[relkretsch_np1_gfn-1];   
    zeta     = gfs[zeta_gfn-1];
    zeta_lop = gfs[zeta_lop_gfn-1];
    zeta_res = gfs[zeta_res_gfn-1];
    zeta_rhs = gfs[zeta_rhs_gfn-1]; 
    mask    = gfs[mask_gfn-1];
    mask_mg = gfs[mask_mg_gfn-1];
    chr = gfs[chr_gfn-1]; 
    chr_mg = gfs[chr_mg_gfn-1]; 
    chrbdy_freepts_extraporder1 = gfs[chrbdy_freepts_extraporder1_gfn-1];
    chrbdy_freepts_extraporder2 = gfs[chrbdy_freepts_extraporder2_gfn-1];
    chrbdy_fixedpts_extraporder1 = gfs[chrbdy_fixedpts_extraporder1_gfn-1];
    chrbdy_fixedpts_extraporder2 = gfs[chrbdy_fixedpts_extraporder2_gfn-1]; 
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
    leadordcoeff_phi1  = gfs[leadordcoeff_phi1_gfn-1];
    quasiset_tt_ll  = gfs[quasiset_tt_ll_gfn-1];
    quasiset_tchi_ll  = gfs[quasiset_tchi_ll_gfn-1];
    quasiset_txi_ll  = gfs[quasiset_txi_ll_gfn-1];
    quasiset_chichi_ll  = gfs[quasiset_chichi_ll_gfn-1];
    quasiset_chixi_ll  = gfs[quasiset_chixi_ll_gfn-1];
    quasiset_xixi_ll  = gfs[quasiset_xixi_ll_gfn-1];
    quasiset_tracell  = gfs[quasiset_tracell_gfn-1];
    quasiset_massdensityll  = gfs[quasiset_massdensityll_gfn-1];
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
    ireskg = gfs[ireskg_gfn-1];
 
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
                is_inside_(&is_2in1,&intr,ex_r[i2],ex_xc[i2],ex_r[i1],ex_xc[i1],&AMRD_dim); //(unless this is disabled, AH[4] is removed at first time step)
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
    ratio_Lhighres_Llowres=1.5; AMRD_real_param(pfile,"ratio_Lhighres_Llowres",&ratio_Lhighres_Llowres,1);
    resolution_degree=1; AMRD_int_param(pfile,"resolution_degree",&resolution_degree,1);
    max_resolution_degree=3; AMRD_int_param(pfile,"max_resolution_degree",&max_resolution_degree,1);
    reduction_factor=2; AMRD_int_param(pfile,"reduction_factor",&reduction_factor,1);   
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
    output_kretsch=0; AMRD_int_param(pfile,"output_kretsch",&output_kretsch,1); 
    output_metricatAH=0; AMRD_int_param(pfile,"output_metricatAH",&output_metricatAH,1);    
    output_bdyquantities=0; AMRD_int_param(pfile,"output_bdyquantities",&output_bdyquantities,1);
    output_AdS_mass=0; AMRD_int_param(pfile,"output_AdS_mass",&output_AdS_mass,1);
    reduced_ascii=0; AMRD_int_param(pfile,"reduced_ascii",&reduced_ascii,1);
    reduction_factor_ascii=1; AMRD_int_param(pfile,"reduction_factor_ascii",&reduction_factor_ascii,1);
    alltimes_ascii=0; AMRD_int_param(pfile,"alltimes_ascii",&alltimes_ascii,1);
    timestep_ascii=0; AMRD_int_param(pfile,"timestep_ascii",&timestep_ascii,1);
    bdy_extrap_freepts=0; AMRD_int_param(pfile,"bdy_extrap_freepts",&bdy_extrap_freepts,1);
    bdy_extrap_fixedpts=0; AMRD_int_param(pfile,"bdy_extrap_fixedpts",&bdy_extrap_fixedpts,1);
    output_bdy_extraporder1=0; AMRD_int_param(pfile,"output_bdy_extraporder1",&output_bdy_extraporder1,1);
    output_bdy_extraporder2=0; AMRD_int_param(pfile,"output_bdy_extraporder2",&output_bdy_extraporder2,1);
    output_bdy_extraporder3=0; AMRD_int_param(pfile,"output_bdy_extraporder3",&output_bdy_extraporder3,1);
    output_bdy_extraporder4=0; AMRD_int_param(pfile,"output_bdy_extraporder4",&output_bdy_extraporder4,1);  

    //allocate memory for relative Kretschmann scalar at the centre of the grid
    if (output_relkretschcentregrid)
    {
        lrelkretschcentregrid0= malloc(sizeof(real));
        maxrelkretschcentregrid0= malloc(sizeof(real));
        minrelkretschcentregrid0= malloc(sizeof(real));
        relkretschcentregrid0= malloc(sizeof(real));
    }   
    if (output_bdyquantities)
    {
    //allocate memory for fixed_coords, i.e. the values (fixed for all resolutions) of coordinates of points that we use for use boundary extrapolation
        if (bdy_extrap_fixedpts)
        {
            currentres_reduction_factor = round(pow(reduction_factor,max_resolution_degree-1));
            currentres_ratio_Lhighres_Llowres=pow(ratio_Lhighres_Llowres,resolution_degree-1);
            ind_distance_fixedpts=round(currentres_reduction_factor*currentres_ratio_Lhighres_Llowres);
            num_fixed_coords=round(((AMRD_base_shape[0]-1)/(ind_distance_fixedpts)+1));
            fixed_coords= malloc(num_fixed_coords*sizeof(real));    
        }
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
        if (output_metricatAH)
        {
            AH_g0_tt[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_tx[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_ty[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_tz[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_xx[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_xy[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_xz[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_yy[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_yz[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_g0_psi[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
        }   
        if (!AMRD_cp_restart || !found_AH[l]) AH_R[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
        AH_w1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));    
        if (output_metricatAH)
        {
            AH_wtt1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wtx1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wty1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wtz1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wxx1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wxy1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wxz1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wyy1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wyz1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
            AH_wpsi1[l]=(real *)malloc(AH_Nchi[l]*AH_Nphi[l]*sizeof(real));
        }   
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
        {   //       ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh;
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
    {   //   for (i=0; i<Nx; i++)  //   {    //      for (j=0; j<Ny; j++)    //      {   //       for (k=0; k<Nz; k++)  //       {    //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)    //        { //         ind=i+Nx*(j+Ny*k);    //         printf("AdS4D_AMRH_var_clear-PRE init_ghbdot_\n");   //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"  //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));    //         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);   //         printf("gb_xx_nm1[ind]=%lf,gb_xx_n[ind]=%lf,gb_xx_np1[ind]=%lf\n",gb_xx_nm1[ind],gb_xx_n[ind],gb_xx_np1[ind]);  //        }   //       } //      }    //   }  
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
                    &AdS_L,phys_bdy,x,y,z,&dt,chr,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);   //   for (i=0; i<Nx; i++)  //   {    //      for (j=0; j<Ny; j++)    //      {   //       for (k=0; k<Nz; k++)  //       {    //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)    //        { //         ind=i+Nx*(j+Ny*k);    //         printf("AdS4D_AMRH_var_clear-POST init_ghbdot_\n");  //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n" //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));   //         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);  //        }   //       } //      }    //   }  
    }   
    // initialize gbars and nm1, np1 time levels
    if ((background || skip_constraints) && ief_bh_r0==0)
    {   //   for (i=0; i<Nx; i++)  //   {    //      for (j=0; j<Ny; j++)    //      {   //       for (k=0; k<Nz; k++)  //       {    //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)    //        { //         ind=i+Nx*(j+Ny*k);    //         printf("AdS4D_AMRH_var_clear-PRE init_ghb_ads_\n");  //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n" //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));   //         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);  //        }   //       } //      }    //   }  
        init_ghb_ads_(gb_tt,gb_tx,gb_ty,
                    gb_tz,
                    gb_xx,gb_xy,
                    gb_xz,
                    gb_yy,
                    gb_yz,
                    psi,Hb_t,Hb_x,Hb_y,
                    Hb_z,
                    &AdS_L,x,y,z,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype); //   for (i=0; i<Nx; i++)    //   {  //      for (j=0; j<Ny; j++)  //      { //       for (k=0; k<Nz; k++)    //       {  //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)  //        {   //         ind=i+Nx*(j+Ny*k);  //         printf("AdS4D_AMRH_var_clear-POST init_ghb_ads_\n");   //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"  //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));    //         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);   //        }    //       }  //      } //   }   
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
                            phys_bdy,x,y,z,&dt,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype);   //   for (i=0; i<Nx; i++)  //   {    //      for (j=0; j<Ny; j++)    //      {   //       for (k=0; k<Nz; k++)  //       {    //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)    //        { //         ind=i+Nx*(j+Ny*k);    //         printf("AdS4D_AMRH_var_clear-POST init_ads4d_bh_\n");    //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n"   //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])); //         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);    //        } //       }   //      }  //   }    
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
        if (gb_xx_nm1) 
        //"np1,n,nm1" variables only allocated on finest MG level //IT WAS NOT HERE IN PREVIOUS VERSIONS, ASK HANS
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
                    &AdS_L,mask_mg,phys_bdy,x,y,z,chr_mg,&AMRD_ex,&Nx,&Ny,&Nz,&regtype,&rhoa,&rhob);    //   for (i=0; i<Nx; i++)   //   { //      for (j=0; j<Ny; j++) //      {    //       for (k=0; k<Nz; k++)   //       { //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1) //        {  //         ind=i+Nx*(j+Ny*k); //         printf("AdS4D_AMRH_var_clear-POST init_ghb_\n");  //         printf("i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i,x[i]=%lf,y[j]=%lf,z[k]=%lf,rho=%lf\n" //                ,i,j,k,Nx,Ny,Nz,x[i],y[j],z[k],sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]));   ////         printf("gb_tt_nm1[ind]=%lf,gb_tt_n[ind]=%lf,gb_tt_np1[ind]=%lf\n",gb_tt_nm1[ind],gb_tt_n[ind],gb_tt_np1[ind]);    //         printf("gb_xx[ind]=%lf,gb_xx_nm1[ind]=%lf,gb_xx_n[ind]=%lf,gb_xx_np1[ind]=%lf\n",gb_xx[ind],gb_xx_nm1[ind],gb_xx_n[ind],gb_xx_np1[ind]); //        }  //       }    //      }   //   } 
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
                    //   {    //      for (j=0; j<Ny; j++)
                    //      { //       for (k=0; k<Nz; k++)
                    //       {  //        if (sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k])<1)
                    //        {   //         ind=i+Nx*(j+Ny*k);
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
    //compute relative Kretschmann scalar of initial data
    if (gb_xx_nm1)
    {
        if (output_kretsch)
        {
            kretsch_(relkretsch_n,
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
            //NOTICE: relkretsch_np1 is not synchronized yet at this stage, meaning that only 1 process has a non-zero value at x=y=z=0. This is crucial for how we save and print relkretschcentregrid in pre_tstep    
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

    //printf("AdS4D_pre_io_calc is called,lsteps=%i\n",lsteps);
    //fflush(stdout); 

    // compute independent residuals of the AdS4D system
    if (ct!=0)
    {
        //(NOTE: for t>t0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
        // so here, time level n is the most advanced time level)   
        // output independent residual
        if (output_ires)
        {   
            ires_(efe_all_ires,
                efe_tt_ires,efe_tx_ires,efe_ty_ires,
                efe_tz_ires,
                efe_xx_ires,efe_xy_ires,
                efe_xz_ires,
                efe_yy_ires,
                efe_yz_ires,
                efe_psi_ires,
                kg_ires,
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
        //we call the following functions just to save and output the values of the grid functions chrbdy, leadordcoeff_phi1 and quasiset_ll for t=0. These will be later then used to extrapolate the value of bdyphi and quasiset componenents at the AdS boundary in pre_tstep for t=0 (and post_tstep for later times)
        if (output_bdyquantities)
        {   
            calc_leadordcoeff_phi1_(leadordcoeff_phi1,
                                    phi1_np1,phi1_n,phi1_nm1,
                                    x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
            calc_quasiset_ll_(
                        quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                        quasiset_chichi_ll,quasiset_chixi_ll,
                        quasiset_xixi_ll,
                        quasiset_tracell,
                        quasiset_massdensityll,
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
                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
        	if (bdy_extrap_freepts)
        	{   
            	if (output_bdy_extraporder1)
            	{
                	bdy_extrap_order=1; //routine that sets a mask for near bdy points. We will call these "nexttobdypoints". The number of nexttobdypoints is also the number of points at the boundary where we will extrapolate the stress-energy tensor in AdS4D_pre_tstep and AdS4D_post_tstep. We call this number numbdypoints.
                	nexttobdypoints_freepts_(chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,&bdy_extrap_order,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
            	}   
            	if (output_bdy_extraporder2)
            	{
                	bdy_extrap_order=2;
                	nexttobdypoints_freepts_(chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,&bdy_extrap_order,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
            	}   
        	}   
        	if (bdy_extrap_fixedpts)
        	{   
            	//set fixed_coords, i.e. the values (fixed for all resolutions) of coordinates of points that we use for use boundary extrapolation
            	for (i=0;i<num_fixed_coords;i++)
            	{
                	fixed_coords[i]=x[i*ind_distance_fixedpts]; //         printf("i=%i,fixed_coords[i]=%lf\n",i,fixed_coords[i]);
            	}           
            	if (output_bdy_extraporder1)
            	{
                	bdy_extrap_order=1;
                	nexttobdypoints_fixedpts_(chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,&bdy_extrap_order,&ind_distance_fixedpts,&currentres_ratio_Lhighres_Llowres,&num_fixed_coords,fixed_coords,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
            	}   
            	if (output_bdy_extraporder2)
            	{
                	bdy_extrap_order=2;
                	nexttobdypoints_fixedpts_(chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,&bdy_extrap_order,&ind_distance_fixedpts,&currentres_ratio_Lhighres_Llowres,&num_fixed_coords,fixed_coords,x,y,z,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);
            	}   
        	}   
            //reduce boundary points if needed  
            //         int reduced_numbdypoints=10000;    
            //         if (numbdypoints*uniSize>reduced_numbdypoints) numbdypoints=roundl(reduced_numbdypoints/uniSize)+1;  
        } //closes condition on output_bdyquantities    
        if (output_ires)
        {   
            ires_(efe_all_ires,
                efe_tt_ires,efe_tx_ires,efe_ty_ires,
                efe_tz_ires,
                efe_xx_ires,efe_xy_ires,
                efe_xz_ires,
                efe_yy_ires,
                efe_yz_ires,
                efe_psi_ires,
                kg_ires,
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
        }   
    } //closes condition on ct==0   

    if (output_ires)
    {
        // fill in independent residual evaluator test functions
        for (i=0; i<Nx; i++)
        {
            for (j=0; j<Ny; j++)
            {
                for (k=0; k<Nz; k++)
                {
                    ind=i+Nx*(j+Ny*k);
                    rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);    
                    if (chr[ind]==AMRD_ex)
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
                        ireskg[ind]=0;
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
                        ireskg[ind]=kg_ires[ind];
                    }
                }
            } 
        }
    }  //closes (output_ires)


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
    int i,j,k,m,zero_i=0;
    int is,ie,js,je,ks,ke;
    int ltrace=0;
    real ct,zero=0;
    int ind;
    real rho;   

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
                //     MPI_Comm_size(MPI_COMM_WORLD,&uniSize);  
                //     for (m=0; m<uniSize; m++)  
                //     {  
                //      for (i=0; i<Nx; i++)  
                //      { 
                //         for (j=0; j<Ny; j++)  
                //         {  
                //          for (k=Nz-2; k<Nz; k++)   
                //          {  
                //            ind=i+Nx*(j+Ny*k);  
                //            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);    
                //            if ((rho<1)&&(my_rank==m))    
                //            { 
                //             printf("BEFORE g_evo_opt: my_rank=%i,i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i\n,ct=%lf,x=%lf,y=%lf,z=%lf,gb_xx_n[ind]=%lf,gb_xx_nm1[ind]=%lf,gb_xx_np1[ind]=%lf\n",my_rank,i,j,k,Nx,Ny,Nz,ct,x[i],y[j],z[k],gb_xx_n[ind],gb_xx_nm1[ind],gb_xx_np1[ind]);  
                //            }   
                //          }  
                //         }  
                //      } 
                //     } 
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
                //     MPI_Comm_size(MPI_COMM_WORLD,&uniSize);    
                //     for (m=0; m<uniSize; m++)    
                //     {    
                //      for (i=0; i<Nx; i++)    
                //      {   
                //         for (j=0; j<Ny; j++)    
                //         {    
                //          for (k=Nz-2; k<Nz; k++) 
                //          {    
                //            ind=i+Nx*(j+Ny*k);    
                //            rho=sqrt(x[i]*x[i]+y[j]*y[j]+z[k]*z[k]);  
                //            if ((rho<1)&&(my_rank==m))  
                //            {   
                //             printf("AFTER g_evo_opt: my_rank=%i,i=%i,j=%i,k=%i,Nx=%i,Ny=%i,Nz=%i\n,ct=%lf,x=%lf,y=%lf,z=%lf,gb_xx_n[ind]=%lf,gb_xx_nm1[ind]=%lf,gb_xx_np1[ind]=%lf\n",my_rank,i,j,k,Nx,Ny,Nz,ct,x[i],y[j],z[k],gb_xx_n[ind],gb_xx_nm1[ind],gb_xx_np1[ind]); 
                //            }  
                //          } 
                //         } 
                //      }    
                //     }    
                //
                //TO CHECK THE ORDER OF SYNCHRONIZATION BETWEEN DIFFERENT GRID PROCESSES AND CYCLING FROM np1,n,nm1 to n,nm1,np1 (done at the end of app_evolve and before app_pre_io_calc), AND TO CHECK THAT SYNCHRONIZATION ONLY ACTS ON TIME LEVEL 2 OF THE AMR HIERARCHY, I.E. ON gb_tt_n,gb_tx_n,etc., not on gb_tt_nm1,gb_tt_np1,gb_tx_nm1,gb_tx_np1,etc.  
                //The order is: 1. Synchronization of the np1 time level, 2.cycling of the 3 time levels. 
                //         is=1; 
                //         ie=Nx-1;  
                //         js=1;  
                //         je=Ny-1;   
                //         ks=1;   
                //         ke=Nz-1;    
                //  
                //        if (ghost_width[0]>0) is=is+ghost_width[0]-1;   
                //        if (ghost_width[1]>0) ie=ie-(ghost_width[1]-1);  
                //        if (ghost_width[2]>0) js=js+ghost_width[2]-1;   
                //        if (ghost_width[3]>0) je=je-(ghost_width[3]-1);  
                //        if (ghost_width[4]>0) ks=ks+ghost_width[4]-1;   
                //        if (ghost_width[5]>0) ke=ke-(ghost_width[5]-1);  
                //    
                //
                //At the end of this for-cycle we have non-synchronized gb_xx_n, gb_xx_nm1,gb_xx_np1, so that we are able to see which variables are synchronised and how they are cycled before app_pre_io_calc. We can see this by outputting quantities in pre_io_calc or looking at what is saved to disk immediately after pre_io_calc    
                //      for (i=0; i<Nx; i++)    
                //      {   
                //         for (j=0; j<Ny; j++)    
                //         {    
                //          for (k=0; k<Nz; k++)    
                //          {   
                //            ind=i+Nx*(j+Ny*k);   
                //            if ((i<is)||(i>=ie)||(j<js)||(j>=je)||(k<ks)||(k>=ke))   
                //            {    
                //             gb_xx_np1[ind]=0; //goes into gb_xx_n before pre_io_calc 
                //             gb_xx_n[ind]=0;   //goes into gb_xx_nm1 before pre_io_calc    
                //             gb_xx_nm1[ind]=0; //goes into gb_nn_np1 before pre_io_calc   
                //            }    
                //          }   
                //         }   
                //       } 
                //
                //
                //   
    //compute relative Kretschmann scalar after each iteration of evolution equations (after AdS4D_evolve variables in amr_sync are synchronized, so this is where we should compute variables that we want to be synchronized) 
    //NOTICE: before calling pre_io_calc, the function called before saving info to disk, relkretsch_np1 is synchronised and then moved to relkretsch_n (while relkretsch_n is moved to relkretsch_nm1 and relkretsch_nm1 is moved to relkretsch_np1). So we have to save relkretsch for the current time step into relkretsch_np1 to have it synchronised and saved in relkretsch_n
        if (output_kretsch)
        {       
            kretsch_(relkretsch_np1,
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
                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    //NOTICE: relkretsch_np1 is not synchronized yet at this stage, meaning that only 1 process has a non-zero value at x=y=z=0. This is crucial for how we save and print relkretschcentregrid in post_tstep   
        }   
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
                    mask[ind]=excised-1;    //       printf("point not excised");   //       printf("EX_MASK: i=%i,j=%i,k=%i,ind=%i\n",i,j,k,ind); //       printf("EX_MASK: x=%lf,y=%lf,z=%lf,rho=%lf\n",x,y,z,rho);   //       printf("EX_MASK: mask[ind]=%lf\n",mask[ind]);
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
                            ex_r_zp=(ex_r[l][2]*(1-ex_rbuf[l]));    //                 printf("xp=%lf,yp=%lf,zp=%lf\n",xp,yp,zp);   //                 printf("ex_r_xp=%lf,ex_r_yp=%lf,ex_r_zp=%lf\n",ex_r_xp,ex_r_yp,ex_r_zp);    
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

//int mem_alloc_bdyquantities_first=1;
int is_bdy_freepts_extraporder1,ie_bdy_freepts_extraporder1;
int is_bdy_freepts_extraporder2,ie_bdy_freepts_extraporder2;
int is_bdy_fixedpts_extraporder1,ie_bdy_fixedpts_extraporder1;
int is_bdy_fixedpts_extraporder2,ie_bdy_fixedpts_extraporder2;

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
    int count_relkretschcentregrid;
    real rho;
    real rh,mh,rhoh;    
    ct=PAMR_get_time(L);    
    Lf=PAMR_get_max_lev(PAMR_AMRH);
    Lc=PAMR_get_min_lev(PAMR_AMRH);  

    //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1; 

    //printf("AdS4D_pre_tstep is called");
    //fflush(stdout); 

    if (AMRD_state!=AMRD_STATE_EVOLVE) return; // if disable, enable(?) reset_AH_shapes below   
    if (pre_tstep_global_first)
    {
        for (l=0; l<MAX_BHS; l++) { AH_count[l]=found_AH[l]=found_count_AH[l]=0; freq0[l]=AH_freq[l]; }
        pre_tstep_global_first=0;
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
            //we save here the values of the Kretschmann scalar at the centre of the grid at times greater than 0, computed at the end of g_evo_opt //The value of the Kretschmann scalar at the centre of the grid at t=0 is saved in pre_io_calc.  
            if (output_relkretschcentregrid)
            {
                *lrelkretschcentregrid0= *relkretschcentregrid;
            }   
            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize); 
                        vecbdypoints_freepts_extraporder1 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_freepts_extraporder1 = malloc(uniSize*sizeof(int));    
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_freepts_extraporder1,1,MPI_INT,vecbdypoints_freepts_extraporder1,1,MPI_INT,MPI_COMM_WORLD); 
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
                        //printf("numbdypoints_freepts_extraporder1=%i\n",numbdypoints_freepts_extraporder1);
                        //fflush(stdout); 
                        quasiset_tt_freepts_extraporder1              = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_tchi_freepts_extraporder1            = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_txi_freepts_extraporder1             = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chichi_freepts_extraporder1          = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chixi_freepts_extraporder1           = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_xixi_freepts_extraporder1            = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_trace_freepts_extraporder1           = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_massdensity_freepts_extraporder1     = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        bdyphi_freepts_extraporder1                   = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));   
                        xextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        yextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        zextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));   
                        lquasiset_tt0_freepts_extraporder1            = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_tchi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_txi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_chichi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_chixi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_xixi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_trace0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_massdensity0_freepts_extraporder1   = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lbdyphi0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        maxquasiset_tt0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_tchi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_txi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_chichi0_freepts_extraporder1      = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_chixi0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_xixi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_trace0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_massdensity0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxbdyphi0_freepts_extraporder1               = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        minquasiset_tt0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_tchi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_txi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_chichi0_freepts_extraporder1      = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_chixi0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_xixi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_trace0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_massdensity0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minbdyphi0_freepts_extraporder1               = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        quasiset_tt0_freepts_extraporder1             = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_tchi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_txi0_freepts_extraporder1            = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chichi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chixi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_xixi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_trace0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_massdensity0_freepts_extraporder1    = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        AdS_mass0_freepts_extraporder1                = malloc(sizeof(real));
                        bdyphi0_freepts_extraporder1                  = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));    
                        xextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        yextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        zextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   

                        //initialize
                        for (i=0;i<numbdypoints_freepts_extraporder1;i++)
                        {
                            quasiset_tt_freepts_extraporder1             [i] = 0;
                            quasiset_tchi_freepts_extraporder1           [i] = 0;
                            quasiset_txi_freepts_extraporder1            [i] = 0;
                            quasiset_chichi_freepts_extraporder1         [i] = 0;
                            quasiset_chixi_freepts_extraporder1          [i] = 0;
                            quasiset_xixi_freepts_extraporder1           [i] = 0;
                            quasiset_trace_freepts_extraporder1          [i] = 0;
                            quasiset_massdensity_freepts_extraporder1    [i] = 0;
                            bdyphi_freepts_extraporder1                  [i] = 0;   
                            xextrap_freepts_extraporder1                 [i] = 0;
                            yextrap_freepts_extraporder1                 [i] = 0;
                            zextrap_freepts_extraporder1                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_freepts_extraporder1;i++)
                        {   
                            lquasiset_tt0_freepts_extraporder1            [i] = 0;
                            lquasiset_tchi0_freepts_extraporder1          [i] = 0;
                            lquasiset_txi0_freepts_extraporder1           [i] = 0;
                            lquasiset_chichi0_freepts_extraporder1        [i] = 0;
                            lquasiset_chixi0_freepts_extraporder1         [i] = 0;
                            lquasiset_xixi0_freepts_extraporder1          [i] = 0;
                            lquasiset_trace0_freepts_extraporder1         [i] = 0;
                            lquasiset_massdensity0_freepts_extraporder1   [i] = 0;
                            lbdyphi0_freepts_extraporder1                 [i] = 0;             
                            maxquasiset_tt0_freepts_extraporder1          [i] = 0;
                            maxquasiset_tchi0_freepts_extraporder1        [i] = 0;
                            maxquasiset_txi0_freepts_extraporder1         [i] = 0;
                            maxquasiset_chichi0_freepts_extraporder1      [i] = 0;
                            maxquasiset_chixi0_freepts_extraporder1       [i] = 0;
                            maxquasiset_xixi0_freepts_extraporder1        [i] = 0;
                            maxquasiset_trace0_freepts_extraporder1       [i] = 0;
                            maxquasiset_massdensity0_freepts_extraporder1 [i] = 0;
                            maxbdyphi0_freepts_extraporder1               [i] = 0;  
                            minquasiset_tt0_freepts_extraporder1          [i] = 0;
                            minquasiset_tchi0_freepts_extraporder1        [i] = 0;
                            minquasiset_txi0_freepts_extraporder1         [i] = 0;
                            minquasiset_chichi0_freepts_extraporder1      [i] = 0;
                            minquasiset_chixi0_freepts_extraporder1       [i] = 0;
                            minquasiset_xixi0_freepts_extraporder1        [i] = 0;
                            minquasiset_trace0_freepts_extraporder1       [i] = 0;
                            minquasiset_massdensity0_freepts_extraporder1 [i] = 0;
                            minbdyphi0_freepts_extraporder1               [i] = 0;  
                            quasiset_tt0_freepts_extraporder1             [i] = 0;
                            quasiset_tchi0_freepts_extraporder1           [i] = 0;
                            quasiset_txi0_freepts_extraporder1            [i] = 0;
                            quasiset_chichi0_freepts_extraporder1         [i] = 0;
                            quasiset_chixi0_freepts_extraporder1          [i] = 0;
                            quasiset_xixi0_freepts_extraporder1           [i] = 0;
                            quasiset_trace0_freepts_extraporder1          [i] = 0;
                            quasiset_massdensity0_freepts_extraporder1    [i] = 0;
                            bdyphi0_freepts_extraporder1                  [i] = 0;   
                            xextrap0_freepts_extraporder1                 [i] = 0;
                            yextrap0_freepts_extraporder1                 [i] = 0;
                            zextrap0_freepts_extraporder1                 [i] = 0;  
                        }
                        *AdS_mass0_freepts_extraporder1                    = 0; 
                        
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_freepts_extraporder1=0;
                        if (my_rank==0)
                        {
                            ie_bdy_freepts_extraporder1=vecbdypoints_freepts_extraporder1[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_freepts_extraporder1=is_bdy_freepts_extraporder1+vecbdypoints_freepts_extraporder1[j];
                            }
                            ie_bdy_freepts_extraporder1=is_bdy_freepts_extraporder1+vecbdypoints_freepts_extraporder1[my_rank];
                        }      
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_freepts_extraporder1[i]=0;
                        }         
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_freepts_extraporder1[i]=dsplsbdypoints_freepts_extraporder1[i]+vecbdypoints_freepts_extraporder1[j];
                                }
                            }
                        }   
                        
                        xyzextrap_(xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);    
                        
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,xextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,yextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,zextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);    //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0...");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FREE POINTS, FIRST ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_freepts_extraporder1 = malloc(sizeof(real));
                                chiextrap0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                                xiextrap0_freepts_extraporder1  = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real)); 
                                chixiextrap_(rhoextrap0_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1);  
                                basebdy_Nchi_freepts_extraporder1=0;//initialize
                                basebdy_Nxi_freepts_extraporder1=0; //initialize
                                bdyn_(&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1);  
                                rhobdy0_freepts_extraporder1 = malloc(sizeof(real));
                                chibdy0_freepts_extraporder1 = malloc(basebdy_Nchi_freepts_extraporder1*sizeof(real));
                                xibdy0_freepts_extraporder1  = malloc(basebdy_Nxi_freepts_extraporder1*sizeof(real));   
                            }
                        }

   
                        extrap_bdyphi_freepts_(bdyphi_freepts_extraporder1,
                                        leadordcoeff_phi1,
                                        xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,
                                        chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,
                                        &bdy_extrap_order,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_freepts_(quasiset_tt_freepts_extraporder1,quasiset_tchi_freepts_extraporder1,quasiset_txi_freepts_extraporder1,
                                quasiset_chichi_freepts_extraporder1,quasiset_chixi_freepts_extraporder1,
                                quasiset_xixi_freepts_extraporder1,
                                quasiset_trace_freepts_extraporder1,
                                quasiset_massdensity_freepts_extraporder1,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,
                                chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,
                                &bdy_extrap_order,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_freepts_extraporder1; i<ie_bdy_freepts_extraporder1; i++)
                        {
                            lquasiset_tt0_freepts_extraporder1           [i] = quasiset_tt_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_tchi0_freepts_extraporder1         [i] = quasiset_tchi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_txi0_freepts_extraporder1          [i] = quasiset_txi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_chichi0_freepts_extraporder1       [i] = quasiset_chichi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_chixi0_freepts_extraporder1        [i] = quasiset_chixi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_xixi0_freepts_extraporder1         [i] = quasiset_xixi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_trace0_freepts_extraporder1        [i] = quasiset_trace_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_massdensity0_freepts_extraporder1  [i] = quasiset_massdensity_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lbdyphi0_freepts_extraporder1                [i] = bdyphi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_freepts_extraporder2 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_freepts_extraporder2 = malloc(uniSize*sizeof(int));  
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_freepts_extraporder2,1,MPI_INT,vecbdypoints_freepts_extraporder2,1,MPI_INT,MPI_COMM_WORLD); 
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
                        quasiset_tt_freepts_extraporder2              = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_tchi_freepts_extraporder2            = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_txi_freepts_extraporder2             = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chichi_freepts_extraporder2          = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chixi_freepts_extraporder2           = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_xixi_freepts_extraporder2            = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_trace_freepts_extraporder2           = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_massdensity_freepts_extraporder2     = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        bdyphi_freepts_extraporder2                   = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));   
                        xextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        yextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        zextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));   
                        lquasiset_tt0_freepts_extraporder2            = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_tchi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_txi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_chichi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_chixi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_xixi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_trace0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_massdensity0_freepts_extraporder2   = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lbdyphi0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        maxquasiset_tt0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_tchi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_txi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_chichi0_freepts_extraporder2      = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_chixi0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_xixi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_trace0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_massdensity0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxbdyphi0_freepts_extraporder2               = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        minquasiset_tt0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_tchi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_txi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_chichi0_freepts_extraporder2      = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_chixi0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_xixi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_trace0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_massdensity0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minbdyphi0_freepts_extraporder2               = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        quasiset_tt0_freepts_extraporder2             = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_tchi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_txi0_freepts_extraporder2            = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chichi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chixi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_xixi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_trace0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_massdensity0_freepts_extraporder2    = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        AdS_mass0_freepts_extraporder2                = malloc(sizeof(real));
                        bdyphi0_freepts_extraporder2                  = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        xextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        yextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        zextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        //initialize
                        for (i=0;i<numbdypoints_freepts_extraporder2;i++)
                        {
                            quasiset_tt_freepts_extraporder2             [i] = 0;
                            quasiset_tchi_freepts_extraporder2           [i] = 0;
                            quasiset_txi_freepts_extraporder2            [i] = 0;
                            quasiset_chichi_freepts_extraporder2         [i] = 0;
                            quasiset_chixi_freepts_extraporder2          [i] = 0;
                            quasiset_xixi_freepts_extraporder2           [i] = 0;
                            quasiset_trace_freepts_extraporder2          [i] = 0;
                            quasiset_massdensity_freepts_extraporder2    [i] = 0;
                            bdyphi_freepts_extraporder2                  [i] = 0;   
                            xextrap_freepts_extraporder2                 [i] = 0;
                            yextrap_freepts_extraporder2                 [i] = 0;
                            zextrap_freepts_extraporder2                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_freepts_extraporder2;i++)
                        {   
                            lquasiset_tt0_freepts_extraporder2            [i] = 0;
                            lquasiset_tchi0_freepts_extraporder2          [i] = 0;
                            lquasiset_txi0_freepts_extraporder2           [i] = 0;
                            lquasiset_chichi0_freepts_extraporder2        [i] = 0;
                            lquasiset_chixi0_freepts_extraporder2         [i] = 0;
                            lquasiset_xixi0_freepts_extraporder2          [i] = 0;
                            lquasiset_trace0_freepts_extraporder2         [i] = 0;
                            lquasiset_massdensity0_freepts_extraporder2   [i] = 0;
                            lbdyphi0_freepts_extraporder2                 [i] = 0;  
                            maxquasiset_tt0_freepts_extraporder2          [i] = 0;
                            maxquasiset_tchi0_freepts_extraporder2        [i] = 0;
                            maxquasiset_txi0_freepts_extraporder2         [i] = 0;
                            maxquasiset_chichi0_freepts_extraporder2      [i] = 0;
                            maxquasiset_chixi0_freepts_extraporder2       [i] = 0;
                            maxquasiset_xixi0_freepts_extraporder2        [i] = 0;
                            maxquasiset_trace0_freepts_extraporder2       [i] = 0;
                            maxquasiset_massdensity0_freepts_extraporder2 [i] = 0;
                            maxbdyphi0_freepts_extraporder2               [i] = 0;  
                            minquasiset_tt0_freepts_extraporder2          [i] = 0;
                            minquasiset_tchi0_freepts_extraporder2        [i] = 0;
                            minquasiset_txi0_freepts_extraporder2         [i] = 0;
                            minquasiset_chichi0_freepts_extraporder2      [i] = 0;
                            minquasiset_chixi0_freepts_extraporder2       [i] = 0;
                            minquasiset_xixi0_freepts_extraporder2        [i] = 0;
                            minquasiset_trace0_freepts_extraporder2       [i] = 0;
                            minquasiset_massdensity0_freepts_extraporder2 [i] = 0;
                            minbdyphi0_freepts_extraporder2               [i] = 0;  
                            quasiset_tt0_freepts_extraporder2             [i] = 0;
                            quasiset_tchi0_freepts_extraporder2           [i] = 0;
                            quasiset_txi0_freepts_extraporder2            [i] = 0;
                            quasiset_chichi0_freepts_extraporder2         [i] = 0;
                            quasiset_chixi0_freepts_extraporder2          [i] = 0;
                            quasiset_xixi0_freepts_extraporder2           [i] = 0;
                            quasiset_trace0_freepts_extraporder2          [i] = 0;
                            quasiset_massdensity0_freepts_extraporder2    [i] = 0;
                            bdyphi0_freepts_extraporder2                  [i] = 0;  
                            xextrap0_freepts_extraporder2                 [i] = 0;
                            yextrap0_freepts_extraporder2                 [i] = 0;
                            zextrap0_freepts_extraporder2                 [i] = 0;  
                        }
                        *AdS_mass0_freepts_extraporder2                    = 0; 
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_freepts_extraporder2=0;
                        if (my_rank==0)
                        {
                            ie_bdy_freepts_extraporder2=vecbdypoints_freepts_extraporder2[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_freepts_extraporder2=is_bdy_freepts_extraporder2+vecbdypoints_freepts_extraporder2[j];
                            }
                            ie_bdy_freepts_extraporder2=is_bdy_freepts_extraporder2+vecbdypoints_freepts_extraporder2[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_freepts_extraporder2[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_freepts_extraporder2[i]=dsplsbdypoints_freepts_extraporder2[i]+vecbdypoints_freepts_extraporder2[j];
                                }
                            }
                        }   
                        xyzextrap_(xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);    
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,xextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,yextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,zextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);    //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FREE POINTS, SECOND ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_freepts_extraporder2 = malloc(sizeof(real));
                                chiextrap0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                                xiextrap0_freepts_extraporder2  = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real)); 
                                chixiextrap_(rhoextrap0_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2);  
                                basebdy_Nchi_freepts_extraporder2=0;//initialize
                                basebdy_Nxi_freepts_extraporder2=0; //initialize
                                bdyn_(&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2);  
                                rhobdy0_freepts_extraporder2 = malloc(sizeof(real));
                                chibdy0_freepts_extraporder2 = malloc(basebdy_Nchi_freepts_extraporder2*sizeof(real));
                                xibdy0_freepts_extraporder2  = malloc(basebdy_Nxi_freepts_extraporder2*sizeof(real));   
                            }
                        }   
                        extrap_bdyphi_freepts_(bdyphi_freepts_extraporder2,
                                        leadordcoeff_phi1,
                                        xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,
                                        chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,
                                        &bdy_extrap_order,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_freepts_(quasiset_tt_freepts_extraporder2,quasiset_tchi_freepts_extraporder2,quasiset_txi_freepts_extraporder2,
                                quasiset_chichi_freepts_extraporder2,quasiset_chixi_freepts_extraporder2,
                                quasiset_xixi_freepts_extraporder2,
                                quasiset_trace_freepts_extraporder2,
                                quasiset_massdensity_freepts_extraporder2,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,
                                chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,
                                &bdy_extrap_order,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_freepts_extraporder2; i<ie_bdy_freepts_extraporder2; i++)
                        {
                            lquasiset_tt0_freepts_extraporder2           [i] = quasiset_tt_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_tchi0_freepts_extraporder2         [i] = quasiset_tchi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_txi0_freepts_extraporder2          [i] = quasiset_txi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_chichi0_freepts_extraporder2       [i] = quasiset_chichi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_chixi0_freepts_extraporder2        [i] = quasiset_chixi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_xixi0_freepts_extraporder2         [i] = quasiset_xixi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_trace0_freepts_extraporder2        [i] = quasiset_trace_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_massdensity0_freepts_extraporder2  [i] = quasiset_massdensity_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lbdyphi0_freepts_extraporder2                [i] = bdyphi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_extrap_freepts

                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_fixedpts_extraporder1 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_fixedpts_extraporder1 = malloc(uniSize*sizeof(int)); 
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_fixedpts_extraporder1,1,MPI_INT,vecbdypoints_fixedpts_extraporder1,1,MPI_INT,MPI_COMM_WORLD);   
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);    
                        quasiset_tt_fixedpts_extraporder1              = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_tchi_fixedpts_extraporder1            = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_txi_fixedpts_extraporder1             = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chichi_fixedpts_extraporder1          = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chixi_fixedpts_extraporder1           = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_xixi_fixedpts_extraporder1            = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_trace_fixedpts_extraporder1           = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_massdensity_fixedpts_extraporder1     = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        bdyphi_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        xextrap_fixedpts_extraporder1                  = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        yextrap_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        zextrap_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));    
                        lquasiset_tt0_fixedpts_extraporder1            = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_tchi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_txi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_chichi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_chixi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_xixi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_trace0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_massdensity0_fixedpts_extraporder1   = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lbdyphi0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        maxquasiset_tt0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_tchi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_txi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_chichi0_fixedpts_extraporder1      = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_chixi0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_xixi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_trace0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_massdensity0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxbdyphi0_fixedpts_extraporder1               = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        minquasiset_tt0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_tchi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_txi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_chichi0_fixedpts_extraporder1      = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_chixi0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_xixi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_trace0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_massdensity0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minbdyphi0_fixedpts_extraporder1               = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        quasiset_tt0_fixedpts_extraporder1             = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_tchi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_txi0_fixedpts_extraporder1            = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chichi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chixi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_xixi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_trace0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_massdensity0_fixedpts_extraporder1    = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        AdS_mass0_fixedpts_extraporder1                = malloc(sizeof(real));
                        bdyphi0_fixedpts_extraporder1                  = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        xextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        yextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        zextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        //initialize
                        for (i=0;i<numbdypoints_fixedpts_extraporder1;i++)
                        {
                            quasiset_tt_fixedpts_extraporder1             [i] = 0;
                            quasiset_tchi_fixedpts_extraporder1           [i] = 0;
                            quasiset_txi_fixedpts_extraporder1            [i] = 0;
                            quasiset_chichi_fixedpts_extraporder1         [i] = 0;
                            quasiset_chixi_fixedpts_extraporder1          [i] = 0;
                            quasiset_xixi_fixedpts_extraporder1           [i] = 0;
                            quasiset_trace_fixedpts_extraporder1          [i] = 0;
                            quasiset_massdensity_fixedpts_extraporder1    [i] = 0;
                            bdyphi_fixedpts_extraporder1                  [i] = 0;  
                            xextrap_fixedpts_extraporder1                 [i] = 0;
                            yextrap_fixedpts_extraporder1                 [i] = 0;
                            zextrap_fixedpts_extraporder1                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_fixedpts_extraporder1;i++)
                        {   
                            lquasiset_tt0_fixedpts_extraporder1            [i] = 0;
                            lquasiset_tchi0_fixedpts_extraporder1          [i] = 0;
                            lquasiset_txi0_fixedpts_extraporder1           [i] = 0;
                            lquasiset_chichi0_fixedpts_extraporder1        [i] = 0;
                            lquasiset_chixi0_fixedpts_extraporder1         [i] = 0;
                            lquasiset_xixi0_fixedpts_extraporder1          [i] = 0;
                            lquasiset_trace0_fixedpts_extraporder1         [i] = 0;
                            lquasiset_massdensity0_fixedpts_extraporder1   [i] = 0;
                            lbdyphi0_fixedpts_extraporder1                 [i] = 0; 
                            maxquasiset_tt0_fixedpts_extraporder1          [i] = 0;
                            maxquasiset_tchi0_fixedpts_extraporder1        [i] = 0;
                            maxquasiset_txi0_fixedpts_extraporder1         [i] = 0;
                            maxquasiset_chichi0_fixedpts_extraporder1      [i] = 0;
                            maxquasiset_chixi0_fixedpts_extraporder1       [i] = 0;
                            maxquasiset_xixi0_fixedpts_extraporder1        [i] = 0;
                            maxquasiset_trace0_fixedpts_extraporder1       [i] = 0;
                            maxquasiset_massdensity0_fixedpts_extraporder1 [i] = 0;
                            maxbdyphi0_fixedpts_extraporder1               [i] = 0; 
                            minquasiset_tt0_fixedpts_extraporder1          [i] = 0;
                            minquasiset_tchi0_fixedpts_extraporder1        [i] = 0;
                            minquasiset_txi0_fixedpts_extraporder1         [i] = 0;
                            minquasiset_chichi0_fixedpts_extraporder1      [i] = 0;
                            minquasiset_chixi0_fixedpts_extraporder1       [i] = 0;
                            minquasiset_xixi0_fixedpts_extraporder1        [i] = 0;
                            minquasiset_trace0_fixedpts_extraporder1       [i] = 0;
                            minquasiset_massdensity0_fixedpts_extraporder1 [i] = 0;
                            minbdyphi0_fixedpts_extraporder1               [i] = 0; 
                            quasiset_tt0_fixedpts_extraporder1             [i] = 0;
                            quasiset_tchi0_fixedpts_extraporder1           [i] = 0;
                            quasiset_txi0_fixedpts_extraporder1            [i] = 0;
                            quasiset_chichi0_fixedpts_extraporder1         [i] = 0;
                            quasiset_chixi0_fixedpts_extraporder1          [i] = 0;
                            quasiset_xixi0_fixedpts_extraporder1           [i] = 0;
                            quasiset_trace0_fixedpts_extraporder1          [i] = 0;
                            quasiset_massdensity0_fixedpts_extraporder1    [i] = 0;
                            bdyphi0_fixedpts_extraporder1                  [i] = 0; 
                            xextrap0_fixedpts_extraporder1                 [i] = 0;
                            yextrap0_fixedpts_extraporder1                 [i] = 0;
                            zextrap0_fixedpts_extraporder1                 [i] = 0; 
                        }
                        *AdS_mass0_fixedpts_extraporder1                    = 0;    
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_fixedpts_extraporder1=0;
                        if (my_rank==0)
                        {
                            ie_bdy_fixedpts_extraporder1=vecbdypoints_fixedpts_extraporder1[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_fixedpts_extraporder1=is_bdy_fixedpts_extraporder1+vecbdypoints_fixedpts_extraporder1[j];
                            }
                            ie_bdy_fixedpts_extraporder1=is_bdy_fixedpts_extraporder1+vecbdypoints_fixedpts_extraporder1[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_fixedpts_extraporder1[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_fixedpts_extraporder1[i]=dsplsbdypoints_fixedpts_extraporder1[i]+vecbdypoints_fixedpts_extraporder1[j];
                                }
                            }
                        }

                        xyzextrap_(xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);   
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,xextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,yextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,zextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);   //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FIXED POINTS, FIRST ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_fixedpts_extraporder1 = malloc(sizeof(real));
                                chiextrap0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                                xiextrap0_fixedpts_extraporder1  = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));   
                                chixiextrap_(rhoextrap0_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1);   
                                basebdy_Nchi_fixedpts_extraporder1=0;//initialize
                                basebdy_Nxi_fixedpts_extraporder1=0; //initialize
                                bdyn_(&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1); 
                                rhobdy0_fixedpts_extraporder1 = malloc(sizeof(real));
                                chibdy0_fixedpts_extraporder1 = malloc(basebdy_Nchi_fixedpts_extraporder1*sizeof(real));
                                xibdy0_fixedpts_extraporder1  = malloc(basebdy_Nxi_fixedpts_extraporder1*sizeof(real)); 
                            }
                        }   
                        extrap_bdyphi_fixedpts_(bdyphi_fixedpts_extraporder1,
                                        leadordcoeff_phi1,
                                        xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,
                                        chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,
                                        &bdy_extrap_order,
                                        &ind_distance_fixedpts,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_fixedpts_(quasiset_tt_fixedpts_extraporder1,quasiset_tchi_fixedpts_extraporder1,quasiset_txi_fixedpts_extraporder1,
                                quasiset_chichi_fixedpts_extraporder1,quasiset_chixi_fixedpts_extraporder1,
                                quasiset_xixi_fixedpts_extraporder1,
                                quasiset_trace_fixedpts_extraporder1,
                                quasiset_massdensity_fixedpts_extraporder1,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,
                                chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,
                                &bdy_extrap_order,
                                &ind_distance_fixedpts,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_fixedpts_extraporder1; i<ie_bdy_fixedpts_extraporder1; i++)
                        {
                            lquasiset_tt0_fixedpts_extraporder1           [i] = quasiset_tt_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_tchi0_fixedpts_extraporder1         [i] = quasiset_tchi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_txi0_fixedpts_extraporder1          [i] = quasiset_txi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_chichi0_fixedpts_extraporder1       [i] = quasiset_chichi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_chixi0_fixedpts_extraporder1        [i] = quasiset_chixi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_xixi0_fixedpts_extraporder1         [i] = quasiset_xixi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_trace0_fixedpts_extraporder1        [i] = quasiset_trace_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_massdensity0_fixedpts_extraporder1  [i] = quasiset_massdensity_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lbdyphi0_fixedpts_extraporder1                [i] = bdyphi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];   //           *lAdS_mass0_fixedpts_extraporder1                 = *AdS_mass_fixedpts_extraporder1;
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_fixedpts_extraporder2 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_fixedpts_extraporder2 = malloc(uniSize*sizeof(int)); 
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_fixedpts_extraporder2,1,MPI_INT,vecbdypoints_fixedpts_extraporder2,1,MPI_INT,MPI_COMM_WORLD);   
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);    
                        quasiset_tt_fixedpts_extraporder2              = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_tchi_fixedpts_extraporder2            = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_txi_fixedpts_extraporder2             = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chichi_fixedpts_extraporder2          = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chixi_fixedpts_extraporder2           = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_xixi_fixedpts_extraporder2            = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_trace_fixedpts_extraporder2           = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_massdensity_fixedpts_extraporder2     = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        bdyphi_fixedpts_extraporder2                   = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        xextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        yextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        zextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        lquasiset_tt0_fixedpts_extraporder2            = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_tchi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_txi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_chichi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_chixi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_xixi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_trace0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_massdensity0_fixedpts_extraporder2   = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lbdyphi0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        maxquasiset_tt0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_tchi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_txi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_chichi0_fixedpts_extraporder2      = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_chixi0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_xixi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_trace0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_massdensity0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxbdyphi0_fixedpts_extraporder2               = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        minquasiset_tt0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_tchi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_txi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_chichi0_fixedpts_extraporder2      = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_chixi0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_xixi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_trace0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_massdensity0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minbdyphi0_fixedpts_extraporder2               = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        quasiset_tt0_fixedpts_extraporder2             = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_tchi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_txi0_fixedpts_extraporder2            = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chichi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chixi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_xixi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_trace0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_massdensity0_fixedpts_extraporder2    = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        AdS_mass0_fixedpts_extraporder2                = malloc(sizeof(real));
                        bdyphi0_fixedpts_extraporder2                  = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        xextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        yextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        zextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        //initialize
                        for (i=0;i<numbdypoints_fixedpts_extraporder2;i++)
                        {
                            quasiset_tt_fixedpts_extraporder2             [i] = 0;
                            quasiset_tchi_fixedpts_extraporder2           [i] = 0;
                            quasiset_txi_fixedpts_extraporder2            [i] = 0;
                            quasiset_chichi_fixedpts_extraporder2         [i] = 0;
                            quasiset_chixi_fixedpts_extraporder2          [i] = 0;
                            quasiset_xixi_fixedpts_extraporder2           [i] = 0;
                            quasiset_trace_fixedpts_extraporder2          [i] = 0;
                            quasiset_massdensity_fixedpts_extraporder2    [i] = 0;
                            bdyphi_fixedpts_extraporder2                  [i] = 0;  
                            xextrap_fixedpts_extraporder2                 [i] = 0;
                            yextrap_fixedpts_extraporder2                 [i] = 0;
                            zextrap_fixedpts_extraporder2                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_fixedpts_extraporder2;i++)
                        {   
                            lquasiset_tt0_fixedpts_extraporder2            [i] = 0;
                            lquasiset_tchi0_fixedpts_extraporder2          [i] = 0;
                            lquasiset_txi0_fixedpts_extraporder2           [i] = 0;
                            lquasiset_chichi0_fixedpts_extraporder2        [i] = 0;
                            lquasiset_chixi0_fixedpts_extraporder2         [i] = 0;
                            lquasiset_xixi0_fixedpts_extraporder2          [i] = 0;
                            lquasiset_trace0_fixedpts_extraporder2         [i] = 0;
                            lquasiset_massdensity0_fixedpts_extraporder2   [i] = 0;
                            lbdyphi0_fixedpts_extraporder2                 [i] = 0; 
                            maxquasiset_tt0_fixedpts_extraporder2          [i] = 0;
                            maxquasiset_tchi0_fixedpts_extraporder2        [i] = 0;
                            maxquasiset_txi0_fixedpts_extraporder2         [i] = 0;
                            maxquasiset_chichi0_fixedpts_extraporder2      [i] = 0;
                            maxquasiset_chixi0_fixedpts_extraporder2       [i] = 0;
                            maxquasiset_xixi0_fixedpts_extraporder2        [i] = 0;
                            maxquasiset_trace0_fixedpts_extraporder2       [i] = 0;
                            maxquasiset_massdensity0_fixedpts_extraporder2 [i] = 0;
                            maxbdyphi0_fixedpts_extraporder2               [i] = 0; 
                            minquasiset_tt0_fixedpts_extraporder2          [i] = 0;
                            minquasiset_tchi0_fixedpts_extraporder2        [i] = 0;
                            minquasiset_txi0_fixedpts_extraporder2         [i] = 0;
                            minquasiset_chichi0_fixedpts_extraporder2      [i] = 0;
                            minquasiset_chixi0_fixedpts_extraporder2       [i] = 0;
                            minquasiset_xixi0_fixedpts_extraporder2        [i] = 0;
                            minquasiset_trace0_fixedpts_extraporder2       [i] = 0;
                            minquasiset_massdensity0_fixedpts_extraporder2 [i] = 0;
                            minbdyphi0_fixedpts_extraporder2               [i] = 0; 
                            quasiset_tt0_fixedpts_extraporder2             [i] = 0;
                            quasiset_tchi0_fixedpts_extraporder2           [i] = 0;
                            quasiset_txi0_fixedpts_extraporder2            [i] = 0;
                            quasiset_chichi0_fixedpts_extraporder2         [i] = 0;
                            quasiset_chixi0_fixedpts_extraporder2          [i] = 0;
                            quasiset_xixi0_fixedpts_extraporder2           [i] = 0;
                            quasiset_trace0_fixedpts_extraporder2          [i] = 0;
                            quasiset_massdensity0_fixedpts_extraporder2    [i] = 0;
                            bdyphi0_fixedpts_extraporder2                  [i] = 0; 
                            xextrap0_fixedpts_extraporder2                 [i] = 0;
                            yextrap0_fixedpts_extraporder2                 [i] = 0;
                            zextrap0_fixedpts_extraporder2                 [i] = 0; 
                        }
                        *AdS_mass0_fixedpts_extraporder2                    = 0;    
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_fixedpts_extraporder2=0;
                        if (my_rank==0)
                        {
                            ie_bdy_fixedpts_extraporder2=vecbdypoints_fixedpts_extraporder2[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_fixedpts_extraporder2=is_bdy_fixedpts_extraporder2+vecbdypoints_fixedpts_extraporder2[j];
                            }
                            ie_bdy_fixedpts_extraporder2=is_bdy_fixedpts_extraporder2+vecbdypoints_fixedpts_extraporder2[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_fixedpts_extraporder2[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_fixedpts_extraporder2[i]=dsplsbdypoints_fixedpts_extraporder2[i]+vecbdypoints_fixedpts_extraporder2[j];
                                }
                            }
                        }   

                        xyzextrap_(xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);   
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,xextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,yextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,zextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);   //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FIXED POINTS, SECOND ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_fixedpts_extraporder2 = malloc(sizeof(real));
                                chiextrap0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                                xiextrap0_fixedpts_extraporder2  = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));   
                                chixiextrap_(rhoextrap0_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2);   
                                basebdy_Nchi_fixedpts_extraporder2=0;//initialize
                                basebdy_Nxi_fixedpts_extraporder2=0; //initialize
                                bdyn_(&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2); 
                                rhobdy0_fixedpts_extraporder2 = malloc(sizeof(real));
                                chibdy0_fixedpts_extraporder2 = malloc(basebdy_Nchi_fixedpts_extraporder2*sizeof(real));
                                xibdy0_fixedpts_extraporder2  = malloc(basebdy_Nxi_fixedpts_extraporder2*sizeof(real)); 
                            }
                        }   
                        extrap_bdyphi_fixedpts_(bdyphi_fixedpts_extraporder2,
                                        leadordcoeff_phi1,
                                        xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,
                                        chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,
                                        &bdy_extrap_order,
                                        &ind_distance_fixedpts,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_fixedpts_(quasiset_tt_fixedpts_extraporder2,quasiset_tchi_fixedpts_extraporder2,quasiset_txi_fixedpts_extraporder2,
                                quasiset_chichi_fixedpts_extraporder2,quasiset_chixi_fixedpts_extraporder2,
                                quasiset_xixi_fixedpts_extraporder2,
                                quasiset_trace_fixedpts_extraporder2,
                                quasiset_massdensity_fixedpts_extraporder2,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,
                                chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,
                                &bdy_extrap_order,
                                &ind_distance_fixedpts,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_fixedpts_extraporder2; i<ie_bdy_fixedpts_extraporder2; i++)
                        {
                            lquasiset_tt0_fixedpts_extraporder2           [i] = quasiset_tt_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_tchi0_fixedpts_extraporder2         [i] = quasiset_tchi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_txi0_fixedpts_extraporder2          [i] = quasiset_txi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_chichi0_fixedpts_extraporder2       [i] = quasiset_chichi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_chixi0_fixedpts_extraporder2        [i] = quasiset_chixi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_xixi0_fixedpts_extraporder2         [i] = quasiset_xixi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_trace0_fixedpts_extraporder2        [i] = quasiset_trace_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_massdensity0_fixedpts_extraporder2  [i] = quasiset_massdensity_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lbdyphi0_fixedpts_extraporder2                [i] = bdyphi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_fixedpts_extrap  
            } //closes condition on output_bdyquantities    
            valid=PAMR_next_g();
        }



        if (lsteps==0)
        {   
            if (output_kretsch && output_relkretschcentregrid)
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
                    FILE *fp;
                    sprintf(name,"%st_relkretschcentregrid.txt",AMRD_save_tag);
                    fp = fopen(name, "a+");
                    fprintf(fp,"%24.16e %24.16e \n",ct,*relkretschcentregrid0);
                    fclose(fp);
                }
            } 

            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder1,          maxquasiset_tt0_freepts_extraporder1,          basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder1,        maxquasiset_tchi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder1,         maxquasiset_txi0_freepts_extraporder1,         basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder1,      maxquasiset_chichi0_freepts_extraporder1,      basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder1,       maxquasiset_chixi0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder1,        maxquasiset_xixi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder1,       maxquasiset_trace0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder1, maxquasiset_massdensity0_freepts_extraporder1, basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder1,               maxbdyphi0_freepts_extraporder1,               basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);    
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder1,          minquasiset_tt0_freepts_extraporder1,          basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder1,        minquasiset_tchi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder1,         minquasiset_txi0_freepts_extraporder1,         basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder1,      minquasiset_chichi0_freepts_extraporder1,      basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder1,       minquasiset_chixi0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder1,        minquasiset_xixi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder1,       minquasiset_trace0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder1, minquasiset_massdensity0_freepts_extraporder1, basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder1,               minbdyphi0_freepts_extraporder1,               basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);    
                        for (i=0; i<basenumbdypoints_freepts_extraporder1; i++)
                        { 
                            if (uniSize>1)
                            {
                                    quasiset_tt0_freepts_extraporder1          [i] = maxquasiset_tt0_freepts_extraporder1          [i] + minquasiset_tt0_freepts_extraporder1          [i];
                                    quasiset_tchi0_freepts_extraporder1        [i] = maxquasiset_tchi0_freepts_extraporder1        [i] + minquasiset_tchi0_freepts_extraporder1        [i];
                                    quasiset_txi0_freepts_extraporder1         [i] = maxquasiset_txi0_freepts_extraporder1         [i] + minquasiset_txi0_freepts_extraporder1         [i];
                                    quasiset_chichi0_freepts_extraporder1      [i] = maxquasiset_chichi0_freepts_extraporder1      [i] + minquasiset_chichi0_freepts_extraporder1      [i];
                                    quasiset_chixi0_freepts_extraporder1       [i] = maxquasiset_chixi0_freepts_extraporder1       [i] + minquasiset_chixi0_freepts_extraporder1       [i];
                                    quasiset_xixi0_freepts_extraporder1        [i] = maxquasiset_xixi0_freepts_extraporder1        [i] + minquasiset_xixi0_freepts_extraporder1        [i];
                                    quasiset_trace0_freepts_extraporder1       [i] = maxquasiset_trace0_freepts_extraporder1       [i] + minquasiset_trace0_freepts_extraporder1       [i];
                                    quasiset_massdensity0_freepts_extraporder1 [i] = maxquasiset_massdensity0_freepts_extraporder1 [i] + minquasiset_massdensity0_freepts_extraporder1 [i];
                                    bdyphi0_freepts_extraporder1               [i] = maxbdyphi0_freepts_extraporder1               [i] + minbdyphi0_freepts_extraporder1               [i];
                                }
                                else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                                {
                                    quasiset_tt0_freepts_extraporder1          [i] = maxquasiset_tt0_freepts_extraporder1          [i];
                                    quasiset_tchi0_freepts_extraporder1        [i] = maxquasiset_tchi0_freepts_extraporder1        [i];
                                    quasiset_txi0_freepts_extraporder1         [i] = maxquasiset_txi0_freepts_extraporder1         [i];
                                    quasiset_chichi0_freepts_extraporder1      [i] = maxquasiset_chichi0_freepts_extraporder1      [i];
                                    quasiset_chixi0_freepts_extraporder1       [i] = maxquasiset_chixi0_freepts_extraporder1       [i];
                                    quasiset_xixi0_freepts_extraporder1        [i] = maxquasiset_xixi0_freepts_extraporder1        [i];
                                    quasiset_trace0_freepts_extraporder1       [i] = maxquasiset_trace0_freepts_extraporder1       [i];
                                    quasiset_massdensity0_freepts_extraporder1 [i] = maxquasiset_massdensity0_freepts_extraporder1 [i];
                                    bdyphi0_freepts_extraporder1               [i] = maxbdyphi0_freepts_extraporder1               [i];
                                }  
                        }       

                        if (my_rank==0)
                        {     
                            FILE *fp;  
    
                            if (alltimes_ascii)
                            {     
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
                                    
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }        
                                    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }      
        
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }
    
                            } //closes if(alltimes_ascii) condition 
    
                            if (timestep_ascii)
                            {     
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder1[j],yextrap0_freepts_extraporder1[j],zextrap0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder1[j],yextrap0_freepts_extraporder1[j],zextrap0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                            } //closes if (timestep_ascii) condition          
                        } //closes if (my_rank==0) condition            //the following bit computes and prints AdS_mass0 (see below) if we're running on only 1 process
                    
                        if (output_AdS_mass)
                        { 
                            if (uniSize>1)
                            {
                                if (my_rank==0) printf("\nTHE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT RELIABLE...NOT COMPUTING AdS MASS\n");
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FREE POINTS, FIRST ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");      
                                *rhobdy0_freepts_extraporder1=1;
                                chibdy_xibdy_(chibdy0_freepts_extraporder1,xibdy0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1,&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1);    
                                doubleintegralonsphere_(AdS_mass0_freepts_extraporder1,quasiset_massdensity0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,rhobdy0_freepts_extraporder1,chibdy0_freepts_extraporder1,xibdy0_freepts_extraporder1,&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1);          
                                FILE *fp;
                                sprintf(name,"AdSbdy_freepts_extraporder1_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_freepts_extraporder1);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder2,          maxquasiset_tt0_freepts_extraporder2,          basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder2,        maxquasiset_tchi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder2,         maxquasiset_txi0_freepts_extraporder2,         basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder2,      maxquasiset_chichi0_freepts_extraporder2,      basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder2,       maxquasiset_chixi0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder2,        maxquasiset_xixi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder2,       maxquasiset_trace0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder2, maxquasiset_massdensity0_freepts_extraporder2, basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder2,               maxbdyphi0_freepts_extraporder2,               basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder2,          minquasiset_tt0_freepts_extraporder2,          basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder2,        minquasiset_tchi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder2,         minquasiset_txi0_freepts_extraporder2,         basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder2,      minquasiset_chichi0_freepts_extraporder2,      basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder2,       minquasiset_chixi0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder2,        minquasiset_xixi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder2,       minquasiset_trace0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder2, minquasiset_massdensity0_freepts_extraporder2, basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder2,               minbdyphi0_freepts_extraporder2,               basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
                        for (i=0; i<basenumbdypoints_freepts_extraporder2; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_freepts_extraporder2          [i] = maxquasiset_tt0_freepts_extraporder2          [i] + minquasiset_tt0_freepts_extraporder2          [i];
                                quasiset_tchi0_freepts_extraporder2        [i] = maxquasiset_tchi0_freepts_extraporder2        [i] + minquasiset_tchi0_freepts_extraporder2        [i];
                                quasiset_txi0_freepts_extraporder2         [i] = maxquasiset_txi0_freepts_extraporder2         [i] + minquasiset_txi0_freepts_extraporder2         [i];
                                quasiset_chichi0_freepts_extraporder2      [i] = maxquasiset_chichi0_freepts_extraporder2      [i] + minquasiset_chichi0_freepts_extraporder2      [i];
                                quasiset_chixi0_freepts_extraporder2       [i] = maxquasiset_chixi0_freepts_extraporder2       [i] + minquasiset_chixi0_freepts_extraporder2       [i];
                                quasiset_xixi0_freepts_extraporder2        [i] = maxquasiset_xixi0_freepts_extraporder2        [i] + minquasiset_xixi0_freepts_extraporder2        [i];
                                quasiset_trace0_freepts_extraporder2       [i] = maxquasiset_trace0_freepts_extraporder2       [i] + minquasiset_trace0_freepts_extraporder2       [i];
                                quasiset_massdensity0_freepts_extraporder2 [i] = maxquasiset_massdensity0_freepts_extraporder2 [i] + minquasiset_massdensity0_freepts_extraporder2 [i];
                                bdyphi0_freepts_extraporder2               [i] = maxbdyphi0_freepts_extraporder2               [i] + minbdyphi0_freepts_extraporder2               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_freepts_extraporder2          [i] = maxquasiset_tt0_freepts_extraporder2          [i];
                                quasiset_tchi0_freepts_extraporder2        [i] = maxquasiset_tchi0_freepts_extraporder2        [i];
                                quasiset_txi0_freepts_extraporder2         [i] = maxquasiset_txi0_freepts_extraporder2         [i];
                                quasiset_chichi0_freepts_extraporder2      [i] = maxquasiset_chichi0_freepts_extraporder2      [i];
                                quasiset_chixi0_freepts_extraporder2       [i] = maxquasiset_chixi0_freepts_extraporder2       [i];
                                quasiset_xixi0_freepts_extraporder2        [i] = maxquasiset_xixi0_freepts_extraporder2        [i];
                                quasiset_trace0_freepts_extraporder2       [i] = maxquasiset_trace0_freepts_extraporder2       [i];
                                quasiset_massdensity0_freepts_extraporder2 [i] = maxquasiset_massdensity0_freepts_extraporder2 [i];
                                bdyphi0_freepts_extraporder2               [i] = maxbdyphi0_freepts_extraporder2               [i];
                            }
                        }   

                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                // save quasiset_ll as ascii
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder2[j],yextrap0_freepts_extraporder2[j],zextrap0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder2[j],yextrap0_freepts_extraporder2[j],zextrap0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                // save quasiset_ll as ascii
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FREE POINTS, SECOND ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                *rhobdy0_freepts_extraporder2=1;
                                chibdy_xibdy_(chibdy0_freepts_extraporder2,xibdy0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2,&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2);
                                doubleintegralonsphere_(AdS_mass0_freepts_extraporder2,quasiset_massdensity0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,rhobdy0_freepts_extraporder2,chibdy0_freepts_extraporder2,xibdy0_freepts_extraporder2,&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2); 
                                FILE *fp;
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_freepts_extraporder2);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_freepts_extrap  

                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder1,          maxquasiset_tt0_fixedpts_extraporder1,          basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder1,        maxquasiset_tchi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder1,         maxquasiset_txi0_fixedpts_extraporder1,         basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder1,      maxquasiset_chichi0_fixedpts_extraporder1,      basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder1,       maxquasiset_chixi0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder1,        maxquasiset_xixi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder1,       maxquasiset_trace0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder1, maxquasiset_massdensity0_fixedpts_extraporder1, basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder1,               maxbdyphi0_fixedpts_extraporder1,               basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder1,          minquasiset_tt0_fixedpts_extraporder1,          basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder1,        minquasiset_tchi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder1,         minquasiset_txi0_fixedpts_extraporder1,         basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder1,      minquasiset_chichi0_fixedpts_extraporder1,      basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder1,       minquasiset_chixi0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder1,        minquasiset_xixi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder1,       minquasiset_trace0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder1, minquasiset_massdensity0_fixedpts_extraporder1, basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder1,               minbdyphi0_fixedpts_extraporder1,               basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
                        for (i=0; i<basenumbdypoints_fixedpts_extraporder1; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_fixedpts_extraporder1          [i] = maxquasiset_tt0_fixedpts_extraporder1          [i] + minquasiset_tt0_fixedpts_extraporder1          [i];
                                quasiset_tchi0_fixedpts_extraporder1        [i] = maxquasiset_tchi0_fixedpts_extraporder1        [i] + minquasiset_tchi0_fixedpts_extraporder1        [i];
                                quasiset_txi0_fixedpts_extraporder1         [i] = maxquasiset_txi0_fixedpts_extraporder1         [i] + minquasiset_txi0_fixedpts_extraporder1         [i];
                                quasiset_chichi0_fixedpts_extraporder1      [i] = maxquasiset_chichi0_fixedpts_extraporder1      [i] + minquasiset_chichi0_fixedpts_extraporder1      [i];
                                quasiset_chixi0_fixedpts_extraporder1       [i] = maxquasiset_chixi0_fixedpts_extraporder1       [i] + minquasiset_chixi0_fixedpts_extraporder1       [i];
                                quasiset_xixi0_fixedpts_extraporder1        [i] = maxquasiset_xixi0_fixedpts_extraporder1        [i] + minquasiset_xixi0_fixedpts_extraporder1        [i];
                                quasiset_trace0_fixedpts_extraporder1       [i] = maxquasiset_trace0_fixedpts_extraporder1       [i] + minquasiset_trace0_fixedpts_extraporder1       [i];
                                quasiset_massdensity0_fixedpts_extraporder1 [i] = maxquasiset_massdensity0_fixedpts_extraporder1 [i] + minquasiset_massdensity0_fixedpts_extraporder1 [i];
                                bdyphi0_fixedpts_extraporder1               [i] = maxbdyphi0_fixedpts_extraporder1               [i] + minbdyphi0_fixedpts_extraporder1               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_fixedpts_extraporder1          [i] = maxquasiset_tt0_fixedpts_extraporder1          [i];
                                quasiset_tchi0_fixedpts_extraporder1        [i] = maxquasiset_tchi0_fixedpts_extraporder1        [i];
                                quasiset_txi0_fixedpts_extraporder1         [i] = maxquasiset_txi0_fixedpts_extraporder1         [i];
                                quasiset_chichi0_fixedpts_extraporder1      [i] = maxquasiset_chichi0_fixedpts_extraporder1      [i];
                                quasiset_chixi0_fixedpts_extraporder1       [i] = maxquasiset_chixi0_fixedpts_extraporder1       [i];
                                quasiset_xixi0_fixedpts_extraporder1        [i] = maxquasiset_xixi0_fixedpts_extraporder1        [i];
                                quasiset_trace0_fixedpts_extraporder1       [i] = maxquasiset_trace0_fixedpts_extraporder1       [i];
                                quasiset_massdensity0_fixedpts_extraporder1 [i] = maxquasiset_massdensity0_fixedpts_extraporder1 [i];
                                bdyphi0_fixedpts_extraporder1               [i] = maxbdyphi0_fixedpts_extraporder1               [i];
                            }
                        }       

                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder1[j],yextrap0_fixedpts_extraporder1[j],zextrap0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder1[j],yextrap0_fixedpts_extraporder1[j],zextrap0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FIXED POINTS, FIRST ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                *rhobdy0_fixedpts_extraporder1=1;
                                chibdy_xibdy_(chibdy0_fixedpts_extraporder1,xibdy0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1,&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1);
                                doubleintegralonsphere_(AdS_mass0_fixedpts_extraporder1,quasiset_massdensity0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,rhobdy0_fixedpts_extraporder1,chibdy0_fixedpts_extraporder1,xibdy0_fixedpts_extraporder1,&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1);  
                                FILE *fp;
                                sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_fixedpts_extraporder1);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder2,          maxquasiset_tt0_fixedpts_extraporder2,          basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder2,        maxquasiset_tchi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder2,         maxquasiset_txi0_fixedpts_extraporder2,         basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder2,      maxquasiset_chichi0_fixedpts_extraporder2,      basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder2,       maxquasiset_chixi0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder2,        maxquasiset_xixi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder2,       maxquasiset_trace0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder2, maxquasiset_massdensity0_fixedpts_extraporder2, basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder2,               maxbdyphi0_fixedpts_extraporder2,               basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder2,          minquasiset_tt0_fixedpts_extraporder2,          basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder2,        minquasiset_tchi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder2,         minquasiset_txi0_fixedpts_extraporder2,         basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder2,      minquasiset_chichi0_fixedpts_extraporder2,      basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder2,       minquasiset_chixi0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder2,        minquasiset_xixi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder2,       minquasiset_trace0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder2, minquasiset_massdensity0_fixedpts_extraporder2, basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder2,               minbdyphi0_fixedpts_extraporder2,               basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
                        for (i=0; i<basenumbdypoints_fixedpts_extraporder2; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_fixedpts_extraporder2          [i] = maxquasiset_tt0_fixedpts_extraporder2          [i] + minquasiset_tt0_fixedpts_extraporder2          [i];
                                quasiset_tchi0_fixedpts_extraporder2        [i] = maxquasiset_tchi0_fixedpts_extraporder2        [i] + minquasiset_tchi0_fixedpts_extraporder2        [i];
                                quasiset_txi0_fixedpts_extraporder2         [i] = maxquasiset_txi0_fixedpts_extraporder2         [i] + minquasiset_txi0_fixedpts_extraporder2         [i];
                                quasiset_chichi0_fixedpts_extraporder2      [i] = maxquasiset_chichi0_fixedpts_extraporder2      [i] + minquasiset_chichi0_fixedpts_extraporder2      [i];
                                quasiset_chixi0_fixedpts_extraporder2       [i] = maxquasiset_chixi0_fixedpts_extraporder2       [i] + minquasiset_chixi0_fixedpts_extraporder2       [i];
                                quasiset_xixi0_fixedpts_extraporder2        [i] = maxquasiset_xixi0_fixedpts_extraporder2        [i] + minquasiset_xixi0_fixedpts_extraporder2        [i];
                                quasiset_trace0_fixedpts_extraporder2       [i] = maxquasiset_trace0_fixedpts_extraporder2       [i] + minquasiset_trace0_fixedpts_extraporder2       [i];
                                quasiset_massdensity0_fixedpts_extraporder2 [i] = maxquasiset_massdensity0_fixedpts_extraporder2 [i] + minquasiset_massdensity0_fixedpts_extraporder2 [i];
                                bdyphi0_fixedpts_extraporder2               [i] = maxbdyphi0_fixedpts_extraporder2               [i] + minbdyphi0_fixedpts_extraporder2               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_fixedpts_extraporder2          [i] = maxquasiset_tt0_fixedpts_extraporder2          [i];
                                quasiset_tchi0_fixedpts_extraporder2        [i] = maxquasiset_tchi0_fixedpts_extraporder2        [i];
                                quasiset_txi0_fixedpts_extraporder2         [i] = maxquasiset_txi0_fixedpts_extraporder2         [i];
                                quasiset_chichi0_fixedpts_extraporder2      [i] = maxquasiset_chichi0_fixedpts_extraporder2      [i];
                                quasiset_chixi0_fixedpts_extraporder2       [i] = maxquasiset_chixi0_fixedpts_extraporder2       [i];
                                quasiset_xixi0_fixedpts_extraporder2        [i] = maxquasiset_xixi0_fixedpts_extraporder2        [i];
                                quasiset_trace0_fixedpts_extraporder2       [i] = maxquasiset_trace0_fixedpts_extraporder2       [i];
                                quasiset_massdensity0_fixedpts_extraporder2 [i] = maxquasiset_massdensity0_fixedpts_extraporder2 [i];
                                bdyphi0_fixedpts_extraporder2               [i] = maxbdyphi0_fixedpts_extraporder2               [i];
                            }
                        }   
                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder2[j],yextrap0_fixedpts_extraporder2[j],zextrap0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder2[j],yextrap0_fixedpts_extraporder2[j],zextrap0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FIXED POINTS, SECOND ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                *rhobdy0_fixedpts_extraporder2=1;
                                chibdy_xibdy_(chibdy0_fixedpts_extraporder2,xibdy0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2,&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2);
                                doubleintegralonsphere_(AdS_mass0_fixedpts_extraporder2,quasiset_massdensity0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,rhobdy0_fixedpts_extraporder2,chibdy0_fixedpts_extraporder2,xibdy0_fixedpts_extraporder2,&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2);  
                                FILE *fp;
                                sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_fixedpts_extraporder2);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_fixedpts_extrap  
            }//closes output_bdyquantities if-condition
        }//closes lsteps==0 if-condition
    }// closes L==Lc && ct==0 if-condition  


    //free memory for boundary quantities
    if (L==Lc && ct==0)
    {   
        valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
        while(valid)
        {   
            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {      
                        free(vecbdypoints_freepts_extraporder1);
                        free(dsplsbdypoints_freepts_extraporder1);
                        free(quasiset_tt_freepts_extraporder1);
                        free(quasiset_tchi_freepts_extraporder1);
                        free(quasiset_txi_freepts_extraporder1);
                        free(quasiset_chichi_freepts_extraporder1);
                        free(quasiset_chixi_freepts_extraporder1);
                        free(quasiset_xixi_freepts_extraporder1);
                        free(quasiset_trace_freepts_extraporder1);
                        free(quasiset_massdensity_freepts_extraporder1);
                        free(bdyphi_freepts_extraporder1);  
                        free(xextrap_freepts_extraporder1);
                        free(yextrap_freepts_extraporder1);
                        free(zextrap_freepts_extraporder1);   
                        free(lquasiset_tt0_freepts_extraporder1);        
                        free(lquasiset_tchi0_freepts_extraporder1);
                        free(lquasiset_txi0_freepts_extraporder1);
                        free(lquasiset_chichi0_freepts_extraporder1);
                        free(lquasiset_chixi0_freepts_extraporder1);
                        free(lquasiset_xixi0_freepts_extraporder1);
                        free(lquasiset_trace0_freepts_extraporder1);
                        free(lquasiset_massdensity0_freepts_extraporder1);
                        free(lbdyphi0_freepts_extraporder1);      
                        free(maxquasiset_tt0_freepts_extraporder1);        
                        free(maxquasiset_tchi0_freepts_extraporder1);
                        free(maxquasiset_txi0_freepts_extraporder1);
                        free(maxquasiset_chichi0_freepts_extraporder1);
                        free(maxquasiset_chixi0_freepts_extraporder1);
                        free(maxquasiset_xixi0_freepts_extraporder1);
                        free(maxquasiset_trace0_freepts_extraporder1);
                        free(maxquasiset_massdensity0_freepts_extraporder1);
                        free(maxbdyphi0_freepts_extraporder1);    
                        free(minquasiset_tt0_freepts_extraporder1);
                        free(minquasiset_tchi0_freepts_extraporder1);
                        free(minquasiset_txi0_freepts_extraporder1);
                        free(minquasiset_chichi0_freepts_extraporder1);
                        free(minquasiset_chixi0_freepts_extraporder1);
                        free(minquasiset_xixi0_freepts_extraporder1);
                        free(minquasiset_trace0_freepts_extraporder1);
                        free(minquasiset_massdensity0_freepts_extraporder1);
                        free(minbdyphi0_freepts_extraporder1);    
                        free(quasiset_tt0_freepts_extraporder1);
                        free(quasiset_tchi0_freepts_extraporder1);
                        free(quasiset_txi0_freepts_extraporder1);
                        free(quasiset_chichi0_freepts_extraporder1);
                        free(quasiset_chixi0_freepts_extraporder1);
                        free(quasiset_xixi0_freepts_extraporder1);
                        free(quasiset_trace0_freepts_extraporder1);
                        free(quasiset_massdensity0_freepts_extraporder1);
                        free(bdyphi0_freepts_extraporder1);
                        free(AdS_mass0_freepts_extraporder1);     
                        free(xextrap0_freepts_extraporder1);
                        free(yextrap0_freepts_extraporder1);
                        free(zextrap0_freepts_extraporder1);          
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_freepts_extraporder1);
                            free(chiextrap0_freepts_extraporder1);
                            free(xiextrap0_freepts_extraporder1);             
                            free(rhobdy0_freepts_extraporder1);
                            free(chibdy0_freepts_extraporder1);
                            free(xibdy0_freepts_extraporder1);    
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {   
                        free(vecbdypoints_freepts_extraporder2);
                        free(dsplsbdypoints_freepts_extraporder2);
                        free(quasiset_tt_freepts_extraporder2);
                        free(quasiset_tchi_freepts_extraporder2);
                        free(quasiset_txi_freepts_extraporder2);
                        free(quasiset_chichi_freepts_extraporder2);
                        free(quasiset_chixi_freepts_extraporder2);
                        free(quasiset_xixi_freepts_extraporder2);
                        free(quasiset_trace_freepts_extraporder2);
                        free(quasiset_massdensity_freepts_extraporder2);
                        free(bdyphi_freepts_extraporder2);  
                        free(xextrap_freepts_extraporder2);
                        free(yextrap_freepts_extraporder2);
                        free(zextrap_freepts_extraporder2); 
                        free(lquasiset_tt0_freepts_extraporder2);
                        free(lquasiset_tchi0_freepts_extraporder2);
                        free(lquasiset_txi0_freepts_extraporder2);
                        free(lquasiset_chichi0_freepts_extraporder2);
                        free(lquasiset_chixi0_freepts_extraporder2);
                        free(lquasiset_xixi0_freepts_extraporder2);
                        free(lquasiset_trace0_freepts_extraporder2);
                        free(lquasiset_massdensity0_freepts_extraporder2);
                        free(lbdyphi0_freepts_extraporder2);    
                        free(maxquasiset_tt0_freepts_extraporder2);
                        free(maxquasiset_tchi0_freepts_extraporder2);
                        free(maxquasiset_txi0_freepts_extraporder2);
                        free(maxquasiset_chichi0_freepts_extraporder2);
                        free(maxquasiset_chixi0_freepts_extraporder2);
                        free(maxquasiset_xixi0_freepts_extraporder2);
                        free(maxquasiset_trace0_freepts_extraporder2);
                        free(maxquasiset_massdensity0_freepts_extraporder2);
                        free(maxbdyphi0_freepts_extraporder2);  
                        free(minquasiset_tt0_freepts_extraporder2);
                        free(minquasiset_tchi0_freepts_extraporder2);
                        free(minquasiset_txi0_freepts_extraporder2);
                        free(minquasiset_chichi0_freepts_extraporder2);
                        free(minquasiset_chixi0_freepts_extraporder2);
                        free(minquasiset_xixi0_freepts_extraporder2);
                        free(minquasiset_trace0_freepts_extraporder2);
                        free(minquasiset_massdensity0_freepts_extraporder2);
                        free(minbdyphi0_freepts_extraporder2);  
                        free(quasiset_tt0_freepts_extraporder2);
                        free(quasiset_tchi0_freepts_extraporder2);
                        free(quasiset_txi0_freepts_extraporder2);
                        free(quasiset_chichi0_freepts_extraporder2);
                        free(quasiset_chixi0_freepts_extraporder2);
                        free(quasiset_xixi0_freepts_extraporder2);
                        free(quasiset_trace0_freepts_extraporder2);
                        free(quasiset_massdensity0_freepts_extraporder2);
                        free(bdyphi0_freepts_extraporder2);
                        free(AdS_mass0_freepts_extraporder2);   
                        free(xextrap0_freepts_extraporder2);
                        free(yextrap0_freepts_extraporder2);
                        free(zextrap0_freepts_extraporder2);    
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_freepts_extraporder2);
                            free(chiextrap0_freepts_extraporder2);
                            free(xiextrap0_freepts_extraporder2);   
                            free(rhobdy0_freepts_extraporder2);
                            free(chibdy0_freepts_extraporder2);
                            free(xibdy0_freepts_extraporder2);  
                        }   
                    } //closes condition on output_bdy_extraporder2 
                }//closes condition on bdy_freepts_extrap       
                
                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {   
                        free(vecbdypoints_fixedpts_extraporder1);
                        free(dsplsbdypoints_fixedpts_extraporder1);
                        free(quasiset_tt_fixedpts_extraporder1);
                        free(quasiset_tchi_fixedpts_extraporder1);
                        free(quasiset_txi_fixedpts_extraporder1);
                        free(quasiset_chichi_fixedpts_extraporder1);
                        free(quasiset_chixi_fixedpts_extraporder1);
                        free(quasiset_xixi_fixedpts_extraporder1);
                        free(quasiset_trace_fixedpts_extraporder1);
                        free(quasiset_massdensity_fixedpts_extraporder1);
                        free(bdyphi_fixedpts_extraporder1); 
                        free(xextrap_fixedpts_extraporder1);
                        free(yextrap_fixedpts_extraporder1);
                        free(zextrap_fixedpts_extraporder1);    
                        free(lquasiset_tt0_fixedpts_extraporder1);
                        free(lquasiset_tchi0_fixedpts_extraporder1);
                        free(lquasiset_txi0_fixedpts_extraporder1);
                        free(lquasiset_chichi0_fixedpts_extraporder1);
                        free(lquasiset_chixi0_fixedpts_extraporder1);
                        free(lquasiset_xixi0_fixedpts_extraporder1);
                        free(lquasiset_trace0_fixedpts_extraporder1);
                        free(lquasiset_massdensity0_fixedpts_extraporder1);
                        free(lbdyphi0_fixedpts_extraporder1);   
                        free(maxquasiset_tt0_fixedpts_extraporder1);
                        free(maxquasiset_tchi0_fixedpts_extraporder1);
                        free(maxquasiset_txi0_fixedpts_extraporder1);
                        free(maxquasiset_chichi0_fixedpts_extraporder1);
                        free(maxquasiset_chixi0_fixedpts_extraporder1);
                        free(maxquasiset_xixi0_fixedpts_extraporder1);
                        free(maxquasiset_trace0_fixedpts_extraporder1);
                        free(maxquasiset_massdensity0_fixedpts_extraporder1);
                        free(maxbdyphi0_fixedpts_extraporder1); 
                        free(minquasiset_tt0_fixedpts_extraporder1);
                        free(minquasiset_tchi0_fixedpts_extraporder1);
                        free(minquasiset_txi0_fixedpts_extraporder1);
                        free(minquasiset_chichi0_fixedpts_extraporder1);
                        free(minquasiset_chixi0_fixedpts_extraporder1);
                        free(minquasiset_xixi0_fixedpts_extraporder1);
                        free(minquasiset_trace0_fixedpts_extraporder1);
                        free(minquasiset_massdensity0_fixedpts_extraporder1);
                        free(minbdyphi0_fixedpts_extraporder1); 
                        free(quasiset_tt0_fixedpts_extraporder1);
                        free(quasiset_tchi0_fixedpts_extraporder1);
                        free(quasiset_txi0_fixedpts_extraporder1);
                        free(quasiset_chichi0_fixedpts_extraporder1);
                        free(quasiset_chixi0_fixedpts_extraporder1);
                        free(quasiset_xixi0_fixedpts_extraporder1);
                        free(quasiset_trace0_fixedpts_extraporder1);
                        free(quasiset_massdensity0_fixedpts_extraporder1);
                        free(bdyphi0_fixedpts_extraporder1);
                        free(AdS_mass0_fixedpts_extraporder1);  
                        free(xextrap0_fixedpts_extraporder1);
                        free(yextrap0_fixedpts_extraporder1);
                        free(zextrap0_fixedpts_extraporder1);   
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_fixedpts_extraporder1);
                            free(chiextrap0_fixedpts_extraporder1);
                            free(xiextrap0_fixedpts_extraporder1);  
                            free(rhobdy0_fixedpts_extraporder1);
                            free(chibdy0_fixedpts_extraporder1);
                            free(xibdy0_fixedpts_extraporder1); 
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {   
                        free(vecbdypoints_fixedpts_extraporder2);
                        free(dsplsbdypoints_fixedpts_extraporder2);
                        free(quasiset_tt_fixedpts_extraporder2);
                        free(quasiset_tchi_fixedpts_extraporder2);
                        free(quasiset_txi_fixedpts_extraporder2);
                        free(quasiset_chichi_fixedpts_extraporder2);
                        free(quasiset_chixi_fixedpts_extraporder2);
                        free(quasiset_xixi_fixedpts_extraporder2);
                        free(quasiset_trace_fixedpts_extraporder2);
                        free(quasiset_massdensity_fixedpts_extraporder2);
                        free(bdyphi_fixedpts_extraporder2); 
                        free(xextrap_fixedpts_extraporder2);
                        free(yextrap_fixedpts_extraporder2);
                        free(zextrap_fixedpts_extraporder2);    
                        free(lquasiset_tt0_fixedpts_extraporder2);
                        free(lquasiset_tchi0_fixedpts_extraporder2);
                        free(lquasiset_txi0_fixedpts_extraporder2);
                        free(lquasiset_chichi0_fixedpts_extraporder2);
                        free(lquasiset_chixi0_fixedpts_extraporder2);
                        free(lquasiset_xixi0_fixedpts_extraporder2);
                        free(lquasiset_trace0_fixedpts_extraporder2);
                        free(lquasiset_massdensity0_fixedpts_extraporder2);
                        free(lbdyphi0_fixedpts_extraporder2);   
                        free(maxquasiset_tt0_fixedpts_extraporder2);
                        free(maxquasiset_tchi0_fixedpts_extraporder2);
                        free(maxquasiset_txi0_fixedpts_extraporder2);
                        free(maxquasiset_chichi0_fixedpts_extraporder2);
                        free(maxquasiset_chixi0_fixedpts_extraporder2);
                        free(maxquasiset_xixi0_fixedpts_extraporder2);
                        free(maxquasiset_trace0_fixedpts_extraporder2);
                        free(maxquasiset_massdensity0_fixedpts_extraporder2);
                        free(maxbdyphi0_fixedpts_extraporder2); 
                        free(minquasiset_tt0_fixedpts_extraporder2);
                        free(minquasiset_tchi0_fixedpts_extraporder2);
                        free(minquasiset_txi0_fixedpts_extraporder2);
                        free(minquasiset_chichi0_fixedpts_extraporder2);
                        free(minquasiset_chixi0_fixedpts_extraporder2);
                        free(minquasiset_xixi0_fixedpts_extraporder2);
                        free(minquasiset_trace0_fixedpts_extraporder2);
                        free(minquasiset_massdensity0_fixedpts_extraporder2);
                        free(minbdyphi0_fixedpts_extraporder2); 
                        free(quasiset_tt0_fixedpts_extraporder2);
                        free(quasiset_tchi0_fixedpts_extraporder2);
                        free(quasiset_txi0_fixedpts_extraporder2);
                        free(quasiset_chichi0_fixedpts_extraporder2);
                        free(quasiset_chixi0_fixedpts_extraporder2);
                        free(quasiset_xixi0_fixedpts_extraporder2);
                        free(quasiset_trace0_fixedpts_extraporder2);
                        free(quasiset_massdensity0_fixedpts_extraporder2);
                        free(bdyphi0_fixedpts_extraporder2);
                        free(AdS_mass0_fixedpts_extraporder2);  
                        free(xextrap0_fixedpts_extraporder2);
                        free(yextrap0_fixedpts_extraporder2);
                        free(zextrap0_fixedpts_extraporder2);   
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_fixedpts_extraporder2);
                            free(chiextrap0_fixedpts_extraporder2);
                            free(xiextrap0_fixedpts_extraporder2);  
                            free(rhobdy0_fixedpts_extraporder2);
                            free(chibdy0_fixedpts_extraporder2);
                            free(xibdy0_fixedpts_extraporder2); 
                        }   
                    } //closes condition on output_bdy_extraporder2 
                }//closes condition on bdy_fixedpts_extrap  
            } //closes condition on output_bdyquantities
            valid=PAMR_next_g();
        }  
    } //closes condition on (L==Lc)&&(ct==0)    


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
            AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0,output_metricatAH); // AH finder 
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
                    if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1,output_metricatAH)))
                    { 
                        // shrink old initial-guess surface
                        if (my_rank==0 && AMRD_evo_trace>=1) printf("... still can't find one (min_resid=%lf)\n"
                            "... shrinking old initial-guess surface by a factor %lf\n",AH_min_resid1,1/AH_reset_scale[l]);
                        for (i=0; i<AH_Nchi[l]*AH_Nphi[l]; i++) AH_R[l][i]=prev_AH_R[i]/AH_reset_scale[l];  
                        if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2,output_metricatAH)))
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
                                if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid2,output_metricatAH)))
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
                                if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid1,output_metricatAH)))
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
                                if (!(AH[l]=find_apph(&M,&J,&c_equat,&c_polar,found_AH[l],&AH_min_resid0,output_metricatAH)))
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
                AH_bbox[3]=2*M_PI;  //               if (AH_xc[l][1]<dy) {AH_bbox[3]=M_PI;} else {AH_bbox[3]=2*M_PI;} //this line comes from codes where y in [0,1]
                int rank=2;
                sprintf(name,"%sAH_R_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_R[l]);
                sprintf(name,"%sAH_theta_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_theta[l]);    
                if (output_metricatAH)
                {   //               for (i=0;i<(AH_Nchi[l]*AH_Nphi[l]);i++) {printf("i=%i,AH_g0_tt[l][i]=%lf\n",i,AH_g0_tt[l][i]);}
                sprintf(name,"%sAH_g0_tt_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_tt[l]);
                sprintf(name,"%sAH_g0_tx_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_tx[l]);
                sprintf(name,"%sAH_g0_ty_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_ty[l]);
                sprintf(name,"%sAH_g0_tz_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_tz[l]);
                sprintf(name,"%sAH_g0_xx_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_xx[l]);
                sprintf(name,"%sAH_g0_xy_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_xy[l]);
                sprintf(name,"%sAH_g0_xz_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_xz[l]);
                sprintf(name,"%sAH_g0_yy_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_yy[l]);
                sprintf(name,"%sAH_g0_yz_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_yz[l]);
                sprintf(name,"%sAH_g0_psi_%i",AMRD_save_tag,l);
                gft_out_bbox(name,ct,AH_shape,rank,AH_bbox,AH_g0_psi[l]);
                }   
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
    {   //     remove_redundant_AH();
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
            rhoh=(-1 + sqrt(1 + pow(rh,2)))/rh; //         ex_r[0][0]=ex_r[0][1]=ex_r[0][2]=rhoh*(1-ex_rbuf[0]);
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

int post_tstep_global_first=1;

void AdS4D_post_tstep(int L)
{
    char name[256];
    int itrace=1,valid;
    static int local_first = 1; 
    real ct;
    int n,i,j,k,ind,j_red,Lf,Lc;
    int is,ie,js,je,ks,ke;  

    //printf("AdS4D_post_tstep is called");
    //fflush(stdout); 

    ct = PAMR_get_time(L);  
    Lf=PAMR_get_max_lev(PAMR_AMRH);
    Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; elise Lc=1;


    // qs objects at t>0, when at coarsest level L=Lc
    if (L==Lc)
    {
        int lsteps=AMRD_lsteps[Lc-1];
        int ivecNt=AMRD_steps/AMRD_save_ivec0[3]+1; //+1 to include t=0
        real lmass,mass;    
        valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
        while(valid)
        {   
            ldptr();    
            //we save here the values of the Kretschmann scalar at the centre of the grid at times greater than 0, computed at the end of g_evo_opt //The value of the Kretschmann scalar at the centre of the grid at t=0 is saved in pre_io_calc.  
            if (output_relkretschcentregrid)
            {
                *lrelkretschcentregrid0= *relkretschcentregrid;
            }   
            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize); 
                        vecbdypoints_freepts_extraporder1 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_freepts_extraporder1 = malloc(uniSize*sizeof(int));    
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_freepts_extraporder1,1,MPI_INT,vecbdypoints_freepts_extraporder1,1,MPI_INT,MPI_COMM_WORLD); 
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
                        //printf("numbdypoints_freepts_extraporder1=%i\n",numbdypoints_freepts_extraporder1);
                        //fflush(stdout); 
                        quasiset_tt_freepts_extraporder1              = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_tchi_freepts_extraporder1            = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_txi_freepts_extraporder1             = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chichi_freepts_extraporder1          = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chixi_freepts_extraporder1           = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_xixi_freepts_extraporder1            = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_trace_freepts_extraporder1           = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_massdensity_freepts_extraporder1     = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        bdyphi_freepts_extraporder1                   = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));   
                        xextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        yextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));
                        zextrap_freepts_extraporder1                  = malloc((numbdypoints_freepts_extraporder1)*sizeof(real));   
                        lquasiset_tt0_freepts_extraporder1            = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_tchi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_txi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_chichi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_chixi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_xixi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_trace0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lquasiset_massdensity0_freepts_extraporder1   = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        lbdyphi0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        maxquasiset_tt0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_tchi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_txi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_chichi0_freepts_extraporder1      = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_chixi0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_xixi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_trace0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxquasiset_massdensity0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        maxbdyphi0_freepts_extraporder1               = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        minquasiset_tt0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_tchi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_txi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_chichi0_freepts_extraporder1      = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_chixi0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_xixi0_freepts_extraporder1        = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_trace0_freepts_extraporder1       = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minquasiset_massdensity0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        minbdyphi0_freepts_extraporder1               = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        quasiset_tt0_freepts_extraporder1             = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_tchi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_txi0_freepts_extraporder1            = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chichi0_freepts_extraporder1         = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_chixi0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_xixi0_freepts_extraporder1           = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_trace0_freepts_extraporder1          = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        quasiset_massdensity0_freepts_extraporder1    = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        AdS_mass0_freepts_extraporder1                = malloc(sizeof(real));
                        bdyphi0_freepts_extraporder1                  = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));    
                        xextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        yextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                        zextrap0_freepts_extraporder1                 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));   
                        //initialize
                        for (i=0;i<numbdypoints_freepts_extraporder1;i++)
                        {
                            quasiset_tt_freepts_extraporder1             [i] = 0;
                            quasiset_tchi_freepts_extraporder1           [i] = 0;
                            quasiset_txi_freepts_extraporder1            [i] = 0;
                            quasiset_chichi_freepts_extraporder1         [i] = 0;
                            quasiset_chixi_freepts_extraporder1          [i] = 0;
                            quasiset_xixi_freepts_extraporder1           [i] = 0;
                            quasiset_trace_freepts_extraporder1          [i] = 0;
                            quasiset_massdensity_freepts_extraporder1    [i] = 0;
                            bdyphi_freepts_extraporder1                  [i] = 0;   
                            xextrap_freepts_extraporder1                 [i] = 0;
                            yextrap_freepts_extraporder1                 [i] = 0;
                            zextrap_freepts_extraporder1                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_freepts_extraporder1;i++)
                        {   
                            lquasiset_tt0_freepts_extraporder1            [i] = 0;
                            lquasiset_tchi0_freepts_extraporder1          [i] = 0;
                            lquasiset_txi0_freepts_extraporder1           [i] = 0;
                            lquasiset_chichi0_freepts_extraporder1        [i] = 0;
                            lquasiset_chixi0_freepts_extraporder1         [i] = 0;
                            lquasiset_xixi0_freepts_extraporder1          [i] = 0;
                            lquasiset_trace0_freepts_extraporder1         [i] = 0;
                            lquasiset_massdensity0_freepts_extraporder1   [i] = 0;
                            lbdyphi0_freepts_extraporder1                 [i] = 0;             
                            maxquasiset_tt0_freepts_extraporder1          [i] = 0;
                            maxquasiset_tchi0_freepts_extraporder1        [i] = 0;
                            maxquasiset_txi0_freepts_extraporder1         [i] = 0;
                            maxquasiset_chichi0_freepts_extraporder1      [i] = 0;
                            maxquasiset_chixi0_freepts_extraporder1       [i] = 0;
                            maxquasiset_xixi0_freepts_extraporder1        [i] = 0;
                            maxquasiset_trace0_freepts_extraporder1       [i] = 0;
                            maxquasiset_massdensity0_freepts_extraporder1 [i] = 0;
                            maxbdyphi0_freepts_extraporder1               [i] = 0;  
                            minquasiset_tt0_freepts_extraporder1          [i] = 0;
                            minquasiset_tchi0_freepts_extraporder1        [i] = 0;
                            minquasiset_txi0_freepts_extraporder1         [i] = 0;
                            minquasiset_chichi0_freepts_extraporder1      [i] = 0;
                            minquasiset_chixi0_freepts_extraporder1       [i] = 0;
                            minquasiset_xixi0_freepts_extraporder1        [i] = 0;
                            minquasiset_trace0_freepts_extraporder1       [i] = 0;
                            minquasiset_massdensity0_freepts_extraporder1 [i] = 0;
                            minbdyphi0_freepts_extraporder1               [i] = 0;  
                            quasiset_tt0_freepts_extraporder1             [i] = 0;
                            quasiset_tchi0_freepts_extraporder1           [i] = 0;
                            quasiset_txi0_freepts_extraporder1            [i] = 0;
                            quasiset_chichi0_freepts_extraporder1         [i] = 0;
                            quasiset_chixi0_freepts_extraporder1          [i] = 0;
                            quasiset_xixi0_freepts_extraporder1           [i] = 0;
                            quasiset_trace0_freepts_extraporder1          [i] = 0;
                            quasiset_massdensity0_freepts_extraporder1    [i] = 0;
                            bdyphi0_freepts_extraporder1                  [i] = 0;   
                            xextrap0_freepts_extraporder1                 [i] = 0;
                            yextrap0_freepts_extraporder1                 [i] = 0;
                            zextrap0_freepts_extraporder1                 [i] = 0;  
                        }
                        *AdS_mass0_freepts_extraporder1                    = 0; 
                        
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_freepts_extraporder1=0;
                        if (my_rank==0)
                        {
                            ie_bdy_freepts_extraporder1=vecbdypoints_freepts_extraporder1[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_freepts_extraporder1=is_bdy_freepts_extraporder1+vecbdypoints_freepts_extraporder1[j];
                            }
                            ie_bdy_freepts_extraporder1=is_bdy_freepts_extraporder1+vecbdypoints_freepts_extraporder1[my_rank];
                        }      
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_freepts_extraporder1[i]=0;
                        }         
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_freepts_extraporder1[i]=dsplsbdypoints_freepts_extraporder1[i]+vecbdypoints_freepts_extraporder1[j];
                                }
                            }
                        }   
                        
                        xyzextrap_(xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);    
                        
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,xextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,yextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_freepts_extraporder1,numbdypoints_freepts_extraporder1,MPI_DOUBLE,zextrap0_freepts_extraporder1,vecbdypoints_freepts_extraporder1,dsplsbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);    //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0...");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FREE POINTS, FIRST ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_freepts_extraporder1 = malloc(sizeof(real));
                                chiextrap0_freepts_extraporder1 = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real));
                                xiextrap0_freepts_extraporder1  = malloc((basenumbdypoints_freepts_extraporder1)*sizeof(real)); 
                                chixiextrap_(rhoextrap0_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1);  
                                basebdy_Nchi_freepts_extraporder1=0;//initialize
                                basebdy_Nxi_freepts_extraporder1=0; //initialize
                                bdyn_(&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1);  
                                rhobdy0_freepts_extraporder1 = malloc(sizeof(real));
                                chibdy0_freepts_extraporder1 = malloc(basebdy_Nchi_freepts_extraporder1*sizeof(real));
                                xibdy0_freepts_extraporder1  = malloc(basebdy_Nxi_freepts_extraporder1*sizeof(real));   
                            }
                        }   


                        extrap_bdyphi_freepts_(bdyphi_freepts_extraporder1,
                                        leadordcoeff_phi1,
                                        xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,
                                        chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,
                                        &bdy_extrap_order,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_freepts_(quasiset_tt_freepts_extraporder1,quasiset_tchi_freepts_extraporder1,quasiset_txi_freepts_extraporder1,
                                quasiset_chichi_freepts_extraporder1,quasiset_chixi_freepts_extraporder1,
                                quasiset_xixi_freepts_extraporder1,
                                quasiset_trace_freepts_extraporder1,
                                quasiset_massdensity_freepts_extraporder1,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_freepts_extraporder1,yextrap_freepts_extraporder1,zextrap_freepts_extraporder1,
                                chrbdy_freepts_extraporder1,&numbdypoints_freepts_extraporder1,
                                &bdy_extrap_order,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_freepts_extraporder1; i<ie_bdy_freepts_extraporder1; i++)
                        {
                            lquasiset_tt0_freepts_extraporder1           [i] = quasiset_tt_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_tchi0_freepts_extraporder1         [i] = quasiset_tchi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_txi0_freepts_extraporder1          [i] = quasiset_txi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_chichi0_freepts_extraporder1       [i] = quasiset_chichi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_chixi0_freepts_extraporder1        [i] = quasiset_chixi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_xixi0_freepts_extraporder1         [i] = quasiset_xixi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_trace0_freepts_extraporder1        [i] = quasiset_trace_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lquasiset_massdensity0_freepts_extraporder1  [i] = quasiset_massdensity_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                            lbdyphi0_freepts_extraporder1                [i] = bdyphi_freepts_extraporder1[i-is_bdy_freepts_extraporder1];
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_freepts_extraporder2 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_freepts_extraporder2 = malloc(uniSize*sizeof(int));  
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_freepts_extraporder2,1,MPI_INT,vecbdypoints_freepts_extraporder2,1,MPI_INT,MPI_COMM_WORLD); 
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  
                        quasiset_tt_freepts_extraporder2              = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_tchi_freepts_extraporder2            = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_txi_freepts_extraporder2             = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chichi_freepts_extraporder2          = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chixi_freepts_extraporder2           = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_xixi_freepts_extraporder2            = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_trace_freepts_extraporder2           = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_massdensity_freepts_extraporder2     = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        bdyphi_freepts_extraporder2                   = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));   
                        xextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        yextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));
                        zextrap_freepts_extraporder2                  = malloc((numbdypoints_freepts_extraporder2)*sizeof(real));   
                        lquasiset_tt0_freepts_extraporder2            = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_tchi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_txi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_chichi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_chixi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_xixi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_trace0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lquasiset_massdensity0_freepts_extraporder2   = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        lbdyphi0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        maxquasiset_tt0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_tchi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_txi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_chichi0_freepts_extraporder2      = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_chixi0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_xixi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_trace0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxquasiset_massdensity0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        maxbdyphi0_freepts_extraporder2               = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        minquasiset_tt0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_tchi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_txi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_chichi0_freepts_extraporder2      = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_chixi0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_xixi0_freepts_extraporder2        = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_trace0_freepts_extraporder2       = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minquasiset_massdensity0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        minbdyphi0_freepts_extraporder2               = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        quasiset_tt0_freepts_extraporder2             = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_tchi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_txi0_freepts_extraporder2            = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chichi0_freepts_extraporder2         = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_chixi0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_xixi0_freepts_extraporder2           = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_trace0_freepts_extraporder2          = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        quasiset_massdensity0_freepts_extraporder2    = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        AdS_mass0_freepts_extraporder2                = malloc(sizeof(real));
                        bdyphi0_freepts_extraporder2                  = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        xextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        yextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                        zextrap0_freepts_extraporder2                 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));   
                        //initialize
                        for (i=0;i<numbdypoints_freepts_extraporder2;i++)
                        {
                            quasiset_tt_freepts_extraporder2             [i] = 0;
                            quasiset_tchi_freepts_extraporder2           [i] = 0;
                            quasiset_txi_freepts_extraporder2            [i] = 0;
                            quasiset_chichi_freepts_extraporder2         [i] = 0;
                            quasiset_chixi_freepts_extraporder2          [i] = 0;
                            quasiset_xixi_freepts_extraporder2           [i] = 0;
                            quasiset_trace_freepts_extraporder2          [i] = 0;
                            quasiset_massdensity_freepts_extraporder2    [i] = 0;
                            bdyphi_freepts_extraporder2                  [i] = 0;   
                            xextrap_freepts_extraporder2                 [i] = 0;
                            yextrap_freepts_extraporder2                 [i] = 0;
                            zextrap_freepts_extraporder2                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_freepts_extraporder2;i++)
                        {   
                            lquasiset_tt0_freepts_extraporder2            [i] = 0;
                            lquasiset_tchi0_freepts_extraporder2          [i] = 0;
                            lquasiset_txi0_freepts_extraporder2           [i] = 0;
                            lquasiset_chichi0_freepts_extraporder2        [i] = 0;
                            lquasiset_chixi0_freepts_extraporder2         [i] = 0;
                            lquasiset_xixi0_freepts_extraporder2          [i] = 0;
                            lquasiset_trace0_freepts_extraporder2         [i] = 0;
                            lquasiset_massdensity0_freepts_extraporder2   [i] = 0;
                            lbdyphi0_freepts_extraporder2                 [i] = 0;  
                            maxquasiset_tt0_freepts_extraporder2          [i] = 0;
                            maxquasiset_tchi0_freepts_extraporder2        [i] = 0;
                            maxquasiset_txi0_freepts_extraporder2         [i] = 0;
                            maxquasiset_chichi0_freepts_extraporder2      [i] = 0;
                            maxquasiset_chixi0_freepts_extraporder2       [i] = 0;
                            maxquasiset_xixi0_freepts_extraporder2        [i] = 0;
                            maxquasiset_trace0_freepts_extraporder2       [i] = 0;
                            maxquasiset_massdensity0_freepts_extraporder2 [i] = 0;
                            maxbdyphi0_freepts_extraporder2               [i] = 0;  
                            minquasiset_tt0_freepts_extraporder2          [i] = 0;
                            minquasiset_tchi0_freepts_extraporder2        [i] = 0;
                            minquasiset_txi0_freepts_extraporder2         [i] = 0;
                            minquasiset_chichi0_freepts_extraporder2      [i] = 0;
                            minquasiset_chixi0_freepts_extraporder2       [i] = 0;
                            minquasiset_xixi0_freepts_extraporder2        [i] = 0;
                            minquasiset_trace0_freepts_extraporder2       [i] = 0;
                            minquasiset_massdensity0_freepts_extraporder2 [i] = 0;
                            minbdyphi0_freepts_extraporder2               [i] = 0;  
                            quasiset_tt0_freepts_extraporder2             [i] = 0;
                            quasiset_tchi0_freepts_extraporder2           [i] = 0;
                            quasiset_txi0_freepts_extraporder2            [i] = 0;
                            quasiset_chichi0_freepts_extraporder2         [i] = 0;
                            quasiset_chixi0_freepts_extraporder2          [i] = 0;
                            quasiset_xixi0_freepts_extraporder2           [i] = 0;
                            quasiset_trace0_freepts_extraporder2          [i] = 0;
                            quasiset_massdensity0_freepts_extraporder2    [i] = 0;
                            bdyphi0_freepts_extraporder2                  [i] = 0;  
                            xextrap0_freepts_extraporder2                 [i] = 0;
                            yextrap0_freepts_extraporder2                 [i] = 0;
                            zextrap0_freepts_extraporder2                 [i] = 0;  
                        }
                        *AdS_mass0_freepts_extraporder2                    = 0; 
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_freepts_extraporder2=0;
                        if (my_rank==0)
                        {
                            ie_bdy_freepts_extraporder2=vecbdypoints_freepts_extraporder2[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_freepts_extraporder2=is_bdy_freepts_extraporder2+vecbdypoints_freepts_extraporder2[j];
                            }
                            ie_bdy_freepts_extraporder2=is_bdy_freepts_extraporder2+vecbdypoints_freepts_extraporder2[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_freepts_extraporder2[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_freepts_extraporder2[i]=dsplsbdypoints_freepts_extraporder2[i]+vecbdypoints_freepts_extraporder2[j];
                                }
                            }
                        }   
                        xyzextrap_(xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);    
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,xextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,yextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_freepts_extraporder2,numbdypoints_freepts_extraporder2,MPI_DOUBLE,zextrap0_freepts_extraporder2,vecbdypoints_freepts_extraporder2,dsplsbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);    //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FREE POINTS, SECOND ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_freepts_extraporder2 = malloc(sizeof(real));
                                chiextrap0_freepts_extraporder2 = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real));
                                xiextrap0_freepts_extraporder2  = malloc((basenumbdypoints_freepts_extraporder2)*sizeof(real)); 
                                chixiextrap_(rhoextrap0_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2);  
                                basebdy_Nchi_freepts_extraporder2=0;//initialize
                                basebdy_Nxi_freepts_extraporder2=0; //initialize
                                bdyn_(&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2);  
                                rhobdy0_freepts_extraporder2 = malloc(sizeof(real));
                                chibdy0_freepts_extraporder2 = malloc(basebdy_Nchi_freepts_extraporder2*sizeof(real));
                                xibdy0_freepts_extraporder2  = malloc(basebdy_Nxi_freepts_extraporder2*sizeof(real));   
                            }
                        }   
                        extrap_bdyphi_freepts_(bdyphi_freepts_extraporder2,
                                        leadordcoeff_phi1,
                                        xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,
                                        chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,
                                        &bdy_extrap_order,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_freepts_(quasiset_tt_freepts_extraporder2,quasiset_tchi_freepts_extraporder2,quasiset_txi_freepts_extraporder2,
                                quasiset_chichi_freepts_extraporder2,quasiset_chixi_freepts_extraporder2,
                                quasiset_xixi_freepts_extraporder2,
                                quasiset_trace_freepts_extraporder2,
                                quasiset_massdensity_freepts_extraporder2,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_freepts_extraporder2,yextrap_freepts_extraporder2,zextrap_freepts_extraporder2,
                                chrbdy_freepts_extraporder2,&numbdypoints_freepts_extraporder2,
                                &bdy_extrap_order,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_freepts_extraporder2; i<ie_bdy_freepts_extraporder2; i++)
                        {
                            lquasiset_tt0_freepts_extraporder2           [i] = quasiset_tt_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_tchi0_freepts_extraporder2         [i] = quasiset_tchi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_txi0_freepts_extraporder2          [i] = quasiset_txi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_chichi0_freepts_extraporder2       [i] = quasiset_chichi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_chixi0_freepts_extraporder2        [i] = quasiset_chixi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_xixi0_freepts_extraporder2         [i] = quasiset_xixi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_trace0_freepts_extraporder2        [i] = quasiset_trace_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lquasiset_massdensity0_freepts_extraporder2  [i] = quasiset_massdensity_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                            lbdyphi0_freepts_extraporder2                [i] = bdyphi_freepts_extraporder2[i-is_bdy_freepts_extraporder2];
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_extrap_freepts

                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_fixedpts_extraporder1 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_fixedpts_extraporder1 = malloc(uniSize*sizeof(int)); 
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_fixedpts_extraporder1,1,MPI_INT,vecbdypoints_fixedpts_extraporder1,1,MPI_INT,MPI_COMM_WORLD);   
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);    
                        quasiset_tt_fixedpts_extraporder1              = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_tchi_fixedpts_extraporder1            = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_txi_fixedpts_extraporder1             = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chichi_fixedpts_extraporder1          = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chixi_fixedpts_extraporder1           = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_xixi_fixedpts_extraporder1            = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_trace_fixedpts_extraporder1           = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_massdensity_fixedpts_extraporder1     = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        bdyphi_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        xextrap_fixedpts_extraporder1                  = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        yextrap_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));
                        zextrap_fixedpts_extraporder1                   = malloc((numbdypoints_fixedpts_extraporder1)*sizeof(real));    
                        lquasiset_tt0_fixedpts_extraporder1            = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_tchi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_txi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_chichi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_chixi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_xixi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_trace0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lquasiset_massdensity0_fixedpts_extraporder1   = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        lbdyphi0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        maxquasiset_tt0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_tchi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_txi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_chichi0_fixedpts_extraporder1      = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_chixi0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_xixi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_trace0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxquasiset_massdensity0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        maxbdyphi0_fixedpts_extraporder1               = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        minquasiset_tt0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_tchi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_txi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_chichi0_fixedpts_extraporder1      = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_chixi0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_xixi0_fixedpts_extraporder1        = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_trace0_fixedpts_extraporder1       = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minquasiset_massdensity0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        minbdyphi0_fixedpts_extraporder1               = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        quasiset_tt0_fixedpts_extraporder1             = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_tchi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_txi0_fixedpts_extraporder1            = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chichi0_fixedpts_extraporder1         = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_chixi0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_xixi0_fixedpts_extraporder1           = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_trace0_fixedpts_extraporder1          = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        quasiset_massdensity0_fixedpts_extraporder1    = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        AdS_mass0_fixedpts_extraporder1                = malloc(sizeof(real));
                        bdyphi0_fixedpts_extraporder1                  = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        xextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        yextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                        zextrap0_fixedpts_extraporder1                 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real)); 
                        //initialize
                        for (i=0;i<numbdypoints_fixedpts_extraporder1;i++)
                        {
                            quasiset_tt_fixedpts_extraporder1             [i] = 0;
                            quasiset_tchi_fixedpts_extraporder1           [i] = 0;
                            quasiset_txi_fixedpts_extraporder1            [i] = 0;
                            quasiset_chichi_fixedpts_extraporder1         [i] = 0;
                            quasiset_chixi_fixedpts_extraporder1          [i] = 0;
                            quasiset_xixi_fixedpts_extraporder1           [i] = 0;
                            quasiset_trace_fixedpts_extraporder1          [i] = 0;
                            quasiset_massdensity_fixedpts_extraporder1    [i] = 0;
                            bdyphi_fixedpts_extraporder1                  [i] = 0;  
                            xextrap_fixedpts_extraporder1                 [i] = 0;
                            yextrap_fixedpts_extraporder1                 [i] = 0;
                            zextrap_fixedpts_extraporder1                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_fixedpts_extraporder1;i++)
                        {   
                            lquasiset_tt0_fixedpts_extraporder1            [i] = 0;
                            lquasiset_tchi0_fixedpts_extraporder1          [i] = 0;
                            lquasiset_txi0_fixedpts_extraporder1           [i] = 0;
                            lquasiset_chichi0_fixedpts_extraporder1        [i] = 0;
                            lquasiset_chixi0_fixedpts_extraporder1         [i] = 0;
                            lquasiset_xixi0_fixedpts_extraporder1          [i] = 0;
                            lquasiset_trace0_fixedpts_extraporder1         [i] = 0;
                            lquasiset_massdensity0_fixedpts_extraporder1   [i] = 0;
                            lbdyphi0_fixedpts_extraporder1                 [i] = 0; 
                            maxquasiset_tt0_fixedpts_extraporder1          [i] = 0;
                            maxquasiset_tchi0_fixedpts_extraporder1        [i] = 0;
                            maxquasiset_txi0_fixedpts_extraporder1         [i] = 0;
                            maxquasiset_chichi0_fixedpts_extraporder1      [i] = 0;
                            maxquasiset_chixi0_fixedpts_extraporder1       [i] = 0;
                            maxquasiset_xixi0_fixedpts_extraporder1        [i] = 0;
                            maxquasiset_trace0_fixedpts_extraporder1       [i] = 0;
                            maxquasiset_massdensity0_fixedpts_extraporder1 [i] = 0;
                            maxbdyphi0_fixedpts_extraporder1               [i] = 0; 
                            minquasiset_tt0_fixedpts_extraporder1          [i] = 0;
                            minquasiset_tchi0_fixedpts_extraporder1        [i] = 0;
                            minquasiset_txi0_fixedpts_extraporder1         [i] = 0;
                            minquasiset_chichi0_fixedpts_extraporder1      [i] = 0;
                            minquasiset_chixi0_fixedpts_extraporder1       [i] = 0;
                            minquasiset_xixi0_fixedpts_extraporder1        [i] = 0;
                            minquasiset_trace0_fixedpts_extraporder1       [i] = 0;
                            minquasiset_massdensity0_fixedpts_extraporder1 [i] = 0;
                            minbdyphi0_fixedpts_extraporder1               [i] = 0; 
                            quasiset_tt0_fixedpts_extraporder1             [i] = 0;
                            quasiset_tchi0_fixedpts_extraporder1           [i] = 0;
                            quasiset_txi0_fixedpts_extraporder1            [i] = 0;
                            quasiset_chichi0_fixedpts_extraporder1         [i] = 0;
                            quasiset_chixi0_fixedpts_extraporder1          [i] = 0;
                            quasiset_xixi0_fixedpts_extraporder1           [i] = 0;
                            quasiset_trace0_fixedpts_extraporder1          [i] = 0;
                            quasiset_massdensity0_fixedpts_extraporder1    [i] = 0;
                            bdyphi0_fixedpts_extraporder1                  [i] = 0; 
                            xextrap0_fixedpts_extraporder1                 [i] = 0;
                            yextrap0_fixedpts_extraporder1                 [i] = 0;
                            zextrap0_fixedpts_extraporder1                 [i] = 0; 
                        }
                        *AdS_mass0_fixedpts_extraporder1                    = 0;    
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_fixedpts_extraporder1=0;
                        if (my_rank==0)
                        {
                            ie_bdy_fixedpts_extraporder1=vecbdypoints_fixedpts_extraporder1[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_fixedpts_extraporder1=is_bdy_fixedpts_extraporder1+vecbdypoints_fixedpts_extraporder1[j];
                            }
                            ie_bdy_fixedpts_extraporder1=is_bdy_fixedpts_extraporder1+vecbdypoints_fixedpts_extraporder1[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_fixedpts_extraporder1[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_fixedpts_extraporder1[i]=dsplsbdypoints_fixedpts_extraporder1[i]+vecbdypoints_fixedpts_extraporder1[j];
                                }
                            }
                        }

                        xyzextrap_(xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);   
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,xextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,yextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_fixedpts_extraporder1,numbdypoints_fixedpts_extraporder1,MPI_DOUBLE,zextrap0_fixedpts_extraporder1,vecbdypoints_fixedpts_extraporder1,dsplsbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_COMM_WORLD);   //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FIXED POINTS, FIRST ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_fixedpts_extraporder1 = malloc(sizeof(real));
                                chiextrap0_fixedpts_extraporder1 = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));
                                xiextrap0_fixedpts_extraporder1  = malloc((basenumbdypoints_fixedpts_extraporder1)*sizeof(real));   
                                chixiextrap_(rhoextrap0_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1);   
                                basebdy_Nchi_fixedpts_extraporder1=0;//initialize
                                basebdy_Nxi_fixedpts_extraporder1=0; //initialize
                                bdyn_(&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1); 
                                rhobdy0_fixedpts_extraporder1 = malloc(sizeof(real));
                                chibdy0_fixedpts_extraporder1 = malloc(basebdy_Nchi_fixedpts_extraporder1*sizeof(real));
                                xibdy0_fixedpts_extraporder1  = malloc(basebdy_Nxi_fixedpts_extraporder1*sizeof(real)); 
                            }
                        }   
                        extrap_bdyphi_fixedpts_(bdyphi_fixedpts_extraporder1,
                                        leadordcoeff_phi1,
                                        xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,
                                        chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,
                                        &bdy_extrap_order,
                                        &ind_distance_fixedpts,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_fixedpts_(quasiset_tt_fixedpts_extraporder1,quasiset_tchi_fixedpts_extraporder1,quasiset_txi_fixedpts_extraporder1,
                                quasiset_chichi_fixedpts_extraporder1,quasiset_chixi_fixedpts_extraporder1,
                                quasiset_xixi_fixedpts_extraporder1,
                                quasiset_trace_fixedpts_extraporder1,
                                quasiset_massdensity_fixedpts_extraporder1,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_fixedpts_extraporder1,yextrap_fixedpts_extraporder1,zextrap_fixedpts_extraporder1,
                                chrbdy_fixedpts_extraporder1,&numbdypoints_fixedpts_extraporder1,
                                &bdy_extrap_order,
                                &ind_distance_fixedpts,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_fixedpts_extraporder1; i<ie_bdy_fixedpts_extraporder1; i++)
                        {
                            lquasiset_tt0_fixedpts_extraporder1           [i] = quasiset_tt_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_tchi0_fixedpts_extraporder1         [i] = quasiset_tchi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_txi0_fixedpts_extraporder1          [i] = quasiset_txi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_chichi0_fixedpts_extraporder1       [i] = quasiset_chichi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_chixi0_fixedpts_extraporder1        [i] = quasiset_chixi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_xixi0_fixedpts_extraporder1         [i] = quasiset_xixi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_trace0_fixedpts_extraporder1        [i] = quasiset_trace_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lquasiset_massdensity0_fixedpts_extraporder1  [i] = quasiset_massdensity_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];
                            lbdyphi0_fixedpts_extraporder1                [i] = bdyphi_fixedpts_extraporder1[i-is_bdy_fixedpts_extraporder1];   //           *lAdS_mass0_fixedpts_extraporder1                 = *AdS_mass_fixedpts_extraporder1;
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        MPI_Comm_size(MPI_COMM_WORLD,&uniSize);
                        vecbdypoints_fixedpts_extraporder2 = malloc(uniSize*sizeof(int));
                        dsplsbdypoints_fixedpts_extraporder2 = malloc(uniSize*sizeof(int)); 
                        //the ith element of vecbdypoints contains the number of nexttobdypoints identified by nexttobdypoints routine for the ith process
                        MPI_Allgather(&numbdypoints_fixedpts_extraporder2,1,MPI_INT,vecbdypoints_fixedpts_extraporder2,1,MPI_INT,MPI_COMM_WORLD);   
                        //basenumbdypoints contains the sum of the number of nexttobdypoints from all processes, i.e. the total number of nexttobdypoints, hence the total number of points at the boundary where we extrapolate the stress-energy tensor
                        MPI_Allreduce(&numbdypoints_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);    
                        quasiset_tt_fixedpts_extraporder2              = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_tchi_fixedpts_extraporder2            = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_txi_fixedpts_extraporder2             = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chichi_fixedpts_extraporder2          = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chixi_fixedpts_extraporder2           = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_xixi_fixedpts_extraporder2            = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_trace_fixedpts_extraporder2           = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_massdensity_fixedpts_extraporder2     = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        bdyphi_fixedpts_extraporder2                   = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        xextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        yextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real));
                        zextrap_fixedpts_extraporder2                  = malloc((numbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        lquasiset_tt0_fixedpts_extraporder2            = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_tchi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_txi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_chichi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_chixi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_xixi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_trace0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lquasiset_massdensity0_fixedpts_extraporder2   = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        lbdyphi0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        maxquasiset_tt0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_tchi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_txi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_chichi0_fixedpts_extraporder2      = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_chixi0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_xixi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_trace0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxquasiset_massdensity0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        maxbdyphi0_fixedpts_extraporder2               = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        minquasiset_tt0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_tchi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_txi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_chichi0_fixedpts_extraporder2      = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_chixi0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_xixi0_fixedpts_extraporder2        = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_trace0_fixedpts_extraporder2       = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minquasiset_massdensity0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        minbdyphi0_fixedpts_extraporder2               = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        quasiset_tt0_fixedpts_extraporder2             = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_tchi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_txi0_fixedpts_extraporder2            = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chichi0_fixedpts_extraporder2         = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_chixi0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_xixi0_fixedpts_extraporder2           = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_trace0_fixedpts_extraporder2          = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        quasiset_massdensity0_fixedpts_extraporder2    = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        AdS_mass0_fixedpts_extraporder2                = malloc(sizeof(real));
                        bdyphi0_fixedpts_extraporder2                  = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        xextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        yextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                        zextrap0_fixedpts_extraporder2                 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real)); 
                        //initialize
                        for (i=0;i<numbdypoints_fixedpts_extraporder2;i++)
                        {
                            quasiset_tt_fixedpts_extraporder2             [i] = 0;
                            quasiset_tchi_fixedpts_extraporder2           [i] = 0;
                            quasiset_txi_fixedpts_extraporder2            [i] = 0;
                            quasiset_chichi_fixedpts_extraporder2         [i] = 0;
                            quasiset_chixi_fixedpts_extraporder2          [i] = 0;
                            quasiset_xixi_fixedpts_extraporder2           [i] = 0;
                            quasiset_trace_fixedpts_extraporder2          [i] = 0;
                            quasiset_massdensity_fixedpts_extraporder2    [i] = 0;
                            bdyphi_fixedpts_extraporder2                  [i] = 0;  
                            xextrap_fixedpts_extraporder2                 [i] = 0;
                            yextrap_fixedpts_extraporder2                 [i] = 0;
                            zextrap_fixedpts_extraporder2                 [i] = 0;
                        }   
                        for (i=0;i<basenumbdypoints_fixedpts_extraporder2;i++)
                        {   
                            lquasiset_tt0_fixedpts_extraporder2            [i] = 0;
                            lquasiset_tchi0_fixedpts_extraporder2          [i] = 0;
                            lquasiset_txi0_fixedpts_extraporder2           [i] = 0;
                            lquasiset_chichi0_fixedpts_extraporder2        [i] = 0;
                            lquasiset_chixi0_fixedpts_extraporder2         [i] = 0;
                            lquasiset_xixi0_fixedpts_extraporder2          [i] = 0;
                            lquasiset_trace0_fixedpts_extraporder2         [i] = 0;
                            lquasiset_massdensity0_fixedpts_extraporder2   [i] = 0;
                            lbdyphi0_fixedpts_extraporder2                 [i] = 0; 
                            maxquasiset_tt0_fixedpts_extraporder2          [i] = 0;
                            maxquasiset_tchi0_fixedpts_extraporder2        [i] = 0;
                            maxquasiset_txi0_fixedpts_extraporder2         [i] = 0;
                            maxquasiset_chichi0_fixedpts_extraporder2      [i] = 0;
                            maxquasiset_chixi0_fixedpts_extraporder2       [i] = 0;
                            maxquasiset_xixi0_fixedpts_extraporder2        [i] = 0;
                            maxquasiset_trace0_fixedpts_extraporder2       [i] = 0;
                            maxquasiset_massdensity0_fixedpts_extraporder2 [i] = 0;
                            maxbdyphi0_fixedpts_extraporder2               [i] = 0; 
                            minquasiset_tt0_fixedpts_extraporder2          [i] = 0;
                            minquasiset_tchi0_fixedpts_extraporder2        [i] = 0;
                            minquasiset_txi0_fixedpts_extraporder2         [i] = 0;
                            minquasiset_chichi0_fixedpts_extraporder2      [i] = 0;
                            minquasiset_chixi0_fixedpts_extraporder2       [i] = 0;
                            minquasiset_xixi0_fixedpts_extraporder2        [i] = 0;
                            minquasiset_trace0_fixedpts_extraporder2       [i] = 0;
                            minquasiset_massdensity0_fixedpts_extraporder2 [i] = 0;
                            minbdyphi0_fixedpts_extraporder2               [i] = 0; 
                            quasiset_tt0_fixedpts_extraporder2             [i] = 0;
                            quasiset_tchi0_fixedpts_extraporder2           [i] = 0;
                            quasiset_txi0_fixedpts_extraporder2            [i] = 0;
                            quasiset_chichi0_fixedpts_extraporder2         [i] = 0;
                            quasiset_chixi0_fixedpts_extraporder2          [i] = 0;
                            quasiset_xixi0_fixedpts_extraporder2           [i] = 0;
                            quasiset_trace0_fixedpts_extraporder2          [i] = 0;
                            quasiset_massdensity0_fixedpts_extraporder2    [i] = 0;
                            bdyphi0_fixedpts_extraporder2                  [i] = 0; 
                            xextrap0_fixedpts_extraporder2                 [i] = 0;
                            yextrap0_fixedpts_extraporder2                 [i] = 0;
                            zextrap0_fixedpts_extraporder2                 [i] = 0; 
                        }
                        *AdS_mass0_fixedpts_extraporder2                    = 0;    
                        //we want the indices from is to ie to identify the bdypoints of each processor starting the count from the last bdypoint of the previous processor
                        is_bdy_fixedpts_extraporder2=0;
                        if (my_rank==0)
                        {
                            ie_bdy_fixedpts_extraporder2=vecbdypoints_fixedpts_extraporder2[0];
                        }
                        else
                        {
                            for (j=0; j<my_rank; j++)
                            {
                                is_bdy_fixedpts_extraporder2=is_bdy_fixedpts_extraporder2+vecbdypoints_fixedpts_extraporder2[j];
                            }
                            ie_bdy_fixedpts_extraporder2=is_bdy_fixedpts_extraporder2+vecbdypoints_fixedpts_extraporder2[my_rank];
                        }   
                        //the ith element of dsplsbdypoints contains the number of nexttobdypoints of the processor i-1. We need this array as displacement array for MPI_Allgatherv below.
                        for (i=0; i<uniSize; i++)
                        {
                            dsplsbdypoints_fixedpts_extraporder2[i]=0;
                        }   
                        for (i=0; i<uniSize; i++)
                        {
                            if (i!=0)
                            {
                                for (j=0; j<i; j++)
                                {
                                    dsplsbdypoints_fixedpts_extraporder2[i]=dsplsbdypoints_fixedpts_extraporder2[i]+vecbdypoints_fixedpts_extraporder2[j];
                                }
                            }
                        }   

                        xyzextrap_(xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,ghost_width);   
                        //x/y/zextrap0 are arrays with xextrap,yextrap,zextrap from all the processors one after the other
                        MPI_Allgatherv(xextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,xextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(yextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,yextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);
                        MPI_Allgatherv(zextrap_fixedpts_extraporder2,numbdypoints_fixedpts_extraporder2,MPI_DOUBLE,zextrap0_fixedpts_extraporder2,vecbdypoints_fixedpts_extraporder2,dsplsbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_COMM_WORLD);   //the following bit allocates memory to compute AdS_mass0 (see below) if we're running on only 1 process
                        if (output_AdS_mass)
                        {
                            if (uniSize>1)
                            {
                                if (my_rank==0) 
                                {
                                    printf("THE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT TRUSTWORTHY...\n NOT ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS\n");
                                    printf("...setting output_AdS_mass to 0");
                                }
                                output_AdS_mass=0;
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("RUNNING ON ONLY 1 PROCESS...ALLOCATING MEMORY FOR COMPUTATION OF AdS MASS FOR FIXED POINTS, SECOND ORDER EXTRAPOLATION ON ONLY 1 PROCESS\n");
                                rhoextrap0_fixedpts_extraporder2 = malloc(sizeof(real));
                                chiextrap0_fixedpts_extraporder2 = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));
                                xiextrap0_fixedpts_extraporder2  = malloc((basenumbdypoints_fixedpts_extraporder2)*sizeof(real));   
                                chixiextrap_(rhoextrap0_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2);   
                                basebdy_Nchi_fixedpts_extraporder2=0;//initialize
                                basebdy_Nxi_fixedpts_extraporder2=0; //initialize
                                bdyn_(&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2); 
                                rhobdy0_fixedpts_extraporder2 = malloc(sizeof(real));
                                chibdy0_fixedpts_extraporder2 = malloc(basebdy_Nchi_fixedpts_extraporder2*sizeof(real));
                                xibdy0_fixedpts_extraporder2  = malloc(basebdy_Nxi_fixedpts_extraporder2*sizeof(real)); 
                            }
                        }   
                        extrap_bdyphi_fixedpts_(bdyphi_fixedpts_extraporder2,
                                        leadordcoeff_phi1,
                                        xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,
                                        chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,
                                        &bdy_extrap_order,
                                        &ind_distance_fixedpts,
                                        x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        extrap_quasiset_fixedpts_(quasiset_tt_fixedpts_extraporder2,quasiset_tchi_fixedpts_extraporder2,quasiset_txi_fixedpts_extraporder2,
                                quasiset_chichi_fixedpts_extraporder2,quasiset_chixi_fixedpts_extraporder2,
                                quasiset_xixi_fixedpts_extraporder2,
                                quasiset_trace_fixedpts_extraporder2,
                                quasiset_massdensity_fixedpts_extraporder2,
                                quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
                                quasiset_chichi_ll,quasiset_chixi_ll,
                                quasiset_xixi_ll,
                                quasiset_tracell,
                                quasiset_massdensityll,
                                xextrap_fixedpts_extraporder2,yextrap_fixedpts_extraporder2,zextrap_fixedpts_extraporder2,
                                chrbdy_fixedpts_extraporder2,&numbdypoints_fixedpts_extraporder2,
                                &bdy_extrap_order,
                                &ind_distance_fixedpts,
                                x,y,z,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,&Nz,phys_bdy,ghost_width);    
                        //distributing the values of the quasiset components of each process over an array lquasiset_ll0 defined globally. This array will be different for each process, in fact it will be zero everywhere except for a certain position (next to the one for the previous processor) containing the values of quasiset_ll of a specific process. This is repeated after each step of the evolution. 
                        for (i=is_bdy_fixedpts_extraporder2; i<ie_bdy_fixedpts_extraporder2; i++)
                        {
                            lquasiset_tt0_fixedpts_extraporder2           [i] = quasiset_tt_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_tchi0_fixedpts_extraporder2         [i] = quasiset_tchi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_txi0_fixedpts_extraporder2          [i] = quasiset_txi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_chichi0_fixedpts_extraporder2       [i] = quasiset_chichi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_chixi0_fixedpts_extraporder2        [i] = quasiset_chixi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_xixi0_fixedpts_extraporder2         [i] = quasiset_xixi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_trace0_fixedpts_extraporder2        [i] = quasiset_trace_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lquasiset_massdensity0_fixedpts_extraporder2  [i] = quasiset_massdensity_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                            lbdyphi0_fixedpts_extraporder2                [i] = bdyphi_fixedpts_extraporder2[i-is_bdy_fixedpts_extraporder2];
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_fixedpts_extrap  
            } //closes condition on output_bdyquantities    
            valid=PAMR_next_g();
        }



        if ((lsteps%AMRD_save_ivec0[3]==0)&&(lsteps!=0))
        {   
            if (output_kretsch && output_relkretschcentregrid)
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
                    FILE *fp;
                    sprintf(name,"%st_relkretschcentregrid.txt",AMRD_save_tag);
                    fp = fopen(name, "a+");
                    fprintf(fp,"%24.16e %24.16e \n",ct,*relkretschcentregrid0);
                    fclose(fp);
                }
            } 

            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder1,          maxquasiset_tt0_freepts_extraporder1,          basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder1,        maxquasiset_tchi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder1,         maxquasiset_txi0_freepts_extraporder1,         basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder1,      maxquasiset_chichi0_freepts_extraporder1,      basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder1,       maxquasiset_chixi0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder1,        maxquasiset_xixi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder1,       maxquasiset_trace0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder1, maxquasiset_massdensity0_freepts_extraporder1, basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder1,               maxbdyphi0_freepts_extraporder1,               basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);    
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder1,          minquasiset_tt0_freepts_extraporder1,          basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder1,        minquasiset_tchi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder1,         minquasiset_txi0_freepts_extraporder1,         basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder1,      minquasiset_chichi0_freepts_extraporder1,      basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder1,       minquasiset_chixi0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder1,        minquasiset_xixi0_freepts_extraporder1,        basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder1,       minquasiset_trace0_freepts_extraporder1,       basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder1, minquasiset_massdensity0_freepts_extraporder1, basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder1,               minbdyphi0_freepts_extraporder1,               basenumbdypoints_freepts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);    
                       for (i=0; i<basenumbdypoints_freepts_extraporder1; i++)
                        { 
                            if (uniSize>1)
                            {
                                quasiset_tt0_freepts_extraporder1          [i] = maxquasiset_tt0_freepts_extraporder1          [i] + minquasiset_tt0_freepts_extraporder1          [i];
                                quasiset_tchi0_freepts_extraporder1        [i] = maxquasiset_tchi0_freepts_extraporder1        [i] + minquasiset_tchi0_freepts_extraporder1        [i];
                                quasiset_txi0_freepts_extraporder1         [i] = maxquasiset_txi0_freepts_extraporder1         [i] + minquasiset_txi0_freepts_extraporder1         [i];
                                quasiset_chichi0_freepts_extraporder1      [i] = maxquasiset_chichi0_freepts_extraporder1      [i] + minquasiset_chichi0_freepts_extraporder1      [i];
                                quasiset_chixi0_freepts_extraporder1       [i] = maxquasiset_chixi0_freepts_extraporder1       [i] + minquasiset_chixi0_freepts_extraporder1       [i];
                                quasiset_xixi0_freepts_extraporder1        [i] = maxquasiset_xixi0_freepts_extraporder1        [i] + minquasiset_xixi0_freepts_extraporder1        [i];
                                quasiset_trace0_freepts_extraporder1       [i] = maxquasiset_trace0_freepts_extraporder1       [i] + minquasiset_trace0_freepts_extraporder1       [i];
                                quasiset_massdensity0_freepts_extraporder1 [i] = maxquasiset_massdensity0_freepts_extraporder1 [i] + minquasiset_massdensity0_freepts_extraporder1 [i];
                                bdyphi0_freepts_extraporder1               [i] = maxbdyphi0_freepts_extraporder1               [i] + minbdyphi0_freepts_extraporder1               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_freepts_extraporder1          [i] = maxquasiset_tt0_freepts_extraporder1          [i];
                                quasiset_tchi0_freepts_extraporder1        [i] = maxquasiset_tchi0_freepts_extraporder1        [i];
                                quasiset_txi0_freepts_extraporder1         [i] = maxquasiset_txi0_freepts_extraporder1         [i];
                                quasiset_chichi0_freepts_extraporder1      [i] = maxquasiset_chichi0_freepts_extraporder1      [i];
                                quasiset_chixi0_freepts_extraporder1       [i] = maxquasiset_chixi0_freepts_extraporder1       [i];
                                quasiset_xixi0_freepts_extraporder1        [i] = maxquasiset_xixi0_freepts_extraporder1        [i];
                                quasiset_trace0_freepts_extraporder1       [i] = maxquasiset_trace0_freepts_extraporder1       [i];
                                quasiset_massdensity0_freepts_extraporder1 [i] = maxquasiset_massdensity0_freepts_extraporder1 [i];
                                bdyphi0_freepts_extraporder1               [i] = maxbdyphi0_freepts_extraporder1               [i];
                            }  
                        }       

                        if (my_rank==0)
                        {     
                            FILE *fp;  
    
                            if (alltimes_ascii)
                            {     
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
                                    
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }        
                                    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }      
        
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
        
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }
    
                            } //closes if(alltimes_ascii) condition 
  
                            if (timestep_ascii)
                            {     
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder1[j],yextrap0_freepts_extraporder1[j],zextrap0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder1[j],yextrap0_freepts_extraporder1[j],zextrap0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);  
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_freepts_extraporder1[j],quasiset_tchi0_freepts_extraporder1[j],quasiset_txi0_freepts_extraporder1[j],quasiset_chichi0_freepts_extraporder1[j],quasiset_chixi0_freepts_extraporder1[j],quasiset_xixi0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                                    sprintf(name,"AdSbdy_freepts_extraporder1_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp);   
    
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_freepts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_freepts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_freepts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }     
    
                            } //closes if (timestep_ascii) condition          
                        } //closes if (my_rank==0) condition            //the following bit computes and prints AdS_mass0 (see below) if we're running on only 1 process
                    
                        if (output_AdS_mass)
                        { 
                            if (uniSize>1)
                            {
                                if (my_rank==0) printf("\nTHE COMPUTATION OF AdS MASS ON MORE THAN 1 PROCESS IS NOT RELIABLE...NOT COMPUTING AdS MASS\n");
                            }
                            else //i.e. we're running on only 1 process
                            {
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FREE POINTS, FIRST ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");      
                                *rhobdy0_freepts_extraporder1=1;
                                chibdy_xibdy_(chibdy0_freepts_extraporder1,xibdy0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,chiextrap0_freepts_extraporder1,xiextrap0_freepts_extraporder1,&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1);    
                                doubleintegralonsphere_(AdS_mass0_freepts_extraporder1,quasiset_massdensity0_freepts_extraporder1,xextrap0_freepts_extraporder1,yextrap0_freepts_extraporder1,zextrap0_freepts_extraporder1,&basenumbdypoints_freepts_extraporder1,rhobdy0_freepts_extraporder1,chibdy0_freepts_extraporder1,xibdy0_freepts_extraporder1,&basebdy_Nchi_freepts_extraporder1,&basebdy_Nxi_freepts_extraporder1);          
                                FILE *fp;
                                sprintf(name,"AdSbdy_freepts_extraporder1_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_freepts_extraporder1);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder2,          maxquasiset_tt0_freepts_extraporder2,          basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder2,        maxquasiset_tchi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder2,         maxquasiset_txi0_freepts_extraporder2,         basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder2,      maxquasiset_chichi0_freepts_extraporder2,      basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder2,       maxquasiset_chixi0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder2,        maxquasiset_xixi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder2,       maxquasiset_trace0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder2, maxquasiset_massdensity0_freepts_extraporder2, basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder2,               maxbdyphi0_freepts_extraporder2,               basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 
                        MPI_Allreduce(lquasiset_tt0_freepts_extraporder2,          minquasiset_tt0_freepts_extraporder2,          basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_freepts_extraporder2,        minquasiset_tchi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_freepts_extraporder2,         minquasiset_txi0_freepts_extraporder2,         basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_freepts_extraporder2,      minquasiset_chichi0_freepts_extraporder2,      basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_freepts_extraporder2,       minquasiset_chixi0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_freepts_extraporder2,        minquasiset_xixi0_freepts_extraporder2,        basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_freepts_extraporder2,       minquasiset_trace0_freepts_extraporder2,       basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_freepts_extraporder2, minquasiset_massdensity0_freepts_extraporder2, basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_freepts_extraporder2,               minbdyphi0_freepts_extraporder2,               basenumbdypoints_freepts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
                        for (i=0; i<basenumbdypoints_freepts_extraporder2; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_freepts_extraporder2          [i] = maxquasiset_tt0_freepts_extraporder2          [i] + minquasiset_tt0_freepts_extraporder2          [i];
                                quasiset_tchi0_freepts_extraporder2        [i] = maxquasiset_tchi0_freepts_extraporder2        [i] + minquasiset_tchi0_freepts_extraporder2        [i];
                                quasiset_txi0_freepts_extraporder2         [i] = maxquasiset_txi0_freepts_extraporder2         [i] + minquasiset_txi0_freepts_extraporder2         [i];
                                quasiset_chichi0_freepts_extraporder2      [i] = maxquasiset_chichi0_freepts_extraporder2      [i] + minquasiset_chichi0_freepts_extraporder2      [i];
                                quasiset_chixi0_freepts_extraporder2       [i] = maxquasiset_chixi0_freepts_extraporder2       [i] + minquasiset_chixi0_freepts_extraporder2       [i];
                                quasiset_xixi0_freepts_extraporder2        [i] = maxquasiset_xixi0_freepts_extraporder2        [i] + minquasiset_xixi0_freepts_extraporder2        [i];
                                quasiset_trace0_freepts_extraporder2       [i] = maxquasiset_trace0_freepts_extraporder2       [i] + minquasiset_trace0_freepts_extraporder2       [i];
                                quasiset_massdensity0_freepts_extraporder2 [i] = maxquasiset_massdensity0_freepts_extraporder2 [i] + minquasiset_massdensity0_freepts_extraporder2 [i];
                                bdyphi0_freepts_extraporder2               [i] = maxbdyphi0_freepts_extraporder2               [i] + minbdyphi0_freepts_extraporder2               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_freepts_extraporder2          [i] = maxquasiset_tt0_freepts_extraporder2          [i];
                                quasiset_tchi0_freepts_extraporder2        [i] = maxquasiset_tchi0_freepts_extraporder2        [i];
                                quasiset_txi0_freepts_extraporder2         [i] = maxquasiset_txi0_freepts_extraporder2         [i];
                                quasiset_chichi0_freepts_extraporder2      [i] = maxquasiset_chichi0_freepts_extraporder2      [i];
                                quasiset_chixi0_freepts_extraporder2       [i] = maxquasiset_chixi0_freepts_extraporder2       [i];
                                quasiset_xixi0_freepts_extraporder2        [i] = maxquasiset_xixi0_freepts_extraporder2        [i];
                                quasiset_trace0_freepts_extraporder2       [i] = maxquasiset_trace0_freepts_extraporder2       [i];
                                quasiset_massdensity0_freepts_extraporder2 [i] = maxquasiset_massdensity0_freepts_extraporder2 [i];
                                bdyphi0_freepts_extraporder2               [i] = maxbdyphi0_freepts_extraporder2               [i];
                            }
                        }   

                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                // save quasiset_ll as ascii
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder2[j],yextrap0_freepts_extraporder2[j],zextrap0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_freepts_extraporder2[j],yextrap0_freepts_extraporder2[j],zextrap0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_freepts_extraporder2[j],j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                // save quasiset_ll as ascii
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                j_red=0;
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_tt0_freepts_extraporder2[j],quasiset_tchi0_freepts_extraporder2[j],quasiset_txi0_freepts_extraporder2[j],quasiset_chichi0_freepts_extraporder2[j],quasiset_chixi0_freepts_extraporder2[j],quasiset_xixi0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_trace0_freepts_extraporder2[j],
                                                j_red);
                                        j_red=j_red+1;
                                    }
                                }
                                fclose(fp);
                            }   
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
                                                j);
                                }
                                fclose(fp); 
                            if (reduced_ascii)
                            {
                                sprintf(name,"AdSbdy_freepts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                fp = fopen(name, "w+");
                                for( j = 0; j < basenumbdypoints_freepts_extraporder2; j++ )
                                {
                                    if ((j%reduction_factor_ascii)==0)
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                ct,
                                                quasiset_massdensity0_freepts_extraporder2[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FREE POINTS, SECOND ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                *rhobdy0_freepts_extraporder2=1;
                                chibdy_xibdy_(chibdy0_freepts_extraporder2,xibdy0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,chiextrap0_freepts_extraporder2,xiextrap0_freepts_extraporder2,&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2);
                                doubleintegralonsphere_(AdS_mass0_freepts_extraporder2,quasiset_massdensity0_freepts_extraporder2,xextrap0_freepts_extraporder2,yextrap0_freepts_extraporder2,zextrap0_freepts_extraporder2,&basenumbdypoints_freepts_extraporder2,rhobdy0_freepts_extraporder2,chibdy0_freepts_extraporder2,xibdy0_freepts_extraporder2,&basebdy_Nchi_freepts_extraporder2,&basebdy_Nxi_freepts_extraporder2); 
                                FILE *fp;
                                sprintf(name,"AdSbdy_freepts_extraporder2_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_freepts_extraporder2);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_freepts_extrap  

                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {
                        bdy_extrap_order=1; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder1,          maxquasiset_tt0_fixedpts_extraporder1,          basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder1,        maxquasiset_tchi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder1,         maxquasiset_txi0_fixedpts_extraporder1,         basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder1,      maxquasiset_chichi0_fixedpts_extraporder1,      basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder1,       maxquasiset_chixi0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder1,        maxquasiset_xixi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder1,       maxquasiset_trace0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder1, maxquasiset_massdensity0_fixedpts_extraporder1, basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder1,               maxbdyphi0_fixedpts_extraporder1,               basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder1,          minquasiset_tt0_fixedpts_extraporder1,          basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder1,        minquasiset_tchi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder1,         minquasiset_txi0_fixedpts_extraporder1,         basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder1,      minquasiset_chichi0_fixedpts_extraporder1,      basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder1,       minquasiset_chixi0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder1,        minquasiset_xixi0_fixedpts_extraporder1,        basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder1,       minquasiset_trace0_fixedpts_extraporder1,       basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder1, minquasiset_massdensity0_fixedpts_extraporder1, basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder1,               minbdyphi0_fixedpts_extraporder1,               basenumbdypoints_fixedpts_extraporder1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
                        for (i=0; i<basenumbdypoints_fixedpts_extraporder1; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_fixedpts_extraporder1          [i] = maxquasiset_tt0_fixedpts_extraporder1          [i] + minquasiset_tt0_fixedpts_extraporder1          [i];
                                quasiset_tchi0_fixedpts_extraporder1        [i] = maxquasiset_tchi0_fixedpts_extraporder1        [i] + minquasiset_tchi0_fixedpts_extraporder1        [i];
                                quasiset_txi0_fixedpts_extraporder1         [i] = maxquasiset_txi0_fixedpts_extraporder1         [i] + minquasiset_txi0_fixedpts_extraporder1         [i];
                                quasiset_chichi0_fixedpts_extraporder1      [i] = maxquasiset_chichi0_fixedpts_extraporder1      [i] + minquasiset_chichi0_fixedpts_extraporder1      [i];
                                quasiset_chixi0_fixedpts_extraporder1       [i] = maxquasiset_chixi0_fixedpts_extraporder1       [i] + minquasiset_chixi0_fixedpts_extraporder1       [i];
                                quasiset_xixi0_fixedpts_extraporder1        [i] = maxquasiset_xixi0_fixedpts_extraporder1        [i] + minquasiset_xixi0_fixedpts_extraporder1        [i];
                                quasiset_trace0_fixedpts_extraporder1       [i] = maxquasiset_trace0_fixedpts_extraporder1       [i] + minquasiset_trace0_fixedpts_extraporder1       [i];
                                quasiset_massdensity0_fixedpts_extraporder1 [i] = maxquasiset_massdensity0_fixedpts_extraporder1 [i] + minquasiset_massdensity0_fixedpts_extraporder1 [i];
                                bdyphi0_fixedpts_extraporder1               [i] = maxbdyphi0_fixedpts_extraporder1               [i] + minbdyphi0_fixedpts_extraporder1               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_fixedpts_extraporder1          [i] = maxquasiset_tt0_fixedpts_extraporder1          [i];
                                quasiset_tchi0_fixedpts_extraporder1        [i] = maxquasiset_tchi0_fixedpts_extraporder1        [i];
                                quasiset_txi0_fixedpts_extraporder1         [i] = maxquasiset_txi0_fixedpts_extraporder1         [i];
                                quasiset_chichi0_fixedpts_extraporder1      [i] = maxquasiset_chichi0_fixedpts_extraporder1      [i];
                                quasiset_chixi0_fixedpts_extraporder1       [i] = maxquasiset_chixi0_fixedpts_extraporder1       [i];
                                quasiset_xixi0_fixedpts_extraporder1        [i] = maxquasiset_xixi0_fixedpts_extraporder1        [i];
                                quasiset_trace0_fixedpts_extraporder1       [i] = maxquasiset_trace0_fixedpts_extraporder1       [i];
                                quasiset_massdensity0_fixedpts_extraporder1 [i] = maxquasiset_massdensity0_fixedpts_extraporder1 [i];
                                bdyphi0_fixedpts_extraporder1               [i] = maxbdyphi0_fixedpts_extraporder1               [i];
                            }
                        }       

                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder1[j],yextrap0_fixedpts_extraporder1[j],zextrap0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder1[j],yextrap0_fixedpts_extraporder1[j],zextrap0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder1[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder1[j],quasiset_tchi0_fixedpts_extraporder1[j],quasiset_txi0_fixedpts_extraporder1[j],quasiset_chichi0_fixedpts_extraporder1[j],quasiset_chixi0_fixedpts_extraporder1[j],quasiset_xixi0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder1[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder1_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder1; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder1[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FIXED POINTS, FIRST ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                *rhobdy0_fixedpts_extraporder1=1;
                                chibdy_xibdy_(chibdy0_fixedpts_extraporder1,xibdy0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,chiextrap0_fixedpts_extraporder1,xiextrap0_fixedpts_extraporder1,&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1);
                                doubleintegralonsphere_(AdS_mass0_fixedpts_extraporder1,quasiset_massdensity0_fixedpts_extraporder1,xextrap0_fixedpts_extraporder1,yextrap0_fixedpts_extraporder1,zextrap0_fixedpts_extraporder1,&basenumbdypoints_fixedpts_extraporder1,rhobdy0_fixedpts_extraporder1,chibdy0_fixedpts_extraporder1,xibdy0_fixedpts_extraporder1,&basebdy_Nchi_fixedpts_extraporder1,&basebdy_Nxi_fixedpts_extraporder1);  
                                FILE *fp;
                                sprintf(name,"AdSbdy_fixedpts_extraporder1_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_fixedpts_extraporder1);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder1  

                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {
                        bdy_extrap_order=2; 
                        // for each n,i point on the outer bdy, save sum{lquasisetll[n,i]}_allprocessors into quasisetll[n,i]
                        //basenumbdypoints is set in AdS4D_post_init
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder2,          maxquasiset_tt0_fixedpts_extraporder2,          basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder2,        maxquasiset_tchi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder2,         maxquasiset_txi0_fixedpts_extraporder2,         basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder2,      maxquasiset_chichi0_fixedpts_extraporder2,      basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder2,       maxquasiset_chixi0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder2,        maxquasiset_xixi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder2,       maxquasiset_trace0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder2, maxquasiset_massdensity0_fixedpts_extraporder2, basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder2,               maxbdyphi0_fixedpts_extraporder2,               basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);  
                        MPI_Allreduce(lquasiset_tt0_fixedpts_extraporder2,          minquasiset_tt0_fixedpts_extraporder2,          basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_tchi0_fixedpts_extraporder2,        minquasiset_tchi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_txi0_fixedpts_extraporder2,         minquasiset_txi0_fixedpts_extraporder2,         basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chichi0_fixedpts_extraporder2,      minquasiset_chichi0_fixedpts_extraporder2,      basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_chixi0_fixedpts_extraporder2,       minquasiset_chixi0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_xixi0_fixedpts_extraporder2,        minquasiset_xixi0_fixedpts_extraporder2,        basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_trace0_fixedpts_extraporder2,       minquasiset_trace0_fixedpts_extraporder2,       basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lquasiset_massdensity0_fixedpts_extraporder2, minquasiset_massdensity0_fixedpts_extraporder2, basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
                        MPI_Allreduce(lbdyphi0_fixedpts_extraporder2,               minbdyphi0_fixedpts_extraporder2,               basenumbdypoints_fixedpts_extraporder2,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);  
                        for (i=0; i<basenumbdypoints_fixedpts_extraporder2; i++)
                        {
                            if (uniSize>1)
                            {
                                quasiset_tt0_fixedpts_extraporder2          [i] = maxquasiset_tt0_fixedpts_extraporder2          [i] + minquasiset_tt0_fixedpts_extraporder2          [i];
                                quasiset_tchi0_fixedpts_extraporder2        [i] = maxquasiset_tchi0_fixedpts_extraporder2        [i] + minquasiset_tchi0_fixedpts_extraporder2        [i];
                                quasiset_txi0_fixedpts_extraporder2         [i] = maxquasiset_txi0_fixedpts_extraporder2         [i] + minquasiset_txi0_fixedpts_extraporder2         [i];
                                quasiset_chichi0_fixedpts_extraporder2      [i] = maxquasiset_chichi0_fixedpts_extraporder2      [i] + minquasiset_chichi0_fixedpts_extraporder2      [i];
                                quasiset_chixi0_fixedpts_extraporder2       [i] = maxquasiset_chixi0_fixedpts_extraporder2       [i] + minquasiset_chixi0_fixedpts_extraporder2       [i];
                                quasiset_xixi0_fixedpts_extraporder2        [i] = maxquasiset_xixi0_fixedpts_extraporder2        [i] + minquasiset_xixi0_fixedpts_extraporder2        [i];
                                quasiset_trace0_fixedpts_extraporder2       [i] = maxquasiset_trace0_fixedpts_extraporder2       [i] + minquasiset_trace0_fixedpts_extraporder2       [i];
                                quasiset_massdensity0_fixedpts_extraporder2 [i] = maxquasiset_massdensity0_fixedpts_extraporder2 [i] + minquasiset_massdensity0_fixedpts_extraporder2 [i];
                                bdyphi0_fixedpts_extraporder2               [i] = maxbdyphi0_fixedpts_extraporder2               [i] + minbdyphi0_fixedpts_extraporder2               [i];
                            }
                            else //if uniSize==1, i.e. there is only 1 process, maxquasiset=minquasiset so we have to take only one of them into consideration
                            {
                                quasiset_tt0_fixedpts_extraporder2          [i] = maxquasiset_tt0_fixedpts_extraporder2          [i];
                                quasiset_tchi0_fixedpts_extraporder2        [i] = maxquasiset_tchi0_fixedpts_extraporder2        [i];
                                quasiset_txi0_fixedpts_extraporder2         [i] = maxquasiset_txi0_fixedpts_extraporder2         [i];
                                quasiset_chichi0_fixedpts_extraporder2      [i] = maxquasiset_chichi0_fixedpts_extraporder2      [i];
                                quasiset_chixi0_fixedpts_extraporder2       [i] = maxquasiset_chixi0_fixedpts_extraporder2       [i];
                                quasiset_xixi0_fixedpts_extraporder2        [i] = maxquasiset_xixi0_fixedpts_extraporder2        [i];
                                quasiset_trace0_fixedpts_extraporder2       [i] = maxquasiset_trace0_fixedpts_extraporder2       [i];
                                quasiset_massdensity0_fixedpts_extraporder2 [i] = maxquasiset_massdensity0_fixedpts_extraporder2 [i];
                                bdyphi0_fixedpts_extraporder2               [i] = maxbdyphi0_fixedpts_extraporder2               [i];
                            }
                        }   
                        if (my_rank==0)
                        {   
                            FILE *fp;   
                            if (alltimes_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_bdyphi_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetll_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisettrace_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint.txt",AMRD_save_tag);
                                    fp = fopen(name, "a+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                            } //closes if(alltimes_ascii) condition 
                            if (timestep_ascii)
                            {   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder2[j],yextrap0_fixedpts_extraporder2[j],zextrap0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_xext_yext_zext_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %i \n",ct,xextrap0_fixedpts_extraporder2[j],yextrap0_fixedpts_extraporder2[j],zextrap0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_bdyphi_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",ct,bdyphi0_fixedpts_extraporder2[j],j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    // save quasiset_ll as ascii
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetll_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    j_red=0;
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %24.16e %24.16e  %24.16e %24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_tt0_fixedpts_extraporder2[j],quasiset_tchi0_fixedpts_extraporder2[j],quasiset_txi0_fixedpts_extraporder2[j],quasiset_chichi0_fixedpts_extraporder2[j],quasiset_chixi0_fixedpts_extraporder2[j],quasiset_xixi0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisettrace_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_trace0_fixedpts_extraporder2[j],
                                                    j_red);
                                            j_red=j_red+1;
                                        }
                                    }
                                    fclose(fp);
                                }   
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
                                                    j);
                                    }
                                    fclose(fp); 
                                if (reduced_ascii)
                                {
                                    sprintf(name,"AdSbdy_fixedpts_extraporder2_reduced_%st_quasisetmassdensity_indbdypoint_tstep%d.txt",AMRD_save_tag,lsteps);
                                    fp = fopen(name, "w+");
                                    for( j = 0; j < basenumbdypoints_fixedpts_extraporder2; j++ )
                                    {
                                        if ((j%reduction_factor_ascii)==0)
                                        {
                                            fprintf(fp,"%24.16e %24.16e %i \n",
                                                    ct,
                                                    quasiset_massdensity0_fixedpts_extraporder2[j],
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
                                printf("\nRUNNING ON ONLY 1 PROCESS...THE NUMERICAL APPROXIMATION OF AdS MASS FOR FIXED POINTS, SECOND ORDER EXTRAPOLATION IS RELIABLE ON 1 PROCESS...COMPUTING AdS MASS\n");    
                                chibdy_xibdy_(chibdy0_fixedpts_extraporder2,xibdy0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,chiextrap0_fixedpts_extraporder2,xiextrap0_fixedpts_extraporder2,&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2);
                                doubleintegralonsphere_(AdS_mass0_fixedpts_extraporder2,quasiset_massdensity0_fixedpts_extraporder2,xextrap0_fixedpts_extraporder2,yextrap0_fixedpts_extraporder2,zextrap0_fixedpts_extraporder2,&basenumbdypoints_fixedpts_extraporder2,rhobdy0_fixedpts_extraporder2,chibdy0_fixedpts_extraporder2,xibdy0_fixedpts_extraporder2,&basebdy_Nchi_fixedpts_extraporder2,&basebdy_Nxi_fixedpts_extraporder2);  
                                FILE *fp;
                                sprintf(name,"AdSbdy_fixedpts_extraporder2_%st_AdSmass.txt",AMRD_save_tag);
                                fp = fopen(name, "a+");
                                    fprintf(fp,"%24.16e %24.16e \n",ct,*AdS_mass0_fixedpts_extraporder2);
                                fclose(fp);
                            }
                        }   
                    }//closes condition on output_bdy_extraporder2  
                }//closes condition on bdy_fixedpts_extrap  
            }//closes output_bdyquantities if-condition
        }//closes if condition on lsteps
    } //closes if condition on L==Lc

    //free memory for boundary quantities
    if (L==Lc)
    {   
        valid=PAMR_init_s_iter(L,PAMR_AMRH,0);
        while(valid)
        {   
            if (output_bdyquantities)
            {   
                //FREE POINTS EXTRAPOLATION
                if (bdy_extrap_freepts)
                {
                    //FREE POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {      
                        free(vecbdypoints_freepts_extraporder1);
                        free(dsplsbdypoints_freepts_extraporder1);
                        free(quasiset_tt_freepts_extraporder1);
                        free(quasiset_tchi_freepts_extraporder1);
                        free(quasiset_txi_freepts_extraporder1);
                        free(quasiset_chichi_freepts_extraporder1);
                        free(quasiset_chixi_freepts_extraporder1);
                        free(quasiset_xixi_freepts_extraporder1);
                        free(quasiset_trace_freepts_extraporder1);
                        free(quasiset_massdensity_freepts_extraporder1);
                        free(bdyphi_freepts_extraporder1);  
                        free(xextrap_freepts_extraporder1);
                        free(yextrap_freepts_extraporder1);
                        free(zextrap_freepts_extraporder1);   
                        free(lquasiset_tt0_freepts_extraporder1);        
                        free(lquasiset_tchi0_freepts_extraporder1);
                        free(lquasiset_txi0_freepts_extraporder1);
                        free(lquasiset_chichi0_freepts_extraporder1);
                        free(lquasiset_chixi0_freepts_extraporder1);
                        free(lquasiset_xixi0_freepts_extraporder1);
                        free(lquasiset_trace0_freepts_extraporder1);
                        free(lquasiset_massdensity0_freepts_extraporder1);
                        free(lbdyphi0_freepts_extraporder1);      
                        free(maxquasiset_tt0_freepts_extraporder1);        
                        free(maxquasiset_tchi0_freepts_extraporder1);
                        free(maxquasiset_txi0_freepts_extraporder1);
                        free(maxquasiset_chichi0_freepts_extraporder1);
                        free(maxquasiset_chixi0_freepts_extraporder1);
                        free(maxquasiset_xixi0_freepts_extraporder1);
                        free(maxquasiset_trace0_freepts_extraporder1);
                        free(maxquasiset_massdensity0_freepts_extraporder1);
                        free(maxbdyphi0_freepts_extraporder1);    
                        free(minquasiset_tt0_freepts_extraporder1);
                        free(minquasiset_tchi0_freepts_extraporder1);
                        free(minquasiset_txi0_freepts_extraporder1);
                        free(minquasiset_chichi0_freepts_extraporder1);
                        free(minquasiset_chixi0_freepts_extraporder1);
                        free(minquasiset_xixi0_freepts_extraporder1);
                        free(minquasiset_trace0_freepts_extraporder1);
                        free(minquasiset_massdensity0_freepts_extraporder1);
                        free(minbdyphi0_freepts_extraporder1);    
                        free(quasiset_tt0_freepts_extraporder1);
                        free(quasiset_tchi0_freepts_extraporder1);
                        free(quasiset_txi0_freepts_extraporder1);
                        free(quasiset_chichi0_freepts_extraporder1);
                        free(quasiset_chixi0_freepts_extraporder1);
                        free(quasiset_xixi0_freepts_extraporder1);
                        free(quasiset_trace0_freepts_extraporder1);
                        free(quasiset_massdensity0_freepts_extraporder1);
                        free(bdyphi0_freepts_extraporder1);
                        free(AdS_mass0_freepts_extraporder1);     
                        free(xextrap0_freepts_extraporder1);
                        free(yextrap0_freepts_extraporder1);
                        free(zextrap0_freepts_extraporder1);          
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_freepts_extraporder1);
                            free(chiextrap0_freepts_extraporder1);
                            free(xiextrap0_freepts_extraporder1);             
                            free(rhobdy0_freepts_extraporder1);
                            free(chibdy0_freepts_extraporder1);
                            free(xibdy0_freepts_extraporder1);    
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FREE POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {   
                        free(vecbdypoints_freepts_extraporder2);
                        free(dsplsbdypoints_freepts_extraporder2);
                        free(quasiset_tt_freepts_extraporder2);
                        free(quasiset_tchi_freepts_extraporder2);
                        free(quasiset_txi_freepts_extraporder2);
                        free(quasiset_chichi_freepts_extraporder2);
                        free(quasiset_chixi_freepts_extraporder2);
                        free(quasiset_xixi_freepts_extraporder2);
                        free(quasiset_trace_freepts_extraporder2);
                        free(quasiset_massdensity_freepts_extraporder2);
                        free(bdyphi_freepts_extraporder2);  
                        free(xextrap_freepts_extraporder2);
                        free(yextrap_freepts_extraporder2);
                        free(zextrap_freepts_extraporder2); 
                        free(lquasiset_tt0_freepts_extraporder2);
                        free(lquasiset_tchi0_freepts_extraporder2);
                        free(lquasiset_txi0_freepts_extraporder2);
                        free(lquasiset_chichi0_freepts_extraporder2);
                        free(lquasiset_chixi0_freepts_extraporder2);
                        free(lquasiset_xixi0_freepts_extraporder2);
                        free(lquasiset_trace0_freepts_extraporder2);
                        free(lquasiset_massdensity0_freepts_extraporder2);
                        free(lbdyphi0_freepts_extraporder2);    
                        free(maxquasiset_tt0_freepts_extraporder2);
                        free(maxquasiset_tchi0_freepts_extraporder2);
                        free(maxquasiset_txi0_freepts_extraporder2);
                        free(maxquasiset_chichi0_freepts_extraporder2);
                        free(maxquasiset_chixi0_freepts_extraporder2);
                        free(maxquasiset_xixi0_freepts_extraporder2);
                        free(maxquasiset_trace0_freepts_extraporder2);
                        free(maxquasiset_massdensity0_freepts_extraporder2);
                        free(maxbdyphi0_freepts_extraporder2);  
                        free(minquasiset_tt0_freepts_extraporder2);
                        free(minquasiset_tchi0_freepts_extraporder2);
                        free(minquasiset_txi0_freepts_extraporder2);
                        free(minquasiset_chichi0_freepts_extraporder2);
                        free(minquasiset_chixi0_freepts_extraporder2);
                        free(minquasiset_xixi0_freepts_extraporder2);
                        free(minquasiset_trace0_freepts_extraporder2);
                        free(minquasiset_massdensity0_freepts_extraporder2);
                        free(minbdyphi0_freepts_extraporder2);  
                        free(quasiset_tt0_freepts_extraporder2);
                        free(quasiset_tchi0_freepts_extraporder2);
                        free(quasiset_txi0_freepts_extraporder2);
                        free(quasiset_chichi0_freepts_extraporder2);
                        free(quasiset_chixi0_freepts_extraporder2);
                        free(quasiset_xixi0_freepts_extraporder2);
                        free(quasiset_trace0_freepts_extraporder2);
                        free(quasiset_massdensity0_freepts_extraporder2);
                        free(bdyphi0_freepts_extraporder2);
                        free(AdS_mass0_freepts_extraporder2);   
                        free(xextrap0_freepts_extraporder2);
                        free(yextrap0_freepts_extraporder2);
                        free(zextrap0_freepts_extraporder2);    
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_freepts_extraporder2);
                            free(chiextrap0_freepts_extraporder2);
                            free(xiextrap0_freepts_extraporder2);   
                            free(rhobdy0_freepts_extraporder2);
                            free(chibdy0_freepts_extraporder2);
                            free(xibdy0_freepts_extraporder2);  
                        }   
                    } //closes condition on output_bdy_extraporder2 
                }//closes condition on bdy_freepts_extrap       
                //FIXED POINTS EXTRAPOLATION
                if (bdy_extrap_fixedpts)
                {
                    //FIXED POINTS, FIRST ORDER EXTRAPOLATION
                    if (output_bdy_extraporder1)
                    {   
                        free(vecbdypoints_fixedpts_extraporder1);
                        free(dsplsbdypoints_fixedpts_extraporder1);
                        free(quasiset_tt_fixedpts_extraporder1);
                        free(quasiset_tchi_fixedpts_extraporder1);
                        free(quasiset_txi_fixedpts_extraporder1);
                        free(quasiset_chichi_fixedpts_extraporder1);
                        free(quasiset_chixi_fixedpts_extraporder1);
                        free(quasiset_xixi_fixedpts_extraporder1);
                        free(quasiset_trace_fixedpts_extraporder1);
                        free(quasiset_massdensity_fixedpts_extraporder1);
                        free(bdyphi_fixedpts_extraporder1); 
                        free(xextrap_fixedpts_extraporder1);
                        free(yextrap_fixedpts_extraporder1);
                        free(zextrap_fixedpts_extraporder1);    
                        free(lquasiset_tt0_fixedpts_extraporder1);
                        free(lquasiset_tchi0_fixedpts_extraporder1);
                        free(lquasiset_txi0_fixedpts_extraporder1);
                        free(lquasiset_chichi0_fixedpts_extraporder1);
                        free(lquasiset_chixi0_fixedpts_extraporder1);
                        free(lquasiset_xixi0_fixedpts_extraporder1);
                        free(lquasiset_trace0_fixedpts_extraporder1);
                        free(lquasiset_massdensity0_fixedpts_extraporder1);
                        free(lbdyphi0_fixedpts_extraporder1);   
                        free(maxquasiset_tt0_fixedpts_extraporder1);
                        free(maxquasiset_tchi0_fixedpts_extraporder1);
                        free(maxquasiset_txi0_fixedpts_extraporder1);
                        free(maxquasiset_chichi0_fixedpts_extraporder1);
                        free(maxquasiset_chixi0_fixedpts_extraporder1);
                        free(maxquasiset_xixi0_fixedpts_extraporder1);
                        free(maxquasiset_trace0_fixedpts_extraporder1);
                        free(maxquasiset_massdensity0_fixedpts_extraporder1);
                        free(maxbdyphi0_fixedpts_extraporder1); 
                        free(minquasiset_tt0_fixedpts_extraporder1);
                        free(minquasiset_tchi0_fixedpts_extraporder1);
                        free(minquasiset_txi0_fixedpts_extraporder1);
                        free(minquasiset_chichi0_fixedpts_extraporder1);
                        free(minquasiset_chixi0_fixedpts_extraporder1);
                        free(minquasiset_xixi0_fixedpts_extraporder1);
                        free(minquasiset_trace0_fixedpts_extraporder1);
                        free(minquasiset_massdensity0_fixedpts_extraporder1);
                        free(minbdyphi0_fixedpts_extraporder1); 
                        free(quasiset_tt0_fixedpts_extraporder1);
                        free(quasiset_tchi0_fixedpts_extraporder1);
                        free(quasiset_txi0_fixedpts_extraporder1);
                        free(quasiset_chichi0_fixedpts_extraporder1);
                        free(quasiset_chixi0_fixedpts_extraporder1);
                        free(quasiset_xixi0_fixedpts_extraporder1);
                        free(quasiset_trace0_fixedpts_extraporder1);
                        free(quasiset_massdensity0_fixedpts_extraporder1);
                        free(bdyphi0_fixedpts_extraporder1);
                        free(AdS_mass0_fixedpts_extraporder1);  
                        free(xextrap0_fixedpts_extraporder1);
                        free(yextrap0_fixedpts_extraporder1);
                        free(zextrap0_fixedpts_extraporder1);   
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_fixedpts_extraporder1);
                            free(chiextrap0_fixedpts_extraporder1);
                            free(xiextrap0_fixedpts_extraporder1);  
                            free(rhobdy0_fixedpts_extraporder1);
                            free(chibdy0_fixedpts_extraporder1);
                            free(xibdy0_fixedpts_extraporder1); 
                        }   
                    }//closes condition on output_bdy_extraporder1  
                    //FIXED POINTS, SECOND ORDER EXTRAPOLATION
                    if (output_bdy_extraporder2)
                    {   
                        free(vecbdypoints_fixedpts_extraporder2);
                        free(dsplsbdypoints_fixedpts_extraporder2);
                        free(quasiset_tt_fixedpts_extraporder2);
                        free(quasiset_tchi_fixedpts_extraporder2);
                        free(quasiset_txi_fixedpts_extraporder2);
                        free(quasiset_chichi_fixedpts_extraporder2);
                        free(quasiset_chixi_fixedpts_extraporder2);
                        free(quasiset_xixi_fixedpts_extraporder2);
                        free(quasiset_trace_fixedpts_extraporder2);
                        free(quasiset_massdensity_fixedpts_extraporder2);
                        free(bdyphi_fixedpts_extraporder2); 
                        free(xextrap_fixedpts_extraporder2);
                        free(yextrap_fixedpts_extraporder2);
                        free(zextrap_fixedpts_extraporder2);    
                        free(lquasiset_tt0_fixedpts_extraporder2);
                        free(lquasiset_tchi0_fixedpts_extraporder2);
                        free(lquasiset_txi0_fixedpts_extraporder2);
                        free(lquasiset_chichi0_fixedpts_extraporder2);
                        free(lquasiset_chixi0_fixedpts_extraporder2);
                        free(lquasiset_xixi0_fixedpts_extraporder2);
                        free(lquasiset_trace0_fixedpts_extraporder2);
                        free(lquasiset_massdensity0_fixedpts_extraporder2);
                        free(lbdyphi0_fixedpts_extraporder2);   
                        free(maxquasiset_tt0_fixedpts_extraporder2);
                        free(maxquasiset_tchi0_fixedpts_extraporder2);
                        free(maxquasiset_txi0_fixedpts_extraporder2);
                        free(maxquasiset_chichi0_fixedpts_extraporder2);
                        free(maxquasiset_chixi0_fixedpts_extraporder2);
                        free(maxquasiset_xixi0_fixedpts_extraporder2);
                        free(maxquasiset_trace0_fixedpts_extraporder2);
                        free(maxquasiset_massdensity0_fixedpts_extraporder2);
                        free(maxbdyphi0_fixedpts_extraporder2); 
                        free(minquasiset_tt0_fixedpts_extraporder2);
                        free(minquasiset_tchi0_fixedpts_extraporder2);
                        free(minquasiset_txi0_fixedpts_extraporder2);
                        free(minquasiset_chichi0_fixedpts_extraporder2);
                        free(minquasiset_chixi0_fixedpts_extraporder2);
                        free(minquasiset_xixi0_fixedpts_extraporder2);
                        free(minquasiset_trace0_fixedpts_extraporder2);
                        free(minquasiset_massdensity0_fixedpts_extraporder2);
                        free(minbdyphi0_fixedpts_extraporder2); 
                        free(quasiset_tt0_fixedpts_extraporder2);
                        free(quasiset_tchi0_fixedpts_extraporder2);
                        free(quasiset_txi0_fixedpts_extraporder2);
                        free(quasiset_chichi0_fixedpts_extraporder2);
                        free(quasiset_chixi0_fixedpts_extraporder2);
                        free(quasiset_xixi0_fixedpts_extraporder2);
                        free(quasiset_trace0_fixedpts_extraporder2);
                        free(quasiset_massdensity0_fixedpts_extraporder2);
                        free(bdyphi0_fixedpts_extraporder2);
                        free(AdS_mass0_fixedpts_extraporder2);  
                        free(xextrap0_fixedpts_extraporder2);
                        free(yextrap0_fixedpts_extraporder2);
                        free(zextrap0_fixedpts_extraporder2);   
                        if (output_AdS_mass)
                        {
                            free(rhoextrap0_fixedpts_extraporder2);
                            free(chiextrap0_fixedpts_extraporder2);
                            free(xiextrap0_fixedpts_extraporder2);  
                            free(rhobdy0_fixedpts_extraporder2);
                            free(chibdy0_fixedpts_extraporder2);
                            free(xibdy0_fixedpts_extraporder2); 
                        }   
                    } //closes condition on output_bdy_extraporder2 
                }//closes condition on bdy_extrap_fixedpts
            } //closes condition on output_bdyquantities
            valid=PAMR_next_g();
        }       
    } //closes condition on (L==Lc)

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

