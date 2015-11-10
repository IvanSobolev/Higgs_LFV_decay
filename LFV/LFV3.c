#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cteqpdf.h>
#include <complex.h>
#include "mynrutil.h"
#include "mynr.h"

#define frand() ((double) rand() / (RAND_MAX+1.0))

extern double ctq6pdf_(int *, double *, double *);
extern int setctq6_(int *);
extern void spline(double *, double *, int, double, double, double *);
extern void splint(double *, double *, double *, int, double, double *);
extern void jacobi(double **a, int n, double d[], double **v, int *nrot);
double integral(double,double,int,int);

double BR(double,double,double,double,double,double,double,double);
double lambda(double,double);
double F1(double,double,double,double);
double F2(double,double,double);
double F3(double,double,double,double);
double F4(double,double,double,double);
double GG1(double,double,int);
double GG2(double,double,int);
double GG3(double,double,int);
/* d -- for "delta-function approximation"*/
double GG1_d(double); 
double GG2_d(double);
double GG3_d(double);
complex f(double);
complex A_Q(double);
complex A_W(double);
complex A_f(double);
double a_s(double);
double c(double);
double c_s(double);
double run_mQ(double,double);
double run_mQ_s(double,double);
double mQ_MS_bar(double); //quark mass in MS-bar scheme


double alpha = 1/127.94; //EM coupling constant on M_Z scale
double m_tau = 1.777; //GeV, tau-lepton mass
double m_mu = 0.106; //GeV, tau-lepton mass
double G_tau = 2.27e-12; //GeV, tau-lepton width
double s2w = 0.23126; //sin^2(teta_w)
double v = 174.0; //Higgs field V.E.V.
double mh = 125.0; //GeV, Higgs mass
double mw = 80.39; //GeV, W-boson mass
double Gw = 2.09; //GeV, W-boson width
double mZ = 91.19;
double GZ = 2.50;
int nf = 5; 
int nf_s = 6;
double GF = 1.166e-05; //Gev^{-2}, Fermi's constant
double BR_max = 4.4e-08; //limit on BR(tau->mu+gamma)
double m_t = 177.1; //GeV, top quark pole mass
double m_b = 4.87; //GeV, bottom quark pole mass
double m_c = 1.64; //GeV, charm quark pole mass


int  main(void) { 
  double a1,a2;
  double G1_W,G2_W,G3_W,G1_Z,G2_Z,G3_Z;
  double G1_W_s,G2_W_s,G3_W_s,G1_Z_s,G2_Z_s,G3_Z_s;
  double delta_V;
  double b_tau,b_tau_s;
  double m_b_Q,m_b_run;
  double tau_t,tau_b,tau_c,tau_W;
  double tau_t_s,tau_b_s,tau_c_s,tau_W_s;
  double MZZ,Mgmgm; //sgoldstino couplings to vector bosons
  double Omega_W,Omega_Z;
  double xi_W,eta_W;
  double xi_Z,eta_Z;
  double a_H,a_b;
  double a_H_s;
  double D_QCD,D_t,D_QCD_s,D_t_s;
  double Gh_gg_SM,Gh_WW_SM,Gh_tau_SM,Gh_ZZ_SM,Gh_bb_SM,Gh_gmgm_SM,Gh_SM; //GeV, partial Higgs widths in SM
  double Br_gg_SM,Br_WW_SM,Br_ZZ_SM,Br_tau_SM,Br_bb_SM,Br_gmgm_SM; //Branching ratios of Higgs in SM
  double Gh_gg_m,Gh_WW_m,Gh_tau_m,Gh_ZZ_m,Gh_bb_m,Gh_gmgm_m,Gh_m; //GeV, partial Higgs widths in our model
  double Br_gg_m,Br_WW_m,Br_ZZ_m,Br_tau_m,Br_bb_m,Br_gmgm_m; //Branching ratios of Higgs in our model
  double mu_gg,mu_WW,mu_ZZ,mu_tau,mu_bb,mu_gmgm;
  double mu_WW_max_A,mu_WW_min_A,mu_WW_max_C,mu_WW_min_C,mu_bb_max,mu_bb_min;
  double mu_ZZ_max_A,mu_ZZ_min_A,mu_ZZ_max_C,mu_ZZ_min_C;
  double mu_gmgm_max_A,mu_gmgm_min_A,mu_gmgm_max_C,mu_gmgm_min_C;
  double mu_tau_max_A,mu_tau_min_A,mu_tau_max_C,mu_tau_min_C;
  double Gsg_gg,Gsg_WW,Gsg_tau,Gsg_ZZ,Gsg_bb,Gsg_gmgm,Gsg_GG,Gsg; //GeV, partial sgoldstino widths
  char c1;
  // we use GeV for energy units
  double mu,beta,M1,M2,M3,X,tn,msL; //msL -- slepton scale
  double mu_min,tn_min,M1_min,M2_min,M3_min,msL_min;
  double mu_max,tn_max,M1_max,M2_max,M3_max,msL_max;  
  double Y,Y_min =0.0008,Y_max = 0.0036; //from hep-ex/1502.07400
  double AF_min,AF_max,AF,AF1,AF2,Amu,Atau; //sqrt(A_{ab}^2+A_{ba}^2)/sqrt(2F)
  double sqF = 4000.0; // sqrt(F)
  double F;
  double ms_min,ms_max;
  double *ms_l,*Br_l,*yd2_l; //limits, ms_l in TeV, Br_l in pb
  double *ms1,*Int_ms,*yd2_ms; //vectors for interpolating function of cross-section
  double yd1,ydN;
  double ms,d_ms; //sgoldstino mass m_s
  double tn2,teta;
  double sqrt_s; //s.o.m. energy
  double s;
  int Np,i,j;
  int Np_rnd = 100000000;
  int N_lim,mode;
  int Nms; //number of points for interpolating function of cross-section
  FILE *output_mu,*output,*input_lim;
  double Int_br,sigma_sg;
  double *Mass_m[5],*vects[5];
  double eigv[5];
  int nrot,flag;
  
  for(i=0;i<=4;i++) {
    Mass_m[i] = (double *) calloc(5,sizeof(double));
    vects[i] = (double *) calloc(5,sizeof(double));
  }
  mu_WW_max_A = 1.31;
  mu_WW_min_A = 0.76;
  mu_WW_max_C = 0.92;
  mu_WW_min_C = 0.54;
  mu_ZZ_max_A = 1.84;
  mu_ZZ_min_A = 1.11;
  mu_ZZ_max_C = 1.32;
  mu_ZZ_min_C = 0.61;
  mu_gmgm_max_A = 1.49;
  mu_gmgm_min_A = 0.8;
  mu_gmgm_max_C = 1.40;
  mu_gmgm_min_C = 0.91;
  mu_tau_max_A = 1.8;
  mu_tau_min_A = 1.0;
  mu_tau_max_C = 1.05;
  mu_tau_min_C = 0.51;
  mu_bb_max = 1.50;
  mu_bb_min = 0.50;
  /* Reading limits on sgoldstino cross_section X Br(sg --> 2*gm)*/
  input_lim = fopen("sg_cs_limit","r");
  //input_lim = fopen("sg_cs_limit_low","r");
  output = fopen("LFV_output","w");
  output_mu = fopen("LFV_output_mu","w");
  //output = fopen("LFV_output_low","w");
  //output_mu = fopen("LFV_output_mu_low","w");
  N_lim = 0;
  while ((c1 = getc(input_lim)) != EOF) {
    if (c1 == '\n')  N_lim++;
  }
  fclose(input_lim);
  ms_l = (double *) calloc(N_lim,sizeof(double));
  Br_l = (double *) calloc(N_lim,sizeof(double));
  yd2_l = (double *) calloc(N_lim,sizeof(double));
  input_lim = fopen("sg_cs_limit","r");
  //input_lim = fopen("sg_cs_limit_low","r");
  for(i=0; i<=N_lim-1; i++) {
    fscanf(input_lim,"%lf %lf\n",ms_l+i,Br_l+i);
  }
  fclose(input_lim);  
  spline(ms_l,Br_l,N_lim,yd1,ydN,yd2_l);  
  ms_min = ms_l[0]*1000.0;
  ms_max = ms_l[N_lim-1]*1000.0;
  /* */
  /*Initializing interpolating function for sgoldstino cross_section*/
  sqrt_s = 8000.0; //s.o.m energy, GeV
  s = sqrt_s*sqrt_s;
  Nms = 400; //number of points
  d_ms = (ms_max-ms_min)/(double)Nms;  
  ms1 = (double *) calloc(Nms+1,sizeof(double));
  Int_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_ms = (double *) calloc(Nms+1,sizeof(double));
  mode = 4;
  Np = 1000;
  for(i=0;i<=Nms;i++) {
    ms1[i] = ms_min+d_ms*i;
    Int_ms[i] = integral(ms1[i],sqrt_s,mode,Np);
  }
  spline(ms1,Int_ms,Nms+1,yd1,ydN,yd2_ms);  
  /* */
  /*Calculating SM Higgs widths*/
  tau_t = 4.0*(m_t*m_t)/(mh*mh);
  tau_b = 4.0*(m_b*m_b)/(mh*mh);
  tau_c = 4.0*(m_c*m_c)/(mh*mh);
  tau_W = 4.0*(mw*mw)/(mh*mh);

  a_H = a_s(mh)/M_PI;
  a_b = a_s(m_b)/M_PI;
  m_b_Q = mQ_MS_bar(m_b);

  xi_W = mw/mh;
  eta_W = Gw/mh;
  xi_Z = mZ/mh;
  eta_Z = GZ/mh;
  Np = 1000;
  G1_W = GG1(xi_W,eta_W,Np);
  G1_Z = GG1(xi_Z,eta_Z,Np);
  G2_W = GG2(xi_W,eta_W,Np);
  G2_Z = GG2(xi_Z,eta_Z,Np);
  G3_W = GG3(xi_W,eta_W,Np);
  G3_Z = GG3(xi_Z,eta_Z,Np);

  // H->WW
  delta_V = 2.0;
  Gh_WW_SM = G1_W*delta_V*GF*mh*mh*mh;
  Gh_WW_SM = Gh_WW_SM/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
  // H->ZZ
  delta_V = 1.0;
  Gh_ZZ_SM = G1_Z*delta_V*GF*mh*mh*mh;
  Gh_ZZ_SM = Gh_ZZ_SM/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
  // H->gg
  Gh_gg_SM = GF*a_s(mh)*a_s(mh)*mh*mh*mh/(36.0*sqrt(2.0)*M_PI*M_PI*M_PI);
  Gh_gg_SM = Gh_gg_SM*pow(cabs(A_Q(tau_t)+A_Q(tau_b)+A_Q(tau_c)),2.0);
  Gh_gg_SM = Gh_gg_SM*(1.0+(95.0/4.0-7.0/6.0*(double)nf)*a_H);
  // H->tau tau
  b_tau = sqrt(1-4.0*m_tau*m_tau/(mh*mh));
  Gh_tau_SM = GF*mh*m_tau*m_tau*b_tau*b_tau*b_tau/(4.0*sqrt(2.0)*M_PI);
  // H->bb
  D_QCD= 1.0+5.67*a_H+(35.94-1.36*(double)nf)*a_H*a_H;
  D_QCD = D_QCD + (164.14-25.77*(double)nf+0.259*(double)nf*(double)nf)*a_H*a_H*a_H;
  D_t = a_H*a_H*(1.57-4.0*log(mh/m_t)/3.0+2.0/9.0*log(run_mQ(m_b,mh)/mh));
  Gh_bb_SM = 3.0*GF*mh*run_mQ(m_b,mh)*run_mQ(m_b,mh)/(4.0*sqrt(2.0)*M_PI)*(D_QCD+D_t);
  // H->gm gm
  Gh_gmgm_SM = pow(cabs(4.0/3.0*A_f(tau_t)+4.0/3.0*A_f(tau_c)+1.0/3.0*A_f(tau_b)+A_W(tau_W)),2.0)*(1-a_H);
  Gh_gmgm_SM = Gh_gmgm_SM*GF*mh*mh*mh*alpha*alpha/(128.0*sqrt(2.0)*M_PI*M_PI*M_PI);
  // total width
  Gh_SM = Gh_WW_SM+Gh_ZZ_SM+Gh_gmgm_SM+Gh_tau_SM+Gh_bb_SM+Gh_gg_SM;
  // Branching ratios
  Br_gg_SM = Gh_gg_SM/Gh_SM;
  Br_WW_SM = Gh_WW_SM/Gh_SM;
  Br_ZZ_SM = Gh_ZZ_SM/Gh_SM;
  Br_gmgm_SM = Gh_gmgm_SM/Gh_SM;
  Br_tau_SM = Gh_tau_SM/Gh_SM;
  Br_bb_SM = Gh_bb_SM/Gh_SM;
  /* */

  //set 2 from 1411.6222
  //  ms_min = 500.0; ms_max = 5000.0;
  mu_min = 100.0; mu_max = 2000.0;
  tn_min = 1.5; tn_max = 50.5;
  M1_min = 100.0; M1_max = sqF;
  M2_min = 200.0; M2_max = sqF;
  M3_min = 1500.0; M3_max = 4000.0;
  msL_min = 300.0; msL_max = sqF;
  AF_min = 0.1; AF_max = 1.0;

  a1 = alpha/s2w;
  a2 = alpha/(1-s2w);
  F = sqF*sqF;

  srand(20);
  for(i=1;i<=Np_rnd;i++) {
    flag = 1;
    msL = msL_min+frand()*(msL_max-msL_min);
    ms = ms_min+frand()*(ms_max-ms_min);
    mu = mu_min+frand()*(mu_max-mu_min);
    mu = -mu;
    M1 = M1_min+frand()*(M1_max-M1_min);
    M2 = M2_min+frand()*(M2_max-M2_min);
    M3 = M3_min+frand()*(M3_max-M3_min);
    tn = tn_min+frand()*(tn_max-tn_min);
    AF1 = AF_min+frand()*(AF_max-AF_min);
    AF2 = AF_min+frand()*(AF_max-AF_min);
    Amu = AF_min+frand()*(AF_max-AF_min);
    Atau = AF_min+frand()*(AF_max-AF_min);

    AF = sqrt(0.5*(AF1*AF1+AF2*AF2));
    
    
    Mass_m[1][1] = m_mu*m_mu+(mZ*mZ*cos(2.0*beta)*(s2w-0.5))+msL*msL;
    Mass_m[2][2] = m_tau*m_tau+(mZ*mZ*cos(2.0*beta)*(s2w-0.5))+msL*msL;
    Mass_m[3][3] = m_mu*m_mu-mZ*mZ*cos(2.0*beta)*s2w+msL*msL;
    Mass_m[4][4] = m_tau*m_tau-mZ*mZ*cos(2.0*beta)*s2w+msL*msL;
    Mass_m[1][2] = Mass_m[2][1] = 0.0;
    Mass_m[3][4] = Mass_m[3][4] = 0.0;
    Mass_m[1][3] = Mass_m[3][1] = m_mu*Amu-m_mu*mu*tn;
    Mass_m[2][4] = Mass_m[4][2] = m_tau*Atau-m_tau*mu*tn;
    Mass_m[1][4] = Mass_m[4][1] = v*cos(beta)*AF1*sqF;
    Mass_m[2][3] = Mass_m[3][2] = v*cos(beta)*AF2*sqF;
    
    jacobi(Mass_m,4,eigv,vects,&nrot);
    
    for(j=1;j<=4;j++) {
      if (eigv[j]<=0) flag=0;
    }
    MZZ = M1*s2w+M2*(1-s2w);
    Mgmgm = M1*(1-s2w)+M2*s2w;
    Omega_W = M2*v/F;
    Omega_Z = MZZ*v/F;
  
    beta = atan(tn);
    X = 2.0*mu*mu*mu*sin(2.0*beta)*v+2*M_PI*v*v*v*(a1*M1+a2*M2)*cos(2.0*beta)*cos(2.0*beta);
    tn2 = 2.0*X/(F*(ms*ms-mh*mh));
    teta = 0.5*atan(tn2);
  
    // H->WW  
    delta_V = 2.0;
    Gh_WW_m = G1_W*cos(teta)*cos(teta)-12.0*Omega_W*G2_W*sin(teta)*cos(teta)+4.0*Omega_W*Omega_W*G3_W*sin(teta)*sin(teta);
    Gh_WW_m = Gh_WW_m*delta_V*GF*mh*mh*mh;
    Gh_WW_m = Gh_WW_m/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
    // H->ZZ
    delta_V = 1.0;
    Gh_ZZ_m = G1_Z*cos(teta)*cos(teta)-12.0*Omega_Z*G2_Z*sin(teta)*cos(teta)+4.0*Omega_Z*Omega_Z*G3_Z*sin(teta)*sin(teta);
    Gh_ZZ_m = Gh_ZZ_m*delta_V*GF*mh*mh*mh;
    Gh_ZZ_m = Gh_ZZ_m/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
    // H->gg
    Gh_gg_m = GF*a_s(mh)*a_s(mh)*mh*mh*mh/(36.0*sqrt(2.0)*M_PI*M_PI*M_PI);
    Gh_gg_m = Gh_gg_m*pow(cabs(cos(teta)*(A_Q(tau_t)+A_Q(tau_b)+A_Q(tau_c))+sin(teta)*6.0*M3*M_PI*v/(a_s(mh)*F)),2.0);
    Gh_gg_m = Gh_gg_m*(1.0+(95.0/4.0-7.0/6.0*(double)nf)*a_H);
    // H->tau tau
    Gh_tau_m = Gh_tau_SM*cos(teta)*cos(teta);
    // H->bb
    Gh_bb_m = Gh_bb_SM*cos(teta)*cos(teta);
    // H->gm gm
    Gh_gmgm_m = pow(cabs(cos(teta)*(4.0/3.0*A_f(tau_t)+4.0/3.0*A_f(tau_c)+1.0/3.0*A_f(tau_b)+A_W(tau_W))+sin(teta)*4.0*Mgmgm*v*M_PI/(alpha*F)),2.0)*(1-a_H);
    Gh_gmgm_m = Gh_gmgm_m*GF*mh*mh*mh*alpha*alpha/(128.0*sqrt(2.0)*M_PI*M_PI*M_PI);
    // total width
    Gh_m = Gh_WW_m+Gh_ZZ_m+Gh_gmgm_m+Gh_tau_m+Gh_bb_m+Gh_gg_m;
    // Branching ratios
    Br_gg_m = Gh_gg_m/Gh_m;
    Br_WW_m = Gh_WW_m/Gh_m;
    Br_ZZ_m = Gh_ZZ_m/Gh_m;
    Br_gmgm_m = Gh_gmgm_m/Gh_m;
    Br_tau_m = Gh_tau_m/Gh_m;
    Br_bb_m = Gh_bb_m/Gh_m;
    
    mu_bb = ((Gh_WW_m+Gh_ZZ_m)*Br_bb_m)/((Gh_WW_SM+Gh_ZZ_SM)*Br_bb_SM);
    mu_WW = (Gh_gg_m*Br_WW_m)/(Gh_gg_SM*Br_WW_SM);
    mu_ZZ = (Gh_gg_m*Br_ZZ_m)/(Gh_gg_SM*Br_ZZ_SM);
    mu_gmgm = (Gh_gg_m*Br_gmgm_m)/(Gh_gg_SM*Br_gmgm_SM);
    mu_tau = ((Gh_WW_m+Gh_ZZ_m+Gh_gg_m)*Br_tau_m)/((Gh_WW_SM+Gh_ZZ_SM+Gh_gg_SM)*Br_tau_SM);
    
    //auxillary integrals over phase space for sgoldstinos
    xi_W = mw/ms;
    xi_Z = mZ/ms;
    eta_W = Gw/ms;
    eta_Z = GZ/ms;
    
    G1_W_s = GG1_d(xi_W);
    G1_Z_s = GG1_d(xi_Z);
    G2_W_s = GG2_d(xi_W);
    G2_Z_s = GG2_d(xi_Z);
    G3_W_s = GG3_d(xi_W);
    G3_Z_s = GG3_d(xi_Z); 
    
    //for low masses
    //Np = 1000;
    /*G1_W_s = GG1(xi_W,eta_W,Np);
    G1_Z_s = GG1(xi_Z,eta_Z,Np);
    G2_W_s = GG2(xi_W,eta_W,Np);
    G2_Z_s = GG2(xi_Z,eta_Z,Np);
    G3_W_s = GG3(xi_W,eta_W,Np);
    G3_Z_s = GG3(xi_Z,eta_Z,Np);  */

    tau_t_s = 4.0*(m_t*m_t)/(ms*ms);
    tau_b_s = 4.0*(m_b*m_b)/(ms*ms);
    tau_c_s = 4.0*(m_c*m_c)/(ms*ms);
    tau_W_s = 4.0*(mw*mw)/(ms*ms);
  
    // s->WW  
    delta_V = 2.0;
    Gsg_WW = G1_W_s*sin(teta)*sin(teta)+12.0*Omega_W*G2_W_s*sin(teta)*cos(teta)+4.0*Omega_W*Omega_W*G3_W_s*cos(teta)*cos(teta);
    Gsg_WW = Gsg_WW*delta_V*GF*ms*ms*ms;
    Gsg_WW = Gsg_WW/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
    // s->ZZ
    delta_V = 1.0;
    Gsg_ZZ = G1_Z_s*sin(teta)*sin(teta)+12.0*Omega_Z*G2_Z_s*sin(teta)*cos(teta)+4.0*Omega_Z*Omega_Z*G3_Z_s*cos(teta)*cos(teta);
    Gsg_ZZ = Gsg_ZZ*delta_V*GF*ms*ms*ms;
    Gsg_ZZ = Gsg_ZZ/(16.0*M_PI*M_PI*M_PI*sqrt(2.0));
    // s->gg
    a_H_s = a_s(ms)/M_PI;
    Gsg_gg = GF*a_s(ms)*a_s(ms)*ms*ms*ms/(36.0*sqrt(2.0)*M_PI*M_PI*M_PI);
    Gsg_gg = Gsg_gg*pow(cabs(sin(teta)*(A_Q(tau_t_s)+A_Q(tau_b_s)+A_Q(tau_c_s))-cos(teta)*6.0*M3*M_PI*v/(a_s(ms)*F)),2.0);
    Gsg_gg = Gsg_gg*(1.0+(95.0/4.0-7.0/6.0*(double)nf_s)*a_H_s);
    // s->tau tau
    b_tau_s = sqrt(1-4.0*m_tau*m_tau/(ms*ms));
    Gsg_tau = GF*ms*m_tau*m_tau*b_tau_s*b_tau_s*b_tau_s/(4.0*sqrt(2.0)*M_PI);
    Gsg_tau = Gsg_tau*sin(teta)*sin(teta);
    // s->bb
    D_QCD_s= 1.0+5.67*a_H_s+(35.94-1.36*(double)nf_s)*a_H_s*a_H_s;
    D_QCD_s = D_QCD_s + (164.14-25.77*(double)nf_s+0.259*(double)nf_s*(double)nf_s)*a_H*a_H*a_H;
    D_t_s = a_H_s*a_H_s*(1.57-4.0*log(ms/m_t)/3.0+2.0/9.0*log(run_mQ_s(m_b,ms)/ms));
    Gsg_bb = 3.0*GF*ms*run_mQ_s(m_b,ms)*run_mQ_s(m_b,ms)/(4.0*sqrt(2.0)*M_PI)*(D_QCD_s+D_t_s);
    Gsg_bb = Gsg_bb*sin(teta)*sin(teta);
    // s->gm gm
    Gsg_gmgm = pow(cabs(sin(teta)*(4.0/3.0*A_f(tau_t_s)+4.0/3.0*A_f(tau_c_s)+1.0/3.0*A_f(tau_b_s)+A_W(tau_W_s))-cos(teta)*4.0*Mgmgm*v*M_PI/(alpha*F)),2.0)*(1-a_H_s);
    Gsg_gmgm = Gsg_gmgm*GF*ms*ms*ms*alpha*alpha/(128.0*sqrt(2.0)*M_PI*M_PI*M_PI);
    // s->GG
    Gsg_GG = pow(ms,5.0)/(32.0*M_PI*F*F);
    // total width
    Gsg = Gsg_WW+Gsg_ZZ+Gsg_gmgm+Gsg_tau+Gsg_bb+Gsg_gg+Gsg_GG;
    /*
      printf("ms=%e\n",ms);
      printf("s(teta)=%e\n",sin(teta));
      printf("G(s)=%e\n",Gsg);
      printf("\n");
      printf("G(H->WW)=%e\n",Gsg_WW);
      printf("Br(H->WW)=%e\n",Gsg_WW/Gsg);
      printf("\n");
      printf("G(H->ZZ)=%e\n",Gsg_ZZ);
      printf("Br(H->ZZ)=%e\n",Gsg_ZZ/Gsg);
      printf("\n");
      printf("G(H->gg)=%e\n",Gsg_gg);
      printf("Br(H->gg)=%e\n",Gsg_gg/Gsg);
      printf("\n");
      printf("G(H->gmgm)=%e\n",Gsg_gmgm);
      printf("Br(H->gmgm)=%e\n",Gsg_gmgm/Gsg);
      printf("\n");
      printf("G(H->tau tau)=%e\n",Gsg_tau);
      printf("Br(H->tau tau)=%e\n",Gsg_tau/Gsg);
      printf("\n");
      printf("G(H->bb)=%e\n",Gsg_bb);
      printf("Br(H->bb)=%e\n",Gsg_bb/Gsg);
    */

    /*sigma_sg = sigma_production x Br(sg --> 2gm)*/
    /*Int_br -- Experimental bound on sigma_sg*/
    splint(ms1,Int_ms,yd2_ms,Nms+1,ms,&sigma_sg);
    splint(ms_l,Br_l,yd2_l,N_lim,ms/1000.0,&Int_br);
    sigma_sg = 0.3894e+09*M_PI*M_PI*Gsg_gg*sigma_sg/(8.0*ms*s);
    sigma_sg = sigma_sg*Gsg_gmgm/Gsg;    
    Y = v*cos(beta)*fabs(sin(teta))*AF/sqF;
    //    fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %lf\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,teta,mu_ZZ);      
    //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %lf\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,teta,mu_tau);
    if (flag == 1) {
      if (BR(M1,M2,msL,beta,AF,ms,mu,sqF)<BR_max) {
	if ((mu_bb>=mu_bb_min) && (mu_bb<=mu_bb_max)) {
	  if ( ((mu_WW>=mu_WW_min_A) && (mu_WW<=mu_WW_max_A)) || ((mu_WW>=mu_WW_min_C) && (mu_WW<=mu_WW_max_C))) {
	    if ( ((mu_ZZ>=mu_ZZ_min_A) && (mu_ZZ<=mu_ZZ_max_A)) || ((mu_ZZ>=mu_ZZ_min_C) && (mu_ZZ<=mu_ZZ_max_C))) {
	      if ( ((mu_gmgm>=mu_gmgm_min_A) && (mu_gmgm<=mu_gmgm_max_A)) || ((mu_gmgm>=mu_gmgm_min_C) && (mu_gmgm<=mu_gmgm_max_C))) {
		if ( ((mu_tau>=mu_tau_min_A) && (mu_tau<=mu_tau_max_A)) || ((mu_tau>=mu_tau_min_C) && (mu_tau<=mu_tau_max_C))) {
		  if (sigma_sg<=Int_br) {
		    fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),1,flag);
		  }
		  else {
		    //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),2,flag);
		  }
		}
		else {
		  //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),3,flag);
		}
	      }
	      else {
		//fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),4,flag); 
	      }
	    }
	    else {
	      //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),5,flag);  	
	    }
	  }
	  else {
	    //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),6,flag);  	
	  }
	}
	else {
	  //fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),7,flag);  	
	}
      }
      else {
	//fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf %d %d\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,fabs(sin(teta)),8,flag);  
      }
    }
    if (flag == 1) {
      if (BR(M1,M2,msL,beta,AF,ms,mu,sqF)<BR_max) {
	if ((mu_bb>=mu_bb_min) && (mu_bb<=mu_bb_max)) {
	  if ( ((mu_WW>=mu_WW_min_A) && (mu_WW<=mu_WW_max_A)) || ((mu_WW>=mu_WW_min_C) && (mu_WW<=mu_WW_max_C))) {
	    if ( ((mu_ZZ>=mu_ZZ_min_A) && (mu_ZZ<=mu_ZZ_max_A)) || ((mu_ZZ>=mu_ZZ_min_C) && (mu_ZZ<=mu_ZZ_max_C))) {
	      if ( ((mu_gmgm>=mu_gmgm_min_A) && (mu_gmgm<=mu_gmgm_max_A)) || ((mu_gmgm>=mu_gmgm_min_C) && (mu_gmgm<=mu_gmgm_max_C))) {
		if ( ((mu_tau>=mu_tau_min_A) && (mu_tau<=mu_tau_max_A)) || ((mu_tau>=mu_tau_min_C) && (mu_tau<=mu_tau_max_C))) {
		  if (sigma_sg<=Int_br) {
		    fprintf(output_mu,"%lf \t %e \t %e \t %e \t%e \t%e \n",ms,mu_WW,mu_ZZ,mu_tau,mu_bb,mu_gmgm);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  fclose(output);
  fclose(output_mu);
  return 0;
}

double BR(double M1,double M2,double msL,double beta,double AF,double ms,double mu,double sqF) {
  double b,a22;
  a22 = (5.0/3.0)*alpha/(1-s2w);
  b = 0.5*alpha*a22*a22*pow(m_tau,3.0)*sqF*AF/G_tau;
  b = b*pow(v*M1*cos(beta)/(48.0*M_PI*pow(msL,4.0)),2.0);
  return b;
}

double lambda(double x,double y) {
  double la;
  la = (1-x*x-y*y)*(1-x*x-y*y);
  la = la - 4.0*x*x*y*y;
  return la;
}

double F1(double x,double y,double xi,double eta) {
  double out;
  out = sqrt(lambda(x,y))*(lambda(x,y)+12.0*x*x*y*y);
  out = y*out/((y*y-xi*xi)*(y*y-xi*xi)+xi*xi*eta*eta);
  return out;
}

double F2(double x,double xi,double eta) {
  double out;
  out = 4.0*x*xi*xi*eta*eta/((x*x-xi*xi)*(x*x-xi*xi)+xi*xi*eta*eta);
  return out;
}

double F3(double x,double y,double xi,double eta) {
  double out;
  out = sqrt(lambda(x,y))*x*x*y*y*y*(1-x*x-y*y);
  out = out/((y*y-xi*xi)*(y*y-xi*xi)+xi*xi*eta*eta);
  return out;
}

double F4(double x,double y,double xi,double eta) {
  double out;
  out = x*x*y*y*y*sqrt(lambda(x,y))*(0.5*(x*x+y*y-1)*(x*x+y*y-1)+x*x*y*y);
  out = out/((y*y-xi*xi)*(y*y-xi*xi)+xi*xi*eta*eta);
  return out;
}

double GG1(double xi,double eta,int N) {
  double out;
  int i,j;
  double x,y,dx,dy;
  double sum;
  dx = 1.0/(double)N;
  out = 0.0;
  for(i=0;i<=N-1;i++) {
    x = (i+0.5)*dx;
    sum = 0.0;
    dy = (1.0-x)/(double)N;
    for(j=0;j<=N-1;j++) {
      y = (j+0.5)*dy;
      sum = sum + F1(x,y,xi,eta)*dy;
    }
    out = out + sum*F2(x,xi,eta)*dx;
  }
  return out;
}

double GG2(double xi,double eta,int N) {
  double out;
  int i,j;
  double x,y,dx,dy;
  double sum;
  dx = 1.0/(double)N;
  out = 0.0;
  for(i=0;i<=N-1;i++) {
    x = (i+0.5)*dx;
    sum = 0.0;
    dy = (1.0-x)/(double)N;
    for(j=0;j<=N-1;j++) {
      y = (j+0.5)*dy;
      sum = sum + F3(x,y,xi,eta)*dy;
    }
    out = out + sum*F2(x,xi,eta)*dx;
  }
  out = out/(xi*xi);
  return out;
}

double GG3(double xi,double eta,int N) {
  double out;
  int i,j;
  double x,y,dx,dy;
  double sum;
  dx = 1.0/(double)N;
  out = 0.0;
  for(i=0;i<=N-1;i++) {
    x = (i+0.5)*dx;
    sum = 0.0;
    dy = (1.0-x)/(double)N;
    for(j=0;j<=N-1;j++) {
      y = (j+0.5)*dy;
      sum = sum + F4(x,y,xi,eta)*dy;
    }
    out = out + sum*F2(x,xi,eta)*dx;
  }
  out = out/(xi*xi*xi*xi);
  return out;
}

double GG1_d(double xi) {
  double out;
  out = M_PI*M_PI*(1.0-4.0*xi*xi+12.0*xi*xi*xi*xi);
  out = out*sqrt(1-4.0*xi*xi);
  return out;
}

double GG2_d(double xi) {
  double out;
  out = M_PI*M_PI*xi*xi*(1-2.0*xi*xi);
  out = out*sqrt(1-4.0*xi*xi);
  return out;
}

double GG3_d(double xi) {
  double out;
  out = 0.5*M_PI*M_PI*(1.0-4.0*xi*xi+6.0*xi*xi*xi*xi);
  out = out*sqrt(1-4.0*xi*xi);
  return out;
}

complex f(double tau) {
  complex out;
  if (tau>=1.0) {
    out = asin(1/sqrt(tau))*asin(1/sqrt(tau));
  }
  else {
    out = -0.25*(log((1.0+sqrt(1.0-tau))/(1.0-sqrt(1.0-tau)))-I*M_PI);
    out = out*(log((1.0+sqrt(1.0-tau))/(1.0-sqrt(1.0-tau)))-I*M_PI);
  }
  return out;
}

complex A_Q(double tau) {
  complex out;
  out = 1.5*tau*(1.0+(1.0-tau)*f(tau));
  return out;
}

complex A_W(double tau) {
  complex out;
  out = -(2.0+3.0*tau+3.0*tau*(2.0-tau)*f(tau));
  return out;
}

complex A_f(double tau) {
  complex out;
  out = 2.0*tau*(1.0+(1.0-tau)*f(tau));
  return out;
}

//running constant
double a_s(double Q) {
  int Nf = 5;
  double b0,b1,b2;
  double t;
  double Lambd = 0.214; //GeV
  double out;
  
  b0 = (33.0-2.0*(double)Nf)/(12.0*M_PI);
  b1 = (153.0-19.0*(double)Nf)/(24.0*M_PI*M_PI);
  b2 = (2857.0-5033.0/9.0*Nf+325.0/27.0*Nf*Nf)/(128.0*M_PI*M_PI*M_PI);
  t = log(Q*Q/(Lambd*Lambd));
  
  out = 1.0-b1*log(t)/(b0*b0*t)+(b1*b1*(log(t)*log(t)-log(t)-1.0)+b0*b2)/(b0*b0*b0*b0*t*t);
  out = out/(b0*t);
  return out;
}

double c(double x) {
  double out;
  out = 1.0+1.175*x+1.501*x*x+0.1725*x*x*x;
  out = out*pow(23.0*x/6.0,12.0/23.0);
  return out;
}

double c_s(double x) {
  double out;
  out = 1.0+1.398*x+1.793*x*x-0.6834*x*x*x;
  out = out*pow(7.0*x/2.0,4.0/7.0);
  return out;
}

double mQ_MS_bar(double mq) {
  double out;
  double a_b = a_s(mq)/M_PI;
  out = mq/(1.0+4.0*a_b/3.0+12.4*a_b*a_b);
  return out;
}

double run_mQ(double mq,double Q) {
  double out;
  out = mQ_MS_bar(mq)*c(a_s(Q)/M_PI)/c(a_s(mq)/M_PI);
  return out;
}

double run_mQ_s(double mq,double Q) {
  double out;
  out = mQ_MS_bar(mq)*c_s(a_s(Q)/M_PI)/c_s(a_s(mq)/M_PI);
  return out;
}

double integral(double ms,double sqrt_s,int Iset,int Np) { 
  int iprtn = 0; //gluon
  int i;
  double xmin,xmax,sum;
  double s,x,x1,dx;

  setctq6_(&Iset);
  
  s = sqrt_s*sqrt_s;
  xmax = 1.0;
  sum = 0.0;
  xmin = ms*ms/s;
  dx = (xmax-xmin)/(double)Np;
  for (i=0;i<=Np-1;i++) {
    x = xmin+(i+0.5)*dx;
    x1 = ms*ms/(x*s);
    if (x>1.0) x = 1.0;
    if (x1>1.0) x1 = 1.0;
    sum = sum + ctq6pdf_(&iprtn,&x,&ms)*ctq6pdf_(&iprtn,&x1,&ms)*dx/x;
  }  
  return sum;
}
