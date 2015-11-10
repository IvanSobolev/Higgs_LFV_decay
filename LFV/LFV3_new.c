#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
extern double ran2_my(int *idum);
extern int srandinter(int seed);
extern float randinter(float a, float b);

double integral(double,double,int,int);

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
double loop_f(double,double); //loop factor, see 1209.1397 for details
double f2n(double); //loop factor from neutralinos; see 1304.2783 for details
double f2lf(double,int);
double f2lg(double,int);
double f2lh(double,int);

double alpha = 1/127.94; //EM coupling constant on M_Z scale
double m_tau = 1.777; //GeV, tau-lepton mass
double m_mu = 0.106; //GeV, mu-lepton mass
double m_e = 0.511e-03; //GeV,  electron mass
double G_tau = 2.27e-12; //GeV, tau-lepton width
double s2w = 0.23126; //sin^2(teta_w)
double v = 174.0; //Higgs field V.E.V.
double mh = 125.0; //GeV, Higgs mass
double mw = 80.39; //GeV, W-boson mass
double Gw = 2.09; //GeV, W-boson width
double mZ = 91.19;//GeV, Z-boson mass
double GZ = 2.50;//GeV, Z-boson width
int nf = 5; 
int nf_s = 6;
double GF = 1.166e-05; //Gev^{-2}, Fermi's constant
double BR_max = 4.4e-08; //limit on BR(tau->mu+gamma)
/*double m_t = 177.1; //GeV, top quark pole mass*/
double m_t = 173.2; //GeV, top quark mass
double m_b = 4.87; //GeV, bottom quark pole mass
double m_c = 1.64; //GeV, charm quark pole mass


int main(void) { 
  double a1,a2;
  double G1_W,G2_W,G3_W,G1_Z,G2_Z,G3_Z;
  double G1_W_s,G2_W_s,G3_W_s,G1_Z_s,G2_Z_s,G3_Z_s;
  double delta_V;
  double b_tau,b_mu,b_tau_s,b_mu_s;
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
  double ghtt,gstt,gsmm,ghmm; //Higgs-tau modified coupling (and the same for sgoldstino)
  double Gh_gg_SM,Gh_WW_SM,Gh_tau_SM,Gh_ZZ_SM,Gh_bb_SM,Gh_gmgm_SM,Gh_SM; //GeV, partial Higgs widths in SM
  double Br_gg_SM,Br_WW_SM,Br_ZZ_SM,Br_tau_SM,Br_bb_SM,Br_gmgm_SM; //Branching ratios of Higgs in SM
  double Gh_gg_m,Gh_WW_m,Gh_tau_m,Gh_ZZ_m,Gh_bb_m,Gh_gmgm_m,Gh_GG_m,Gh_tm_m,Gh_mu_m,Gh_m; //GeV, partial Higgs widths in our model
  double Br_gg_m,Br_WW_m,Br_ZZ_m,Br_tau_m,Br_bb_m,Br_gmgm_m,Br_GG_m,Br_tm_m,Br_mu_m; //Branching ratios of Higgs in our model
  double mu_gg,mu_WW,mu_ZZ,mu_tau,mu_bb,mu_gmgm;
  double mu_WW_max_A,mu_WW_min_A,mu_WW_max_C,mu_WW_min_C,mu_bb_max_C,mu_bb_min_C;
  double mu_ZZ_max_A,mu_ZZ_min_A,mu_ZZ_max_C,mu_ZZ_min_C,mu_bb_max_A,mu_bb_min_A;
  double mu_gmgm_max_A,mu_gmgm_min_A,mu_gmgm_max_C,mu_gmgm_min_C;
  double mu_tau_max_A,mu_tau_min_A,mu_tau_max_C,mu_tau_min_C;
  double Gsg_gg,Gsg_WW,Gsg_tau,Gsg_mu,Gsg_ZZ,Gsg_bb,Gsg_gmgm,Gsg_GG,Gsg_tm,Gsg; //GeV, partial sgoldstino widths
  double Br_gg_s,Br_WW_s,Br_ZZ_s,Br_tau_s,Br_mu_s,Br_bb_s,Br_gmgm_s,Br_tm_s,Br_GG_s; //Branching ratios of sgoldstino in our model
  double BRtmg; // Branching of tau -> mu+gm
  char c1;
  // we use GeV for energy units
  double mu,beta,M1,M2,M3,X,tn,msL; //msL -- slepton scale
  double mu_min,tn_min,M1_min,M2_min,M3_min,msL_min;
  double mu_max,tn_max,M1_max,M2_max,M3_max,msL_max;  
  double Y_s,Y,Y_min =0.0008,Y_max = 0.0036; //from hep-ex/1502.07400
  double BR_tm_min =0.078e-02,BR_tm_max = 1.51e-02; //from hep-ex/1502.07400
  double AF_min,AF_max,AF,AFtm,AFmt,Amu,Atau; //sqrt(A_{ab}^2+A_{ba}^2)/sqrt(2F)
  double AFme,AFem,AFet,AFte,Ae;
  double sqF; /*sqrt(F)*/
  double F;
  double ms_min,ms_max;
  double *ms_l,*Br_l,*yd2_l; //limits, ms_l in TeV, Br_l in pb
  double *ms1,*Int_ms,*yd2_ms; //vectors for interpolating function of cross-section
  double *G1_W_ms,*G1_Z_ms,*G2_W_ms,*G2_Z_ms,*G3_W_ms,*G3_Z_ms; //vectors for interpolating of phase spaces
  double *yd2_G1W,*yd2_G1Z,*yd2_G2W,*yd2_G2Z,*yd2_G3W,*yd2_G3Z; //vectors for interpolating of phase spaces
  double yd1,ydN;
  double ms,d_ms; //sgoldstino mass m_s
  double tn2,teta;
  double sqrt_s; //s.o.m. energy
  double s;
  int is_minus; /*1 -- mu is negative, 0 -- mu is positive*/
  int Np,i,j;
  int Np_rnd = 100000000;
  int N_lim,mode;
  int Nms; //number of points for interpolating function of cross-section
  FILE *output_mu,*output,*input_lim,*output_log,*input;
  FILE *output_BR_h,*output_BR_s;
  double Int_br,sigma_sg;
  double *Mass_m[7],*vects[7];
  double eigv[7];
  double min_mass,min_mass2;
  double Ysmm,Ystt,Ysmt,Ystm,Yhmt,Yhtm,Yhtt,Yhmm;
  double cL,cR; //Wilson coefficients
  double dcWg,dctg,kappa; /*kappa -- coefficient from (A9)*/
  double dcL,dcR;
  double dcLn,dcRn; //neutralino contribution
  double loop;
  int nrot,flag_t;
  int flag[9], not_p[9],ntot,nCMS;
  
  /*not_p[i] -- number of points passed through one of these tests
    0 -- branching ratio of invisible mode
    1 -- slepton matrice is positive defined 
    2 -- tau -> mu+gm decay (with higgs, sgoldstino & neutralino contributions)
    3 -- "h-> bb" signal
    4 -- "h-> WW" signal
    5 -- "h-> ZZ" signal
    6 -- "h-> gm gm" signal
    7 -- "h-> tau tau" signal
    8 -- p p -> s -> gm gm [diphoton resonances searches]
  */
  for(i=0;i<=8;i++) {
    not_p[i] = 0;
  }
  ntot = 0; /*number of points passed through all the tests*/
  nCMS = 0; /*number of points passed through all the tests and inside CMS bounds*/
  for(i=0;i<=6;i++) {
    Mass_m[i] = (double *) calloc(7,sizeof(double));
    vects[i] = (double *) calloc(7,sizeof(double));
  }

  /*Signal strength constraints. "A" -- ATLAS,"C" -- CMS*/

  /*1412.2641
  mu (ggF) = 1.02+0.29-0.26
  mu (VBF) = 1.27+0.53-0.45 
  See p.57 formula (16) in citated paper
  /*We use ggF channel for it is more restrictive*/
  mu_WW_max_A = 1.31;
  mu_WW_min_A = 0.76;

  /*1312.1129*/
  /*2l2nu + 0/1-jet correspond mainly to ggH channel (see Table 22)
    This channel is the most constraining one (see Fig.23)
 */
  mu_WW_max_C = 0.92;
  mu_WW_min_C = 0.54;

  /*See 1408.5191 (11.Summary) for details
  mu = 1.44+0.40-0.33 is combined result for ggH+ttH+bbH+VBH+VH channels
  mu (ggH+ttH+bbH) = 1.7+0.5-0.4
  mu (VBF+VH) = 0.3+1.6-0.9 */
  mu_ZZ_max_A = 1.84;
  mu_ZZ_min_A = 1.11;

  /*1312.5353*/
  /* mu = 0.93+0.26-0.23(stat)+0.13-0.09(syst)
     This a combined result
     We use this one or \mu for 0/1-jet category
     (mu = 0.83+0.31-0.25) where the number of events in
     ggH is dominant (see Table 5) 
   */
  mu_ZZ_max_C = 1.14;
  mu_ZZ_min_C = 0.58;
  /*  mu_ZZ_max_C = 1.32;
  mu_ZZ_min_C = 0.61; */

  /*These are combined results from ATLAS group (see PHYSICAL REVIEW D 90, 112015 (2014))*/
  /*Signals from ggF are 1.32+-0.38*/
  /*Combined signal mu = 1.17+-0.27*/
  /*Signals from other production channels are less constraining*/
  /*  mu_gmgm_max_A = 1.44; 
  mu_gmgm_min_A = 0.8; */
  mu_gmgm_max_A = 1.70; 
  mu_gmgm_min_A = 0.94;
  /*\mu = 1.14+0.26-0.23 is combined signal from 7 TeV and 8 TeV datasets taking into account*/
  /*ggH signal \mu = 1.12+0.37-0.32 is less constraining than combined one*/
  /*the rest of signals are even less stringent*/
  /*see 1407.0558 for details*/
  /*  mu_gmgm_max_C = 1.40;
      mu_gmgm_min_C = 0.91;*/
  mu_gmgm_max_C = 1.49;
  mu_gmgm_min_C = 0.80;

  /*1501.04943, see Figure 10 for combined result (which we actually use)
    For each production channel (ggF, VBF+VH) see Table 14
  */
  mu_tau_max_A = 1.8;
  mu_tau_min_A = 1.0;

  /*1401.5041, 1401.6527 -- the same results for combined signal strength*/
  /*Combined = ggH +  VH + VBF*/
  mu_tau_max_C = 1.05;
  mu_tau_min_C = 0.51;

  /*1409.6212*/
  mu_bb_max_A = 0.91;
  mu_bb_min_A = 0.58;
  /*Papers about H->bb by CMS:*/
  /*1506.01010, 1401.6527, 1310.3687*/
  /*Production channel -- VZ*/
  /*Combined results*/
  /*  mu_bb_max_C = 1.40;
      mu_bb_min_C = 0.60;*/
  /*VH results (see Table 5 from 1506.01010)*/
  mu_bb_max_C = 1.32;
  mu_bb_min_C = 0.46;
  /* Reading limits on sgoldstino cross_section X Br(sg --> 2*gm)*/
  input_lim = fopen("sg_cs_limit","r");
  input = fopen("LFV_input","r");
  //input_lim = fopen("sg_cs_limit_low","r");
  
  output = fopen("LFV_output","w");
  output_mu = fopen("LFV_output_mu","w");
  output_log = fopen("LFV_output_log","w");
  output_BR_h = fopen("LFV_output_BR_h","w");
  output_BR_s = fopen("LFV_output_BR_s","w"); 

  //output = fopen("LFV_output_low","w");
  //output_mu = fopen("LFV_output_mu_low","w");

  fscanf(input,"%lf\n",&sqF);
  fscanf(input,"%d\n",&is_minus);
  fclose(input);
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
  /* 
  G1_W_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G1W = (double *) calloc(Nms+1,sizeof(double));
  G1_Z_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G1Z = (double *) calloc(Nms+1,sizeof(double));
  G2_W_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G2W = (double *) calloc(Nms+1,sizeof(double));
  G2_Z_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G2Z = (double *) calloc(Nms+1,sizeof(double));
  G3_W_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G3W = (double *) calloc(Nms+1,sizeof(double));
  G3_Z_ms = (double *) calloc(Nms+1,sizeof(double));
  yd2_G3Z = (double *) calloc(Nms+1,sizeof(double)); 
*/
  mode = 4;
  Np = 1000;
  for(i=0;i<=Nms;i++) {
    ms1[i] = ms_min+d_ms*i;
    Int_ms[i] = integral(ms1[i],sqrt_s,mode,Np);
    /*    xi_W = mw/ms1[i];
    xi_Z = mZ/ms1[i];
    eta_W = Gw/ms1[i];
    eta_Z = GZ/ms1[i];
    G1_W_ms[i] = GG1(xi_W,eta_W,Np);
    G1_Z_ms[i] = GG1(xi_Z,eta_Z,Np);
    G2_W_ms[i] = GG2(xi_W,eta_W,Np);
    G2_Z_ms[i] = GG2(xi_Z,eta_Z,Np);
    G3_W_ms[i] = GG3(xi_W,eta_W,Np);
    G3_Z_ms[i] = GG3(xi_Z,eta_Z,Np);*/
  }
  spline(ms1,Int_ms,Nms+1,yd1,ydN,yd2_ms);  
  /*  spline(ms1,G1_W_ms,Nms+1,yd1,ydN,yd2_G1W);
  spline(ms1,G1_Z_ms,Nms+1,yd1,ydN,yd2_G1Z);
  spline(ms1,G2_W_ms,Nms+1,yd1,ydN,yd2_G2W);
  spline(ms1,G2_Z_ms,Nms+1,yd1,ydN,yd2_G2Z);
  spline(ms1,G3_W_ms,Nms+1,yd1,ydN,yd2_G3W);
  spline(ms1,G3_Z_ms,Nms+1,yd1,ydN,yd2_G3Z); */
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

  /*2-loop Wilson coefficients for tau->mu +gm with Higgs*/
  dcWg = 3.0*f2lf(0.25*tau_W,10000)+5.0*f2lg(0.25*tau_W,10000)+0.75*f2lg(0.25*tau_W,10000)+0.75*f2lh(0.25*tau_W,10000)+0.5*(f2lf(0.25*tau_W,10000)-f2lg(0.25*tau_W,10000))/(0.25*tau_W);
  dctg = f2lf(0.25*tau_t,100000);
  kappa = 0.5*alpha*GF*v/(M_PI*m_tau);

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

  //set 2 from 1411.6222 (partially)
  mu_min = 100.0; mu_max = 2000.0;
  tn_min = 1.5; tn_max = 50.5;
  M1_min = 100.0; M1_max = sqF;
  M2_min = 200.0; M2_max = sqF;
  M3_min = 1500.0; /*M3_max = sqF; M3_max = 4000.0*/
  if (sqF>4000.0) M3_max = 4000.0; else M3_max = sqF;
  msL_min = 300.0; msL_max = sqF;
  AF_min = 0.1; AF_max = 1.0;

  a1 = alpha/s2w;
  a2 = alpha/(1-s2w);
  F = sqF*sqF;
  
  //fprintf(output,"# ms \t\t Y \t\t AF \t\t mu \t\t tn \t\t M1 \t\t M2 \t\t M3 \t\t msL \t\t sin(teta) min_mass\n");
  fprintf(output,"# ms \t\t BR_htm \t AF \t\t mu \t\t tn \t\t M1 \t\t M2 \t\t M3 \t\t msL \t       sin(teta)     BRtmg\n");
  fprintf(output_mu,"# ms \t\t mu_WW \t\t mu_ZZ \t\t mu_tau \t mu_bb \t\t mu_gmgm \n");
  fprintf(output_BR_h,"# ms \t\t Gh_m \t\t Br_gg_m \t Br_WW_m \t Br_ZZ_m \t Br_gmgm_m \t Br_tau_m \t Br_bb_m \t Br_GG_m \n");
  fprintf(output_BR_s,"# ms \t\t Gsg \t\t Br_gg_s \t Br_WW_s \t Br_ZZ_s \t Br_gmgm_s \t Br_tau_s \t Br_bb_s \t Br_GG_s \t sigma_sg \n");
  srand(20);
  for(i=1;i<=Np_rnd;i++) {
    /*    msL = msL_min+frand()*(msL_max-msL_min);
    ms = ms_min+frand()*(ms_max-ms_min);
    mu = mu_min+frand()*(mu_max-mu_min);
    if(is_minus) mu = -mu;
    M1 = M1_min+frand()*(M1_max-M1_min);
    M2 = M2_min+frand()*(M2_max-M2_min);
    M3 = M3_min+frand()*(M3_max-M3_min);
    tn = tn_min+frand()*(tn_max-tn_min);
    AFtm = AF_min+frand()*(AF_max-AF_min);
    AFmt = AF_min+frand()*(AF_max-AF_min);
    Amu = AF_min+frand()*(AF_max-AF_min);
    Atau = AF_min+frand()*(AF_max-AF_min);
    AF = sqrt(0.5*(AFtm*AFtm+AFmt*AFmt)); */

    msL = msL_min+randinter(0.0,1.0)*(msL_max-msL_min);
    ms = ms_min+randinter(0.0,1.0)*(ms_max-ms_min);
    mu = mu_min+randinter(0.0,1.0)*(mu_max-mu_min);
    if(is_minus) mu = -mu;
    M1 = M1_min+randinter(0.0,1.0)*(M1_max-M1_min);
    M2 = M2_min+randinter(0.0,1.0)*(M2_max-M2_min);
    M3 = M3_min+randinter(0.0,1.0)*(M3_max-M3_min);
    tn = tn_min+randinter(0.0,1.0)*(tn_max-tn_min);
    AFtm = AF_min+randinter(0.0,1.0)*(AF_max-AF_min);
    AFmt = AF_min+randinter(0.0,1.0)*(AF_max-AF_min);
    Amu = AF_min+randinter(0.0,1.0)*(AF_max-AF_min);
    Atau = AF_min+randinter(0.0,1.0)*(AF_max-AF_min);
    AF = sqrt(0.5*(AFtm*AFtm+AFmt*AFmt));


    AFem = AFet = AFme = AFte = Ae = 0.0;
    /*Mass_m is slepton mass matrix.
      We took into account only soft terms responsible for LR and RL mixing between tau and mu
      See 1304.2783 for details. */
    Mass_m[1][1] = m_e*m_e+(mZ*mZ*cos(2.0*beta)*(s2w-0.5))+msL*msL;
    Mass_m[2][2] = m_mu*m_mu+(mZ*mZ*cos(2.0*beta)*(s2w-0.5))+msL*msL;
    Mass_m[3][3] = m_tau*m_tau+(mZ*mZ*cos(2.0*beta)*(s2w-0.5))+msL*msL;

    Mass_m[4][4] = m_e*m_e-mZ*mZ*cos(2.0*beta)*s2w+msL*msL;
    Mass_m[5][5] = m_mu*m_mu-mZ*mZ*cos(2.0*beta)*s2w+msL*msL;
    Mass_m[6][6] = m_tau*m_tau-mZ*mZ*cos(2.0*beta)*s2w+msL*msL;

    Mass_m[1][4] = Mass_m[4][1] = m_e*Ae*sqF - m_e*mu*tn;
    Mass_m[2][5] = Mass_m[5][2] = m_mu*Amu*sqF-m_mu*mu*tn;
    Mass_m[3][6] = Mass_m[6][3] = m_tau*Atau*sqF-m_tau*mu*tn;

    Mass_m[1][5] = Mass_m[5][1] = v*cos(beta)*AFem*sqF;
    Mass_m[1][6] = Mass_m[6][1] = v*cos(beta)*AFet*sqF;
    
    Mass_m[2][4] = Mass_m[4][2] = v*cos(beta)*AFme*sqF;
    Mass_m[2][6] = Mass_m[6][2] = v*cos(beta)*AFmt*sqF;

    Mass_m[3][4] = Mass_m[4][3] = v*cos(beta)*AFte*sqF;
    Mass_m[3][5] = Mass_m[5][3] = v*cos(beta)*AFtm*sqF;

    Mass_m[1][2] = Mass_m[1][3] = 0;
    Mass_m[2][1] = Mass_m[3][1] = 0; 
    Mass_m[2][3] = Mass_m[3][2] = 0;

    Mass_m[4][5] = Mass_m[4][6] = Mass_m[5][4] = Mass_m[6][4] = 0;
    Mass_m[5][6] = Mass_m[6][5] = 0; 

    /*Calculating eigenvalues of mass matrix using Jacobi procedure (see Numerical recipes for details) */
    jacobi(Mass_m,6,eigv,vects,&nrot);
    
    min_mass2 = eigv[1];
    for(j=1;j<=6;j++) {
      if (eigv[j]<min_mass2) min_mass2 = eigv[j];
    }
    //min_mass is a square of a lightest slepton particle.
    //If it is less than 325^2 the point must be excluded.
    //1403.5294
    if (min_mass2 > 0) min_mass = sqrt(min_mass2);
    beta = atan(tn);
    X = 2.0*mu*mu*mu*sin(2.0*beta)*v+2*M_PI*v*v*v*(a1*M1+a2*M2)*cos(2.0*beta)*cos(2.0*beta);
    tn2 = 2.0*X/(F*(ms*ms-mh*mh));
    teta = 0.5*atan(tn2);
  
    /*Auxillary loop coefficients. See 1209.1397 for details.*/
    Yhmm = m_mu*cos(teta)/(sqrt(2.0)*v)-v*cos(beta)*Amu*sin(teta)/(sqrt(2.0)*sqF);
    Yhtt = m_tau*cos(teta)/(sqrt(2.0)*v)-v*cos(beta)*Atau*sin(teta)/(sqrt(2.0)*sqF);
    Ysmm = m_mu*sin(teta)/(v*sqrt(2.0))+v*cos(beta)*Amu*cos(teta)/(sqrt(2.0)*sqF);
    Ystt = m_tau*sin(teta)/(v*sqrt(2.0))+v*cos(beta)*Atau*cos(teta)/(sqrt(2.0)*sqF);
    Ysmt = v*cos(beta)*AFmt*cos(teta)/(sqrt(2.0)*sqF);
    Ystm = v*cos(beta)*AFtm*cos(teta)/(sqrt(2.0)*sqF);
    Yhmt = -v*cos(beta)*AFmt*sin(teta)/(sqrt(2.0)*sqF);
    Yhtm = -v*cos(beta)*AFtm*sin(teta)/(sqrt(2.0)*sqF);
    /*1-loop coefficients for Higgs and sgoldstino*/
    cL = Yhtt*Yhtm*loop_f(mh,m_tau)+Yhmm*Yhtm*loop_f(mh,m_mu)+Ystt*Ystm*loop_f(ms,m_tau)*mh*mh/(ms*ms)+Ysmm*Ystm*loop_f(ms,m_mu)*mh*mh/(ms*ms);
    cR = Yhtt*Yhmt*loop_f(mh,m_tau)+Yhmm*Yhmt*loop_f(mh,m_mu)+Ystt*Ysmt*loop_f(ms,m_tau)*mh*mh/(ms*ms)+Ysmm*Ysmt*loop_f(ms,m_mu)*mh*mh/(ms*ms);
    /*2-loop coefficients (Higgs only)*/
    dcL = (-6.0*kappa*(2.0/3.0)*(2.0/3.0)*Yhtm*cos(teta)*dctg+kappa*Yhtm*dcWg)*12*mh*mh;
    dcR = (-6.0*kappa*(2.0/3.0)*(2.0/3.0)*Yhmt*cos(teta)*dctg+kappa*Yhmt*dcWg)*12*mh*mh;
    /*neutralino loops*/
    dcLn = M_PI*(5.0/3.0)*a2*sqF*AFmt*v*cos(beta)*(M1/m_tau)*2.0*f2n(M1*M1/(msL*msL))/pow(msL,4.0)*12*mh*mh;
    dcRn = M_PI*(5.0/3.0)*a2*sqF*AFtm*v*cos(beta)*(M1/m_tau)*2.0*f2n(M1*M1/(msL*msL))/pow(msL,4.0)*12*mh*mh;
    /*Full Wilson coefficients*/
    cL = cL + dcL + dcLn;
    cR = cR + dcR + dcRn;
    /*sgoldstino couplings to gauge bosons*/
    MZZ = M1*s2w+M2*(1-s2w);
    Mgmgm = M1*(1-s2w)+M2*s2w;
    Omega_W = M2*v/F;
    Omega_Z = MZZ*v/F;

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
    b_tau = sqrt(1-4.0*m_tau*m_tau/(mh*mh));
    ghtt = m_tau/(sqrt(2.0)*v)*cos(teta)-Atau*v*cos(beta)*sin(teta)/(sqrt(2.0)*sqF);
    Gh_tau_m = ghtt*ghtt*mh*b_tau*b_tau*b_tau/(8.0*M_PI);
    // H->mu mu
    b_mu = sqrt(1-4.0*m_mu*m_mu/(mh*mh));
    ghmm = m_mu/(sqrt(2.0)*v)*cos(teta)-Amu*v*cos(beta)*sin(teta)/(sqrt(2.0)*sqF);
    Gh_mu_m = ghmm*ghmm*mh*b_mu*b_mu*b_mu/(8.0*M_PI);
    // H->bb
    Gh_bb_m = Gh_bb_SM*cos(teta)*cos(teta);
    // H->gm gm
    Gh_gmgm_m = pow(cabs(cos(teta)*(4.0/3.0*A_f(tau_t)+4.0/3.0*A_f(tau_c)+1.0/3.0*A_f(tau_b)+A_W(tau_W))+sin(teta)*4.0*Mgmgm*v*M_PI/(alpha*F)),2.0)*(1-a_H);
    Gh_gmgm_m = Gh_gmgm_m*GF*mh*mh*mh*alpha*alpha/(128.0*sqrt(2.0)*M_PI*M_PI*M_PI);
    // H-> GG
    Gh_GG_m = pow(ms,4.0)*mh/(32.0*M_PI*F*F)*sin(teta)*sin(teta); 
    // H->tau+mu
    Y = v*cos(beta)*fabs(sin(teta))*AF/sqF;
    Gh_tm_m = mh*Y*Y/(8.0*M_PI);
    // total width
    Gh_m = Gh_WW_m+Gh_ZZ_m+Gh_gmgm_m+Gh_tau_m+Gh_bb_m+Gh_gg_m+Gh_GG_m+Gh_tm_m+Gh_mu_m;
    // Branching ratios
    Br_gg_m = Gh_gg_m/Gh_m;
    Br_WW_m = Gh_WW_m/Gh_m;
    Br_ZZ_m = Gh_ZZ_m/Gh_m;
    Br_gmgm_m = Gh_gmgm_m/Gh_m;
    Br_tau_m = Gh_tau_m/Gh_m;
    Br_bb_m = Gh_bb_m/Gh_m;
    Br_GG_m = Gh_GG_m/Gh_m;
    Br_tm_m = Gh_tm_m/Gh_m;
    Br_mu_m = Gh_mu_m/Gh_m;
    // Signal strengths
    mu_bb = ((Gh_WW_m+Gh_ZZ_m)*Br_bb_m)/((Gh_WW_SM+Gh_ZZ_SM)*Br_bb_SM);
    mu_WW = (Gh_gg_m*Br_WW_m)/(Gh_gg_SM*Br_WW_SM);
    mu_ZZ = ((Gh_WW_m+Gh_ZZ_m+Gh_gg_m)*Br_ZZ_m)/((Gh_WW_SM+Gh_ZZ_SM+Gh_gg_SM)*Br_ZZ_SM);
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
    
    /*    splint(ms1,G1_W_ms,yd2_G1W,Nms+1,ms,&G1_W_s);
    splint(ms1,G1_Z_ms,yd2_G1Z,Nms+1,ms,&G1_Z_s);
    splint(ms1,G2_W_ms,yd2_G2W,Nms+1,ms,&G2_W_s);
    splint(ms1,G2_Z_ms,yd2_G2Z,Nms+1,ms,&G2_Z_s);
    splint(ms1,G3_W_ms,yd2_G3W,Nms+1,ms,&G3_W_s);
    splint(ms1,G3_Z_ms,yd2_G3Z,Nms+1,ms,&G3_Z_s); */
  
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
    gstt = m_tau/(sqrt(2.0)*v)*sin(teta)+Atau*v*cos(beta)*cos(teta)/(sqrt(2.0)*sqF);
    Gsg_tau = gstt*gstt*ms*b_tau_s*b_tau_s*b_tau_s/(8.0*M_PI);
    // s->mu mu
    b_mu_s = sqrt(1-4.0*m_mu*m_mu/(ms*ms));
    gsmm = m_mu/(sqrt(2.0)*v)*sin(teta)+Amu*v*cos(beta)*cos(teta)/(sqrt(2.0)*sqF);
    Gsg_mu = gsmm*gsmm*ms*b_mu_s*b_mu_s*b_mu_s/(8.0*M_PI);
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
    Gsg_GG = pow(ms,5.0)/(32.0*M_PI*F*F)*cos(teta)*cos(teta);
    // s->tau+mu
    Y_s = v*cos(beta)*fabs(cos(teta))*AF/sqF;
    Gsg_tm = mh*Y*Y/(8.0*M_PI);

    // total width
    Gsg = Gsg_WW+Gsg_ZZ+Gsg_gmgm+Gsg_tau+Gsg_mu+Gsg_bb+Gsg_gg+Gsg_GG+Gsg_tm;
    // branching ratios of sgoldstino
    Br_gg_s = Gsg_gg/Gsg;
    Br_WW_s = Gsg_WW/Gsg;
    Br_ZZ_s = Gsg_ZZ/Gsg;
    Br_gmgm_s = Gsg_gmgm/Gsg;
    Br_tau_s = Gsg_tau/Gsg;
    Br_tm_s = Gsg_tm/Gsg;
    Br_mu_s = Gsg_mu/Gsg;
    Br_bb_s = Gsg_bb/Gsg;
    Br_GG_s = Gsg_GG/Gsg;
  
    /*sigma_sg = sigma_production x Br(sg --> 2gm)*/
    /*Int_br -- Experimental bound on sigma_sg*/
    splint(ms1,Int_ms,yd2_ms,Nms+1,ms,&sigma_sg);
    splint(ms_l,Br_l,yd2_l,N_lim,ms/1000.0,&Int_br);
    sigma_sg = 0.3894e+09*M_PI*M_PI*Gsg_gg*sigma_sg/(8.0*ms*s);
    sigma_sg = sigma_sg*Br_gmgm_s;

    BRtmg = alpha*pow(m_tau,5.0)/(9216*pow(mh*M_PI,4.0))*(cL*cL+cR*cR)/G_tau;

    flag[0] = ( Br_GG_m < 0.25 ); /*1509.00672*/
    flag[1] = ( min_mass2 > 325*325 );
    flag[2] = ( BRtmg < BR_max );
    flag[3] = ( ((mu_bb>=mu_bb_min_A) && (mu_bb<=mu_bb_max_A)) || ((mu_bb>=mu_bb_min_C) && (mu_bb<=mu_bb_max_C)) );
    flag[4] = ( ((mu_WW>=mu_WW_min_A) && (mu_WW<=mu_WW_max_A)) || ((mu_WW>=mu_WW_min_C) && (mu_WW<=mu_WW_max_C)));
    flag[5] = ( ((mu_ZZ>=mu_ZZ_min_A) && (mu_ZZ<=mu_ZZ_max_A)) || ((mu_ZZ>=mu_ZZ_min_C) && (mu_ZZ<=mu_ZZ_max_C)));
    flag[6] = ( ((mu_gmgm>=mu_gmgm_min_A) && (mu_gmgm<=mu_gmgm_max_A)) || ((mu_gmgm>=mu_gmgm_min_C) && (mu_gmgm<=mu_gmgm_max_C)));
    flag[7] = ( ((mu_tau>=mu_tau_min_A) && (mu_tau<=mu_tau_max_A)) || ((mu_tau>=mu_tau_min_C) && (mu_tau<=mu_tau_max_C)));
    flag[8] = ( sigma_sg<=Int_br );

    for (j=0;j<=8;j++) {
      not_p[j]+=flag[j];
    }
    flag_t = 1;
    for (j=0;j<=8;j++) {
      flag_t = flag_t*flag[j];
    };
    if (flag_t == 1) {
      ntot++;
      //      if ( (Y>Y_min) && (Y<Y_max) ) nCMS++;
      if ( (Br_tm_m>BR_tm_min) && (Br_tm_m<BR_tm_max) ) nCMS++;
      //      fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf   %lf\n",ms,Y,AF,mu,tn,M1,M2,M3,msL,sin(teta),min_mass);
      fprintf(output,"%lf \t %e \t %lf \t %lf \t %lf \t %lf \t%lf \t%lf \t%lf \t%lf   %e\n",ms,Br_tm_m,AF,mu,tn,M1,M2,M3,msL,sin(teta),BRtmg);
      fprintf(output_mu,"%lf \t %e \t %e \t %e \t%e \t%e \n",ms,mu_WW,mu_ZZ,mu_tau,mu_bb,mu_gmgm);
      fprintf(output_BR_h,"%lf \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",ms,Gh_m,Br_gg_m,Br_WW_m,Br_ZZ_m,Br_gmgm_m,Br_tau_m,Br_bb_m,Br_GG_m);
      fprintf(output_BR_s,"%lf \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n",ms,Gsg,Br_gg_s,Br_WW_s,Br_ZZ_s,Br_gmgm_s,Br_tau_s,Br_bb_s,Br_GG_s,sigma_sg);
    }
  }
  
  fprintf(output_log,"# sqF = %lf\n",sqF);
  if (is_minus) {
    fprintf(output_log,"# mu is negative\n");
  }
  else {
    fprintf(output_log,"# mu is positive\n");
  }
  fprintf(output_log,"# Total number of points used \t %d\n",Np_rnd);
  fprintf(output_log,"# Total percent of points that passed through all the tests \t %e\n",(double)ntot/(double)Np_rnd);
  fprintf(output_log,"# Total number of points that passed through all the tests \t %d\n",ntot);
  fprintf(output_log,"# Number of points that passed through every particular test: \n");
  fprintf(output_log,"# Invisible branching ratio \t %lf\n",(double)not_p[0]/(double)Np_rnd);
  fprintf(output_log,"# Wrong slepton masses \t %lf\n", (double)not_p[1]/(double)Np_rnd);
  fprintf(output_log,"# tau-> mu+gm \t %lf\n",(double)not_p[2]/(double)Np_rnd);
  fprintf(output_log,"# h->bb signal \t %lf\n",(double)not_p[3]/(double)Np_rnd);
  fprintf(output_log,"# h->WW signal \t %lf\n",(double)not_p[4]/(double)Np_rnd);
  fprintf(output_log,"# h->ZZ signal \t %lf\n",(double)not_p[5]/(double)Np_rnd);
  fprintf(output_log,"# h->gm gm signal \t %lf\n",(double)not_p[6]/(double)Np_rnd);
  fprintf(output_log,"# h-> tau tau signal \t %lf\n",(double)not_p[7]/(double)Np_rnd);
  fprintf(output_log,"# Diphoton resonance searches\t %lf\n",(double)not_p[8]/(double)Np_rnd);
  fprintf(output_log,"# Passed all the tests and inside CMS bound \t %e\n",(double)nCMS/(double)Np_rnd);
  fprintf(output_log,"# Passed all the tests and inside CMS bound \t %d\n",nCMS);
  
  fclose(output_log);
  fclose(output);
  fclose(output_mu);
  fclose(output_BR_h);
  fclose(output_BR_s);
  return 0;
}

/*2loop functions from sgoldstino/higgs running in a loop*/
double f2lf(double z, int Np) { /*(A6) from 1209.1397*/
  int i;
  double sum,x,dx;
  dx = 1.0/(double)Np;
  sum = 0.0;
  for (i=0;i<=Np-1;i++) {
    x = (i+0.5)*dx;
    sum = sum + (1.0-2.0*x*(1.0-x))/(x*(1.0-x)-z)*log(x*(1.0-x)/z)*dx;
  }
  sum = 0.5*z*sum;
  return sum;
}

double f2lg(double z, int Np) { /*(A7) from 1209.1397*/
  int i;
  double sum,x,dx;
  dx = 1.0/(double)Np;
  sum = 0.0;
  for (i=0;i<=Np-1;i++) {
    x = (i+0.5)*dx;
    sum = sum + 1.0/(x*(1.0-x)-z)*log(x*(1.0-x)/z)*dx;
  }
  sum = 0.5*z*sum;
  return sum;  
}

double f2lh(double z, int Np) { /*(A8) from 1209.1397*/
  int i;
  double sum,x,dx;
  dx = 1.0/(double)Np;
  sum = 0.0;
  for (i=0;i<=Np-1;i++) {
    x = (i+0.5)*dx;
    sum = sum + 1.0/(z-x*(1.0-x))*(1.0+z/(z-x*(1.0-x))*log(x*(1.0-x)/z))*dx;
  }
  sum = 0.5*z*sum;
  return sum;    
}



/*loop factor from higgs and sgoldstino running in a loop
see 1209.1397 for details*/
double loop_f(double mass1,double mass2) {
  double loop;
  loop = -4.0 + 3.0*log(mass1*mass1/(mass2*mass2));
  return loop;
}

/*loop factor from neutralinos running in a loop
see 1304.2783 for details */
double f2n(double a) {
  double b;
  double a1;
  a1 = a-1.0;
  if (fabs(a-1.0)<=1e-3) {
    b = 1/24.0-a1/30.0+pow(a1,2.0)/40.0-2.0*pow(a1,3.0)/105.0+5.0*pow(a1,4.0)/336.0-pow(a1,5.0)/84.0;
  }
  else {
    b = (-5.0*a*a+4.0*a+1.0+2.0*a*(a+2.0)*log(a))/(4.0*pow(1-a,4));
  }
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
