

//#include "Pythia8/Pythia.h"

#include<fstream>
#include<iostream>
#include<iomanip>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<cstdlib>

//#include"TH1.h"
//#include"TH2.h"
//#include"TFile.h"

//#include"string.h"

using namespace std;
//using namespace Pythia8;

//#include "fastjet/ClusterSequence.hh"
#include <math.h>
#include <vector>

//using namespace fastjet;

//root for interactive graphics.
//#include "TVirtualPad.h"
//#include "TApplication.h"
//#include "TSystem.h"

//...constant


// for heavy quark radiation table
const int HQener_gn=400;
const int t_gn=75;
const int temp_gn=100;

const double HQener_max=200.0;
const double t_max=15.0;
const double temp_max=0.55;
const double temp_min=0.15;

// for MC initialization of jet partons
const double ipTmin=0.0;
const double ipTmax=70.0;
const double eta_cut=0.5;
const int maxMC=2000000;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate

const int N_p1=100;
const int N_T=50;
const int N_e2=75;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate










//  const double  pi=3.1415926;  //conflict with fastjet  ???
const double    CA=3.0; 
const double    sctr=0.197;	     // 1 GeV^-1 = 0.197 fm
const double    pi=3.1415926; 
//...constant
//const double    pi=3.1415926;
//const double    CA=3.0;
const double    CF=4.0/3.0;
//const double    sctr=0.1973;	     //fm to GeV^-1
const double    hydro_Tc=0.165;
const double    KPamp=2.0;
const double    KPsig=5.0;
const double    KTamp=2.0;
const double    KTsig=0.05*hydro_Tc;
const double    preKT=0.85;













class LBTclass{


public:









//...input with current machine time
//...random number seed (any negative integer)
	  
//   long  NUM1 = -33;
                  
long  NUM1;

//    scattering rate
double Rg[50][100];         //total gluon scattering rate as functions of initial energy and temperature 
double Rg1[50][100];        //gg-gg              CT1
double Rg2[50][100];        //gg-qqbar           CT2
double Rg3[50][100];        //gq-qg              CT3
double Rq[50][100];         //total gluon scattering rate as functions of initial energy and temperature
double Rq3[50][100];        //qg-qg              CT13
double Rq4[50][100];        //qiqj-qiqj          CT4
double Rq5[50][100];        //qiqi-qiqi          CT5
double Rq6[50][100];        //qiqibar-qjqjbar    CT6
double Rq7[50][100];        //qiqibar-qiqibar    CT7
double Rq8[50][100];        //qqbar-gg           CT8

double qhatLQ[50][100];
double qhatG[50][100];
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



double RHQ[50][100];        //total scattering rate for heavy quark
double RHQ11[50][100];      //Qq->Qq
double RHQ12[50][100];      //Qg->Qg
double qhatHQ[50][100];     //qhat of heavy quark
double qhat_over_T3;        //qhat/T^3 for heavy quark as fnc of (T,p)
double D2piT;

// for heavy quark radiation table

double dNg_over_dt_c[75+2][100+1][400+1];
double dNg_over_dt_q[75+2][100+1][400+1];
double dNg_over_dt_g[75+2][100+1][400+1];
double max_dNgfnc_c[75+2][100+1][400+1];
double max_dNgfnc_q[75+2][100+1][400+1];
double max_dNgfnc_g[75+2][100+1][400+1];


double delta_tg;
double delta_temp;
double delta_HQener;

// for MC initialization of jet partons


double initMCX[2000000],initMCY[2000000];




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate













//....Hydro profile
double temphydro[40][80][80][80];
double fractionhydro[40][80][80][80];
double VXhydro[40][80][80][80];
double VYhydro[40][80][80][80];
double VZhydro[40][80][80][80];

//...input information

//...............................................................2019
int Singlestepswitch;
int Force22;
int Force23;
int Force2n;
double radlength;

int ExchangeIDswitch;
//...............................................................2019


int switchsingle;
int switchmedium;

double temp0medium;
double VXmedium;
double VYmedium;
double VZmedium;

double mass;
double px0;
double py0;
double pz0;
double en0;

double Xinitial;
double Yinitial;
double Zinitial;


//...............................................................805


int Reachtauend;



//...............................................................805



int    ncall;
int    nprint;		

int    nrank;                    //NOT NECESSARY

////////////////////////////////////////////////////////////////////////////...2019
int    NRcut;                   //scattering rank
int    NScut;                   //splitting rank

int    NR0cut;                  //scattering rank for negative
int    NS0cut;                  //splitting rank for negative

double dMD;                     //distance between mother and daughter parton
double dMDcut;                  //distance of decoherence

int switchphasecut;

double collenercut;
double collthetacut;

double radenercut;
double radthetacut;
////////////////////////////////////////////////////////////////////////////


double dt;                       //dtime when tauswitch is turned off (0) / dtau when tauswitch is turned on (1)
double ti0;                      //initial time
double timend;                   //end time or tau RENAME
double al;                       //2 dimensional plots range NOT NECESSARY

double ener;                     //initial energy of the jet parton/heavy quark
double amss;                     //initial mass of the jet parton/heavy quark 
double temp;                    
double temp0;                  
double temp00; 
double Ecut;                     //energy cut of the recoiled partons
		
int    Kjet;                     //initial flavor of the jet parton
int	   Kqhat0;                   //Debye screening mass switch RENAME
int    Kalphas;                  //alphas switch		
int    KINT;                     //radiation switch
int    KINT0;                     //radiation switch
int    Kradiation;
int    Kprimary;
int    Kfishbone;
int    tauswitch;                //coordinate switch

int    LBTswitch;                //switch of the Linear Boltzmann Transport

int    switchCoLBT_Hydro;


////////////////////////////////////////////////////////////////////////////...2019
int    switchtwcoll;
////////////////////////////////////////////////////////////////////////////


		
//...time-tau 				
double dtau;		
double tau0;		
double tauend;		
double time0;

//...variables with switch		
double alphas;
double qhat0;                    //Debye mass RENAME
double qhat00;		

int np;                          //number of partons at this time step ti				
		
//...radiation block		
int    icl22;                   
int    icl23;                    //the numerical switch in colljet23
int    iclrad;                   //the numerical switch in radiation 	  
int    isp;                      //the splitting function switch
int    nj;                       //number of leading shower partons 

//...Hydroprofile
double tauhydroend;
double tauhydrofile;
double tauhydro0;
double dtauhydro;
int ntauhydro;
int nxhydro;
int nyhydro;
int netahydro;

//...global variable qhat	
int counth100;
  
//double qhat;	                 //transport parameter

double qhat;	                 //transport parameter

double dng0[101][101];	 //table of dn/dkperp2/dx 	  

double Vtemp[4];

//...time system
double tirad[50000];
double tiscatter[50000];
double tiform[50000];	 //pythia 


double Tint_lrf[50000];          //for heavy quark
double eGluon;
double nGluon;








	  
//....radiated gluon
double radng[50000];	  
double xwm[3];
double wkt2m;
	  
double vf[4];              //flow velocity in tau-eta coordinate
double vfcar[4];           //flow velocity in t-z coordinate

double vp[4];              //position of particle		
double vc0[4];             //flow velocity     

//...dimensions in subroutine colljet and twcoll
double vc[4];              
double pc0[4];
double pc2[4];
double pc3[4];
double pc4[4];		
double p0[4];
double p2[4];
double p3[4];		
double p4[4];

double pc00[4];
double pc30[4];	
		
double pc01[4];
double pb[4];
		
double Pj0[4];
		
double V[4][50000];           //parton position 
double P[4][50000];           //parton 4-momentum
double V0[4][50000];          //negative parton position
double P0[4][50000];          //negative parton 4-momentum

double Prad[4][50000];		
		
int NR[50000];                  //scattering rank

////////////////////////////////////////////////////////////////////////////...2019
int NS[50000];                   //splitting rank

int NR0[50000];                  //scattering rank
int NS0[50000];                  //splitting rank

int Mother[50000];               //ID of initial parton of creation
////////////////////////////////////////////////////////////////////////////
 
int KATT1[50000];               //parton flavor  
int KATT10[50000];              //negative parton flavor

double PGm[4];
double tjp[50000];
double Vfrozen[4][50000];     //parton final 4 coordinate
double Vfrozen0[4][50000];    //negative parton final 4 coordinate 
double Ecmcut;                   //energy cut for free streaming   

//........................................................................................................NCUT

double VV[4][50000];
double VV0[4][50000];
double PP[4][50000];
double PP0[4][50000];
int CAT[50000];
int CAT0[50000];
		
int ncut;
int ncut0;

//........................................................................................................NCUTEND

//////////////
//...test
        int n_coll22;
		
        //int n_coll22=0;		
        int n_coll23;
        int ng_coll23;
        int ng_nrad;
        int n_radiation;
        int ng_radiation;
        int n_gluon;
        int n_sp1;
        int n_sp2;
//      ofstream datEg("./Eg.dat");
//////////////	  


        int ntest22;
        int ntestrad;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate

// Variables for HQ 2->2
double min_p1,min_T,min_e2,max_p1,max_T,max_e2,bin_p1,bin_T,bin_e2;

double distFncB[50][100][75],distFncF[50][100][75],distMaxB[50][100][75],distMaxF[50][100][75];
double distFncBM[50][100],distFncFM[50][100];

//int loopN=10000;
int loopN;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double WT[50000];
double WT0[50000];
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double vcfrozen[4][50000];
double vcfrozen0[4][50000];


double Tfrozen[50000];
double Tfrozen0[50000];


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tempHydro[50000];
double vxHydro[50000];
double vyHydro[50000];
double vzHydro[50000];
double fractionHydro[50000];


























	  
//...functions	  

void LinearBoltzmannTransport(int &n, double &ti);

void trans(double v[4],double p[4]);
void transback(double v[4],double p[4]);
  
float ran0(long *idum);
   
double alphas0(int &Kalphas,double temp0);
double DebyeMass2(int &Kqhat0,double alphas,double temp0);
	  
void lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);	  

void flavor(int &CT,int &KATT0,int &KATT2,int &KATT3,double RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2);
void linear(int KATT,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2,double &RTEg,double &RTEg1,double &RTEg2,double &RTEg3,double &RTEq,double &RTEq3,double &RTEq4,double &RTEq5,double &RTEq6,double &RTEq7,double &RTEq8);
void twflavor(int &CT,int &KATT0,int &KATT2,double E,double T);
void twlinear(int KATT,double E,double T,double &RTEg1,double &RTEg2,double &RTEq6,double &RTEq7,double &RTEq8);
void colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
void twcoll(int CT,double qhat0ud,double v0[4],double p0[4],double p2[4]);
	  
void titau(double ti,double vf[4],double vp[4],double p0[4],double &Vx,double &Vy,double &Veta,double &Xtau);


//void colljet23( double temp, double qhat0ud, double v0[4], double p0[4],double p2[4], double p3[4], double p4[4], double qt, int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);
//void radiation(double qhat0ud,double v0[4],double P2[4],double P3[4],double P4[4],double Pj0[4], int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);

void colljet23( double temp, double qhat0ud, double v0[4], double p0[4],double p2[4], double p3[4], double p4[4], double qt, int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);
void radiation(double qhat0ud,double v0[4],double P2[4],double P3[4],double P4[4],double Pj0[4], int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint);


void rotate(double px,double py,double pz,double p[4],int icc);
void dngintt(double dtint,double tint,double tisint,double tirint,double dtlastrad,double qhat0int,double Elabint,double Ejpint);
double angluon(double dtint,double tint,double tisint,double tirint,double qhat0int,double Elabint,double Ejpint,double dtlrf,double dtlastrad);
double dngrad(double xw0, double wkt20, double Elab0, double Ejp0, double tint0, double tis0, double tir0, double qhat0int0, double dtlastrad);
double dng(double x,double y);
	  
int KPoisson(double alambda);	  

void bulklinear(double tau, double x,double y,double eta,double &ed, double &temp,double &fraction,double &VX,double &VY,double &VZ);






//////////////////////////////////////////////////////////////////////////////////////////////////...NEW



void radiationHQ(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt);  
void collHQ22(int CT,double temp,double qhat,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt);
double Mqc2qc(double s, double t, double M);
double Mgc2gc(double s, double t, double M);

void collHQ23(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int &ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt);
double dNg_over_dxdydt(int parID, double x0g, double y0g, double HQenergy, double HQmass, double temp_med, double Tdiff);
double tau_f(double x0g, double y0g, double HQenergy, double HQmass);
double splittingP(int parID, double z0g);
double lambdas(double kTFnc);
double nflavor(double kTFnc);
double alphasHQ(double kTFnc, double tempFnc);
double nHQgluon(int parID,double dtLRF,double &time_gluon,double &temp_med,double &HQenergy,double &max_Ng);

void read_xyMC(int &numXY);
void jetInitialize(int numXY);

//////////////////////////////////////////////////////////////////////////////////////////////////...NEW



void LBTinitialize();

void LBTclear();





};











