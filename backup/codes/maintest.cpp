#include"LBT.h"

#include "Pythia8/Pythia.h"

#include"TH1.h"
#include"TH2.h"
#include"TFile.h"


#include "fastjet/ClusterSequence.hh"
#include <math.h>
#include <vector>

//................................................................jet grooming
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "fastjet/tools/CASubJetTagger.hh"

#include "fastjet/Selector.hh"

//double Sjet1px, Sjet1py, Sjet1pz, Sjet1E, Sjet1M; 
//double Sjet2px, Sjet2py, Sjet2pz, Sjet2E, Sjet2M;

#include "SoftDrop.hh" // In external code, this should be fastjet/contrib/SoftDrop.hh
//................................................................jet grooming



using namespace fastjet;

//...................................................................................................... fastjet negative modification.

typedef fastjet::JetDefinition::Recombiner Recombiner;
/// Recombiner class that propagates the user index and arranges the
/// recombination accordingly
class NegativeEnergyRecombiner : public  Recombiner {
public:
  NegativeEnergyRecombiner(const int ui) : _ui(ui) {}

  virtual std::string description() const {return "E-scheme Recombiner that checks a flag for a 'negative momentum' particle, and subtracts the 4-momentum in recombinations.";}

  /// recombine pa and pb and put result into pab
  virtual void recombine(const fastjet::PseudoJet & pa, 
                         const fastjet::PseudoJet & pb, 
                         fastjet::PseudoJet & pab) const {

     int ai=1,bi=1;
     
     // If a particle is flagged, restore its real negative energy. 
     // The -1 will flip the full 4-momentum, reversing the convention for 
     // negative energy particles.
     if (pa.user_index() < 0) { ai = -1;}
     if (pb.user_index() < 0) { bi = -1;}

     // recombine particles
     pab = ai*pa+bi*pb;
     
     // if the combination has negative energy, flip the whole 4-momentum and flag it, 
     // so that we have the same convention as for input particles
     if(pab.E() < 0) { 
	pab.set_user_index(_ui); 
        pab.reset_momentum(-pab.px(),
	                   -pab.py(),
	                   -pab.pz(),
   	                   -pab.E());
//        cout << "XX3 " << pab.E()  << " " << tmppa.E() << " " << tmppb.E() << endl;
     } else { pab.set_user_index(0);}

  }

private:
  const int _ui;
};

//...................................................................................................... fastjet negative modification end.

double Jptsmearing;
double Jptinitial;



int main(int argc, char **argv){	  				

//..define LBTclass
	LBTclass *LBT = new LBTclass();

//..LBT initialization(input(output) setting, datafile read in)	
    LBT->LBTinitialize();		

	
     double R=0.4;
	
	
	 char readingeometry_route[1024];	
	
     //readingeometry_route="/home/luotan/August-gammajet/2.76/Hydro/030/readindatafile/Hydroprofile";	

     strcpy(readingeometry_route, "/home/yayun/August-gammajet/2.76/Hydro/030/readindatafile/Hydroprofile");

	 
	 char readinpp_route[1024];	
	
     //readinpp_route="/home/luotan/August-gammajet/2.76/run/pythia6/80/1";
	
     //strcpy(readinpp_route, "/home/luotan/August-gammajet/2.76/run/pythia6/80/multiple/1");
     //strcpy(readinpp_route, "/home/luotan/August-gammajet/2.76/run/pythia6/80/6");

     strcpy(readinpp_route, "/home/yayun/August-gammajet/2.76/run/pythia6/80/10");	 
	 
     //strcpy(readinpp_route, "//home/luotan/final/testrecombination/test2/pythiadata");
	
	
//..datafile of all partons for recombination
    ofstream positiveLBT("./outputdatafile/positive.dat");
    ofstream negativeLBT("./outputdatafile/negative.dat");

    ofstream numpositiveLBT("./outputdatafile/numptpositive.dat");
    ofstream numnegativeLBT("./outputdatafile/numptnegative.dat");	
	
    ofstream inforLBT("./LBTinfor.dat");		
	
	
    TFile *LBTjet= new TFile("./LBTjet.root","Recreate");	 	
	
	int switchroot=1;
    //if(switchroot>0)
    //{	

//................................................................histogram setting.

//..histogram setting for Fastjet
	 
	 double       deltaphi;
	 double       Jpt;
	 double       Jeta;
	
//     double	  JEtarange=10.0;
//     double	  JPhirange=3.1415926;
//     int	      ndJeta=50;
//     int	      ndJphi=50;


     double	      JEtarange=10.0;
     double	      SJEtarange=3.0;	 
     double	      JPhirange=3.1415926;
     double	      SJPhirange=3.1415926;	 
     int	      ndJeta=121;
     int	      ndJphi=121;


	 
     double	      JPtrange=200;
     double	      JErange=200;
     double       JThetarange=3.1415926;
     int	      ndJpt=20;
     int	      ndJe=100;
     int	      ndJtheta=100;
	
     int	      nJtstep=10;
     double       dJtstep=1.0;
		
     int          ndndxpl=100;
     int          ndndxpl0=100;

     double       Jxrange=1.0;
     double       Jx0range=1.0;

     int          ndJR=20;
	 
     double       ndJRsuper=20;
 	 
     double       JRrange=R;	
	 
     double       JRrangesuper=1.0;  	 
	 
//   double       LBT->timend=10.0;     //@@@???
     int          ntistep=10;

//..histogram setting for LBT
	 
     double 	  Rrange=10.0;
     double 	  Xrange=10.0;
     double 	  Yrange=10.0;
     double 	  Zrange=10.0;	 
     int  	      ndr=80;
     int	      ndz=80;
	
     double	      Etarange=4.5;
     double	      Phirange=3.1415926;
     int	      ndeta=50;
     int  	      ndphi=50;
	
     double	      Ptrange=100;
     double	      Erange=100;
     double	      Thetarange=3.1415926;
     int	      ndpt=100;
     int	      nde=100;
     int	      ndtheta=100;
	
     int	      ntstep=10;
     double       dtstep=1.0;		 
 
	 
//................................................................histogram setting end.	


	 
//............................................................................................Histograms setting

//...Book histogram for fastjet.


//.................................................................................803

    TH1D *h100Rjgamma = new TH1D("h100Rjgamma","", 40, 0.0, 200.0);
    TH1D *h100Ngamma = new TH1D("h100Ngamma","", 40, 0.0, 200.0);
    TH1D *h100Njgamma = new TH1D("h100Njgamma","", 40, 0.0, 200.0);
    TH1D *h100Xjgamma = new TH1D("h100Xjgamma","", 40, 0.0, 200.0);	
    //double ww100 = 1/(2.0/16)/LBT->ncall; 

//.................................................................................803

//...gamma jet asymmetry
    TH1D *h100 = new TH1D("h100","", 16, 0.0, 2.0);
    double ww100 = 1/(2.0/16.0)/LBT->ncall;   		


//.................................................................................jet grooming
    TH1D *hzg = new TH1D("hzg","", 10, 0.0, 0.5);
    double wwzg = 1.0/(0.5/10.0)/LBT->ncall;

    TH1D *hmassjpt = new TH1D("hmassjpt","", 27, 0.0, 0.27);
    double wwmassjpt = 1.0/(0.27/27.0)/LBT->ncall;
	
    double wwRg = 1.0/(0.5/10.0)/LBT->ncall;	

    double wwSDtheta = 1.0/(1.0/100.0)/LBT->ncall;

    double wwgroomedmass = 1.0/(Ptrange/100.0)/LBT->ncall; 

    double wwSDlnztheta = 1.0/(5.0/50.0)/(12.0/50.0)/LBT->ncall;
	
    int numgroomedjet[100] = {0};
	
    ntistep = ntistep + 1;

    TH1D *SDzg[ntistep];
    TH1D *SDmassjpt[ntistep];
    TH1D *SDmass[ntistep];
    TH1D *SDRg[ntistep];
    TH1D *SDtheta[ntistep];
	
    TH2D *SDlnztheta[ntistep];
    TH2D *SDEtheta[ntistep];
	
    ntistep = ntistep - 1;


    string h;
    const char *hname;
	
    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	

    h = "SDzg";
    hname = (h).c_str();
    SDzg[ik] = new TH1D(Form("SDzg%d",ik), "", 10, 0.0, 0.5);	
	
    h = "SDmassjpt";
    hname = (h).c_str();
    SDmassjpt[ik] = new TH1D(Form("SDmassjpt%d",ik), "", 27, 0.0, 0.27);	

    h = "SDmass";
    hname = (h).c_str();
    SDmass[ik] = new TH1D(Form("SDmass%d",ik), "", 100, 0.0, Ptrange);	
	
    h = "SDRg";
    hname = (h).c_str();
    SDRg[ik] = new TH1D(Form("SDRg%d",ik), "", 10, 0.0, 0.5);

    h = "SDtheta";
    hname = (h).c_str();
    SDtheta[ik] = new TH1D(Form("SDtheta%d",ik), "", 100, 0.0, 1.0);

    h = "SDlnztheta";
    hname = (h).c_str();	
    SDlnztheta[ik] = new TH2D(Form("SDlnztheta%d",ik), "", 50, 0.0, 5.0, 50, -10.0, 2.0);

    h = "SDEtheta";    	
    SDEtheta[ik] = new TH2D(Form("SDEtheta%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);
	
    }
//.................................................................................jet grooming



///////////////////////////////////////////////////////////////////////////////////////////////////
	

    TH1D *SleadingJetenergy;
    TH1D *SleadingJetpt;
    TH1D *SleadingJetpesudopt;	
    TH1D *SleadingJetelossf;
    TH1D *SleadingJeteloss;
    TH1D *SleadingJettheta;
    TH1D *SleadingJetphi;
    TH1D *SleadingJeteta;

    //string h;

    //const char *hname;

    h = "SleadingJetenergy";
    hname = (h).c_str();
    SleadingJetenergy = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJetpt";
    hname = (h).c_str();
    SleadingJetpt = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJetpesudopt";
    hname = (h).c_str();
    SleadingJetpesudopt = new TH1D(hname, "", 10, -0.5, 9.5);
	
    h = "SleadingJetelossf";
    hname = (h).c_str();
    SleadingJetelossf = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJeteloss";
    hname = (h).c_str();
    SleadingJeteloss = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJettheta";
    hname = (h).c_str();
    SleadingJettheta = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJetphi";
    hname = (h).c_str();
    SleadingJetphi = new TH1D(hname, "", 10, -0.5, 9.5);

    h = "SleadingJeteta";
    hname = (h).c_str();
    SleadingJeteta = new TH1D(hname, "", 10, -0.5, 9.5);


///////////////////////////////////////////////////////////////////////////////////////////////////////


	
    TH1D *leadingJetenergy;
    TH1D *leadingJetpt;
    TH1D *leadingJetpesudopt;	
    TH1D *leadingJetelossf;
    TH1D *leadingJeteloss;
    TH1D *leadingJettheta;
    TH1D *leadingJetphi;
    TH1D *leadingJeteta;

//    string h;

//    const char *hname;

    h = "leadingJetenergy";
    hname = (h).c_str();
    leadingJetenergy = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJetpt";
    hname = (h).c_str();
    leadingJetpt = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJetpesudopt";
    hname = (h).c_str();
    leadingJetpesudopt = new TH1D(hname, "", ntistep, 0.0, LBT->timend);
	
    h = "leadingJetelossf";
    hname = (h).c_str();
    leadingJetelossf = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJeteloss";
    hname = (h).c_str();
    leadingJeteloss = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJettheta";
    hname = (h).c_str();
    leadingJettheta = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJetphi";
    hname = (h).c_str();
    leadingJetphi = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingJeteta";
    hname = (h).c_str();
    leadingJeteta = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

	ntistep=ntistep + 1;
	
    TH2D *dndedtleadingJet0[ntistep];
    TH2D *dndptdtleadingJet0[ntistep];
    TH2D *dndpt2dtleadingJet0[ntistep];
    TH2D *dndthetadtleadingJet0[ntistep];
    TH2D *dndphidtleadingJet0[ntistep];	
    TH2D *dndetadtleadingJet0[ntistep];

	ntistep=ntistep - 1;	

    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {
	
    h = "dndedtleadingJet0";
    hname = (h).c_str();	
    dndedtleadingJet0[ik] = new TH2D(Form("dndedtleadingJet0%d",ik), "", ndJe, 0.0, JErange, ntistep, 0.0, LBT->timend);

    h = "dndptdtleadingJet0";	
    hname = (h).c_str();
    dndptdtleadingJet0[ik] = new TH2D(Form("dndptdtleadingJet0%d",ik), "", ndJpt, 0.0, JPtrange, ntistep, 0.0, LBT->timend);

    h = "dndpt2dtleadingJet0";	
    hname = (h).c_str();		
    dndpt2dtleadingJet0[ik] = new TH2D(Form("dndpt2dtleadingJet0%d",ik), "", ndJpt, 0.0, JPtrange, ntistep, 0.0, LBT->timend);

    h = "dndthetadtleadingJet0";	
    hname = (h).c_str();	
    dndthetadtleadingJet0[ik] = new TH2D(Form("dndthetadtleadingJet0g%d",ik), "", ndJtheta, 0.0, JThetarange, ntistep, 0.0, LBT->timend);

    h = "dndphidtleadingJet0";	
    hname = (h).c_str();	
    dndphidtleadingJet0[ik] = new TH2D(Form("dndphidtleadingJet0%d",ik), "", ndJphi, 0.0, JPhirange, ntistep, 0.0, LBT->timend);

    h = "dndetadtleadingJet0";	
    hname = (h).c_str();
    dndetadtleadingJet0[ik] = new TH2D(Form("dndetadtleadingJet0%d",ik), "", ndJeta, -JEtarange, JEtarange, ntistep, 0.0, LBT->timend);

	}
	
    ntistep = ntistep + 1;

    TH1D *dndeleadingJet[ntistep];
    TH1D *dndptleadingJet[ntistep];
    TH1D *dndpt2leadingJet[ntistep];
    TH1D *dndpesudoptleadingJet[ntistep];
    TH1D *dndpesudopt2leadingJet[ntistep];		
    TH1D *dndthetaleadingJet[ntistep];
    TH1D *dndphileadingJet[ntistep];
    TH1D *dndetaleadingJet[ntistep];
	
    TH2D *dedetadphileadingJet[ntistep];
    TH2D *dndetadphileadingJet[ntistep];
	
    TH1D *dndxplleadingJet[ntistep];
    TH1D *dndxpl0leadingJet[ntistep];

    TH1D *rholeadingJet[ntistep];


    TH1D *rholeadingJetsuper[ntistep];
    TH1D *rholeadingJetsuperCAT[ntistep];	
    TH1D *rholeadingJetsuper0510[ntistep];
    TH1D *rholeadingJetsuper1020[ntistep];	
    TH1D *rholeadingJetsuper2030[ntistep];
    TH1D *rholeadingJetsuper2040[ntistep];
    TH1D *rholeadingJetsuper3040[ntistep];
    TH1D *rholeadingJetsuper4080[ntistep];
    TH1D *rholeadingJetsuper80[ntistep];

    TH2D *SetaphileadingJet[ntistep];
    TH2D *SetaphileadingJetCAT[ntistep];
	
	TH2D *SPTetaphileadingJet[ntistep];
	TH2D *SPTetaphileadingJetCAT[ntistep];	
    TH2D *SPTetaphileadingJet0510[ntistep];	
    TH2D *SPTetaphileadingJet1020[ntistep];	
    TH2D *SPTetaphileadingJet2030[ntistep];
    TH2D *SPTetaphileadingJet2040[ntistep];
    TH2D *SPTetaphileadingJet3040[ntistep];
    TH2D *SPTetaphileadingJet4080[ntistep];
    TH2D *SPTetaphileadingJet80[ntistep];					
	
	
    ntistep = ntistep - 1;

    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	
    h = "dndeleadingJet";
    hname = (h).c_str();
    dndeleadingJet[ik] = new TH1D(Form("dndeleadingJet%d",ik), "", ndJe, 0.0, JErange);

    h = "dndptleadingJet";
    hname = (h).c_str();	
    dndptleadingJet[ik] = new TH1D(Form("dndptleadingJet%d",ik), "", ndJpt, 0.0, JPtrange);	
	
    h = "dndpt2leadingJet";
    hname = (h).c_str();	
    dndpt2leadingJet[ik] = new TH1D(Form("dndpt2leadingJet%d",ik), "", ndJpt, 0.0, JPtrange);	

    h = "dndpesudoptleadingJet";
    hname = (h).c_str();	
    dndpesudoptleadingJet[ik] = new TH1D(Form("dndpesudoptleadingJet%d",ik), "", ndJpt, 0.0, JPtrange);	
	
    h = "dndpesudopt2leadingJet";
    hname = (h).c_str();	
    dndpesudopt2leadingJet[ik] = new TH1D(Form("dndpesudopt2leadingJet%d",ik), "", ndJpt, 0.0, JPtrange);	
	
	
    h = "dndthetaleadingJet";
    hname = (h).c_str();
    dndthetaleadingJet[ik] = new TH1D(Form("dndthetaleadingJet%d",ik), "", ndJtheta, 0.0, JThetarange); 	 

    h = "dndphileadingJet";
    hname = (h).c_str();	
    dndphileadingJet[ik] = new TH1D(Form("dndphileadingJet%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dndetaleadingJet";
    hname = (h).c_str();	
    dndetaleadingJet[ik] = new TH1D(Form("dndetaleadingJet%d",ik), "", ndJeta, -JEtarange, JEtarange);

    h = "dedetadphileadingJet";
    hname = (h).c_str();	
    dedetadphileadingJet[ik] = new TH2D(Form("dedetadphileadingJet%d",ik), "", ndJphi, 0.0, JPhirange, ndJeta, -JEtarange, JEtarange);

    h = "dndetadphileadingJet";
    hname = (h).c_str();	
    dndetadphileadingJet[ik] = new TH2D(Form("dndetadphileadingJet%d",ik), "", ndJphi, 0.0, JPhirange, ndJeta, -JEtarange, JEtarange);

    h = "dndxplleadingJet";
    hname = (h).c_str();	
    dndxplleadingJet[ik] = new TH1D(Form("dndxplleadingJet%d",ik), "", ndndxpl, 0.0, Jxrange);		  

    h = "dndxpl0leadingJet";
    hname = (h).c_str();	
    dndxpl0leadingJet[ik] = new TH1D(Form("dndxpl0leadingJet%d",ik), "", ndndxpl0, 0.0, Jx0range); 

    h = "rholeadingJet";
    hname = (h).c_str();	
    rholeadingJet[ik] = new TH1D(Form("rholeadingJet%d",ik), "", ndJR, 0.0, JRrange);

	h = "rholeadingJetsuper";
    hname = (h).c_str();	
    rholeadingJetsuper[ik] = new TH1D(Form("rholeadingJetsuper%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuperCAT";
    hname = (h).c_str();	
    rholeadingJetsuperCAT[ik] = new TH1D(Form("rholeadingJetsuperCAT%d",ik), "", ndJRsuper, 0.0, JRrangesuper);	
	
	
	h = "rholeadingJetsuper0510";
    hname = (h).c_str();	
    rholeadingJetsuper0510[ik] = new TH1D(Form("rholeadingJetsuper0510%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuper1020";
    hname = (h).c_str();	
    rholeadingJetsuper1020[ik] = new TH1D(Form("rholeadingJetsuper1020%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuper2030";
    hname = (h).c_str();	
    rholeadingJetsuper2030[ik] = new TH1D(Form("rholeadingJetsuper2030%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuper2040";
    hname = (h).c_str();	
    rholeadingJetsuper2040[ik] = new TH1D(Form("rholeadingJetsuper2040%d",ik), "", ndJRsuper, 0.0, JRrangesuper);
	
	h = "rholeadingJetsuper3040";
    hname = (h).c_str();	
    rholeadingJetsuper3040[ik] = new TH1D(Form("rholeadingJetsuper3040%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuper4080";
    hname = (h).c_str();	
    rholeadingJetsuper4080[ik] = new TH1D(Form("rholeadingJetsuper4080%d",ik), "", ndJRsuper, 0.0, JRrangesuper);

	h = "rholeadingJetsuper80";
    hname = (h).c_str();	
    rholeadingJetsuper80[ik] = new TH1D(Form("rholeadingJetsuper80%d",ik), "", ndJRsuper, 0.0, JRrangesuper);
	
	
	
    h = "SetaphileadingJet";
    hname = (h).c_str();	
    SetaphileadingJet[ik] = new TH2D(Form("SetaphileadingJet%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);	

    h = "SetaphileadingJetCAT";
    hname = (h).c_str();	
    SetaphileadingJetCAT[ik] = new TH2D(Form("SetaphileadingJetCAT%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);		
	
	
    h = "SPTetaphileadingJet";
    hname = (h).c_str();	
    SPTetaphileadingJet[ik] = new TH2D(Form("SPTetaphileadingJet%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);	

    h = "SPTetaphileadingJetCAT";
    hname = (h).c_str();	
    SPTetaphileadingJetCAT[ik] = new TH2D(Form("SPTetaphileadingJetCAT%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);	
	
    h = "SPTetaphileadingJet0510";
    hname = (h).c_str();	
    SPTetaphileadingJet0510[ik] = new TH2D(Form("SPTetaphileadingJet0510%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet1020";
    hname = (h).c_str();	
    SPTetaphileadingJet1020[ik] = new TH2D(Form("SPTetaphileadingJet1020%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet2030";
    hname = (h).c_str();	
    SPTetaphileadingJet2030[ik] = new TH2D(Form("SPTetaphileadingJet2030%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet2040";
    hname = (h).c_str();	
    SPTetaphileadingJet2040[ik] = new TH2D(Form("SPTetaphileadingJet2040%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet3040";
    hname = (h).c_str();	
    SPTetaphileadingJet3040[ik] = new TH2D(Form("SPTetaphileadingJet3040%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet4080";
    hname = (h).c_str();	
    SPTetaphileadingJet4080[ik] = new TH2D(Form("SPTetaphileadingJet4080%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

    h = "SPTetaphileadingJet80";
    hname = (h).c_str();	
    SPTetaphileadingJet80[ik] = new TH2D(Form("SPTetaphileadingJet80%d",ik), "", ndJphi, -1.0/2.0*SJPhirange, 3.0/2.0*SJPhirange, ndJeta, -SJEtarange, SJEtarange);

	}
	
	
//..........energy flow

	ntistep=ntistep + 1;	

    TH1D *dndphiJet01[ntistep];			
    TH1D *dndphiJet12[ntistep];			
    TH1D *dndphiJet24[ntistep];
    TH1D *dndphiJet40[ntistep];

    TH1D *dEtdphi01[ntistep];			
    TH1D *dEtdphi12[ntistep];			
    TH1D *dEtdphi24[ntistep];
    TH1D *dEtdphi40[ntistep];
			
    TH1D *dpldphi01[ntistep];			
    TH1D *dpldphi12[ntistep];			
    TH1D *dpldphi24[ntistep];
    TH1D *dpldphi40[ntistep];		

	ntistep=ntistep - 1;		

    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	
    h = "dndphiJet01";
    hname = (h).c_str();
    dndphiJet01[ik] = new TH1D(Form("dndphiJet01%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dndphiJet12";
    hname = (h).c_str();
    dndphiJet12[ik] = new TH1D(Form("dndphiJet12%d",ik), "", ndJphi, 0.0, JPhirange);			
			
    h = "dndphiJet24";
    hname = (h).c_str();
    dndphiJet24[ik] = new TH1D(Form("dndphiJet24%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dndphiJet40";
    hname = (h).c_str();	
    dndphiJet40[ik] = new TH1D(Form("dndphiJet40%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dEtdphi01";
    hname = (h).c_str();	
    dEtdphi01[ik] = new TH1D(Form("dEtdphi01%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dEtdphi12";
    hname = (h).c_str();	
    dEtdphi12[ik] = new TH1D(Form("dEtdphi12%d",ik), "", ndJphi, 0.0, JPhirange);			
			
    h = "dEtdphi24";
    hname = (h).c_str();	
    dEtdphi24[ik] = new TH1D(Form("dEtdphi24%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dEtdphi40";
    hname = (h).c_str();	
    dEtdphi40[ik] = new TH1D(Form("dEtdphi40%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dpldphi01";
    hname = (h).c_str();	
    dpldphi01[ik] = new TH1D(Form("dpldphi01%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dpldphi12";
    hname = (h).c_str();	
    dpldphi12[ik] = new TH1D(Form("dpldphi12%d",ik), "", ndJphi, 0.0, JPhirange);			
			
    h = "dpldphi24";
    hname = (h).c_str();	
    dpldphi24[ik] = new TH1D(Form("dpldphi24%d",ik), "", ndJphi, 0.0, JPhirange);

    h = "dpldphi40";
    hname = (h).c_str();	
    dpldphi40[ik] = new TH1D(Form("dpldphi40%d",ik), "", ndJphi, 0.0, JPhirange);
	}
			
//..........energy flow


	ntistep=ntistep + 1;	

    TH1D *fractionAj[ntistep];

	ntistep=ntistep - 1;		
	
	double AjRange=1.0;
	int nAj=20;

    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	
    h = "fractionAj";
    hname = (h).c_str();	
    fractionAj[ik] = new TH1D(Form("fractionAj%d",ik), "", nAj, 0.0, AjRange);   
	}

	ntistep=ntistep + 1;		
	
    TH1D *ptllAj[ntistep];

	ntistep=ntistep - 1;
	
	double ptllAjRange=0.5;
	int nptllAj=5;

    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	
    h = "ptllAj";
    hname = (h).c_str();	
    ptllAj[ik] = new TH1D(Form("ptllAj%d",ik), "", nptllAj, 0.0, ptllAjRange); 	
	}

	
	
	ntistep=ntistep + 1;	
	
    TH1D *ptlldeltaR[ntistep];
    TH1D *ptlldeltaR0510[ntistep];
    TH1D *ptlldeltaR1020[ntistep];	
    TH1D *ptlldeltaR2030[ntistep];	
    TH1D *ptlldeltaR2040[ntistep];	
    TH1D *ptlldeltaR3040[ntistep];	
    TH1D *ptlldeltaR4080[ntistep];	
    TH1D *ptlldeltaR80[ntistep];

	ntistep=ntistep - 1;	
	
	double deltaRRange=1.0;
	int ndeltaR=10;


    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {
	
    h = "ptlldeltaR";
    hname = (h).c_str();	
    ptlldeltaR[ik] = new TH1D(Form("ptlldeltaR%d",ik), "", ndeltaR, 0.0, deltaRRange); 	

    h = "ptlldeltaR0510";
    hname = (h).c_str();	
    ptlldeltaR0510[ik] = new TH1D(Form("ptlldeltaR0510%d",ik), "", ndeltaR, 0.0, deltaRRange); 

    h = "ptlldeltaR1020";
    hname = (h).c_str();	
    ptlldeltaR1020[ik] = new TH1D(Form("ptlldeltaR1020%d",ik), "", ndeltaR, 0.0, deltaRRange); 
	
    h = "ptlldeltaR2030";
    hname = (h).c_str();	
    ptlldeltaR2030[ik] = new TH1D(Form("ptlldeltaR2030%d",ik), "", ndeltaR, 0.0, deltaRRange); 

    h = "ptlldeltaR2040";
    hname = (h).c_str();	
    ptlldeltaR2040[ik] = new TH1D(Form("ptlldeltaR2040%d",ik), "", ndeltaR, 0.0, deltaRRange); 	
	
    h = "ptlldeltaR3040";
    hname = (h).c_str();	
    ptlldeltaR3040[ik] = new TH1D(Form("ptlldeltaR3040%d",ik), "", ndeltaR, 0.0, deltaRRange); 

    h = "ptlldeltaR4080";
    hname = (h).c_str();	
    ptlldeltaR4080[ik] = new TH1D(Form("ptlldeltaR4080%d",ik), "", ndeltaR, 0.0, deltaRRange); 

    h = "ptlldeltaR80";
    hname = (h).c_str();	
    ptlldeltaR80[ik] = new TH1D(Form("ptlldeltaR80%d",ik), "", ndeltaR, 0.0, deltaRRange); 	

    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH1D *gammajetfraction;

    h = "gammajetfraction";
    hname = (h).c_str();	
    gammajetfraction = new TH1D(hname, "", 10, 0.0, 10.0); 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	
//...Histogram weight for fastjet
	
    double wwJ10 = 1.0 / LBT->ncall;

    double wwJdndetadphi = ndJeta * ndJphi / JEtarange / JPhirange / LBT->ncall;

    double wwJdndptde = ndJpt * ndJe / JPtrange / JErange / LBT->ncall;

    double wwJdndthetade = ndJtheta * ndJe / JThetarange / JErange / LBT->ncall;

    double wwJdndphide = ndJphi * ndJe / JPhirange / JErange / LBT->ncall;

    double wwJdnde = ndJe / JErange / LBT->ncall;

    double wwJdndpt = ndJpt / JPtrange / LBT->ncall;

    double wwJdndtheta = ndJtheta / JThetarange / LBT->ncall;

    double wwJdndphi = ndJphi / JPhirange / LBT->ncall;

    double wwJdndeta = ndJeta / JEtarange / LBT->ncall;	

    double wwdR = ndJR / JRrange / LBT->ncall;

    double wwdRsuper = ndJRsuper / JRrangesuper / LBT->ncall;
	
    double wwdndxpl = ndndxpl / Jxrange /LBT->ncall;
	
    double wwdndxpl0 = ndndxpl0 / Jx0range /LBT->ncall;
	
//............................................................................................Histograms setting end.





//............................................................................................Histograms for LBT transport

//................................................................histogram weight	
    double ww10 = 1.0/LBT->ncall;	
	
    double wwdndrdz = ndr * ndz / Rrange / Zrange / LBT->ncall;

    double wwdndetadphi = ndeta * ndphi / Etarange / Phirange / LBT->ncall;

	int nmix0=1;
	
    //double wwSdndetadphi = ndJeta * ndJphi / 2.0*SJEtarange / 2.0*SJPhirange / LBT->ncall / nmix0;		

    double wwSdndetadphi = 1.0/ LBT->ncall / nmix0;			
	
    double wwdndptde = ndpt * nde / Ptrange / Erange / LBT->ncall;

    double wwdndthetade = ndtheta * (nde*10) / Thetarange / Erange / LBT->ncall;

    double wwdndphide = ndphi * nde / Phirange / Erange / LBT->ncall;

    double wwdnde = nde / Erange / LBT->ncall;

    double wwdndpt = ndpt / Ptrange / LBT->ncall;

    double wwdndtheta = ndtheta / Thetarange / LBT->ncall;

    double wwdndphi = ndphi / Phirange / LBT->ncall;

    double wwdndeta = ndeta / Etarange / LBT->ncall;	
	
	double wwdndetade = ndeta * nde / Etarange / Erange / LBT->ncall;
//................................................................histogram weight end
	






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    TH1D *Sleadingpartonenergy;
    TH1D *Sleadingpartonpt;
    TH1D *Sleadingpartonpesudopt;	
    TH1D *Sleadingpartonpesudopt2;
    TH1D *Sleadingpartonqhat;	
    TH1D *Sleadingelossf;
    TH1D *Sleadingeloss;
    TH1D *Sleadingtheta;
    TH1D *Sleadingphi;
    TH1D *Sleadingeta;

    h = "Sleadingpartonenergy";
    hname = (h).c_str();
    Sleadingpartonenergy = new TH1D(hname, "", 100, -0.05, 9.95);

    h = "Sleadingpartonpt";
    hname = (h).c_str();
    Sleadingpartonpt = new TH1D(hname, "", 100, -0.05, 9.95);

    h = "Sleadingpartonpesudopt";
    hname = (h).c_str();
    Sleadingpartonpesudopt = new TH1D(hname, "", 100, -0.05, 9.95);

    h = "Sleadingpartonpesudopt2";
    hname = (h).c_str();
    Sleadingpartonpesudopt2 = new TH1D(hname, "", 100, -0.05, 9.95);

    h = "Sleadingpartonqhat";
    hname = (h).c_str();
    Sleadingpartonqhat = new TH1D(hname, "", 100, -0.05, 9.95);
	
    h = "Sleadingelossf";
    hname = (h).c_str();
    Sleadingelossf = new TH1D(hname, "", 100, -0.05, 9.95);	
	
    h = "Sleadingeloss";
    hname = (h).c_str();
    Sleadingeloss = new TH1D(hname, "", 100, -0.05, 9.95);	
	
    h = "Sleadingtheta";
    hname = (h).c_str();
    Sleadingtheta = new TH1D(hname, "", 100, -0.05, 9.95);	
	
    h = "Sleadingphi";
    hname = (h).c_str();
    Sleadingphi = new TH1D(hname, "", 100, -0.05, 9.95);
	
    h = "Sleadingeta";
    hname = (h).c_str();
    Sleadingeta = new TH1D(hname, "", 100, -0.05, 9.95);




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////





	
//................................................................histogram booking	
    TH1D *leadingpartonenergy;
    TH1D *leadingpartonpt;
    TH1D *leadingpartonpesudopt;	
    TH1D *leadingpartonpesudopt2;	
    TH1D *leadingelossf;
    TH1D *leadingeloss;
    TH1D *leadingtheta;
    TH1D *leadingphi;
    TH1D *leadingeta;

    h = "leadingpartonenergy";
    hname = (h).c_str();
    leadingpartonenergy = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingpartonpt";
    hname = (h).c_str();
    leadingpartonpt = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingpartonpesudopt";
    hname = (h).c_str();
    leadingpartonpesudopt = new TH1D(hname, "", ntistep, 0.0, LBT->timend);

    h = "leadingpartonpesudopt2";
    hname = (h).c_str();
    leadingpartonpesudopt2 = new TH1D(hname, "", ntistep, 0.0, LBT->timend);
	
    h = "leadingelossf";
    hname = (h).c_str();
    leadingelossf = new TH1D(hname, "", ntistep, 0.0, LBT->timend);	
	
    h = "leadingeloss";
    hname = (h).c_str();
    leadingeloss = new TH1D(hname, "", ntistep, 0.0, LBT->timend);	
	
    h = "leadingtheta";
    hname = (h).c_str();
    leadingtheta = new TH1D(hname, "", ntistep, 0.0, LBT->timend);	
	
    h = "leadingphi";
    hname = (h).c_str();
    leadingphi = new TH1D(hname, "", ntistep, 0.0, LBT->timend);
	
    h = "leadingeta";
    hname = (h).c_str();
    leadingeta = new TH1D(hname, "", ntistep, 0.0, LBT->timend);		

    TH2D *dndedtleading0;
    TH2D *dndptdtleading0;
    TH2D *dndpt2dtleading0;			
    TH2D *dndthetadtleading0;
    TH2D *dndphidtleading0;	
    TH2D *dndetadtleading0;

    h = "dndedtleading0";
    hname = (h).c_str();	
    dndedtleading0 = new TH2D(hname, "", nde, 0.0, Erange, ntistep, 0.0, LBT->timend);
	
    h = "dndptdtleading0";	
    hname = (h).c_str();
    dndptdtleading0 = new TH2D(hname, "", ndpt, 0.0, Ptrange, ntistep, 0.0, LBT->timend);
	
    h = "dndpt2dtleading0";	
    hname = (h).c_str();		
    dndpt2dtleading0 = new TH2D(hname, "", ndpt, 0.0, Ptrange, ntistep, 0.0, LBT->timend);
	
    h = "dndthetadtleading0";	
    hname = (h).c_str();	
    dndthetadtleading0 = new TH2D(hname, "", ndtheta, 0.0, Thetarange, ntistep, 0.0, LBT->timend);
	
    h = "dndphidtleading0";	
    hname = (h).c_str();	
    dndphidtleading0 = new TH2D(hname, "", ndphi, -Phirange, Phirange, ntistep, 0.0, LBT->timend);
	
    h = "dndetadtleading0";	
    hname = (h).c_str();	
    dndetadtleading0 = new TH2D(hname, "", ndeta, -Etarange, Etarange, ntistep, 0.0, LBT->timend);






///////////////////////////////////////////////////////////////...321

    TH1D *CATdndtheta;
    TH1D *CATdnde;	
    TH1D *CATdndpt;
    TH1D *CATdndpt2;	
    TH1D *CATdndeta;
    TH1D *CATdndphi;

    TH1D *CATdndtheta0010;
    TH1D *CATdndtheta1020;
    TH1D *CATdndtheta2040;	
    TH1D *CATdndtheta4080;	
    TH1D *CATdndtheta80;
	
	
    TH2D *CATdndthetade;	

    TH1D *RADdndtheta;
    TH1D *RADdnde;	
    TH1D *RADdndpt;
    TH1D *RADdndpt2;	
    TH1D *RADdndeta;
    TH1D *RADdndphi;

    TH1D *RADdndtheta0010;
    TH1D *RADdndtheta1020;
    TH1D *RADdndtheta2040;	
    TH1D *RADdndtheta4080;	
    TH1D *RADdndtheta80;
	

    TH2D *RADdndthetade;
	

    TH1D *NEGdndtheta;
    TH1D *NEGdnde;	
    TH1D *NEGdndpt;
    TH1D *NEGdndpt2;	
    TH1D *NEGdndeta;
    TH1D *NEGdndphi;

    TH1D *NEGdndtheta0010;
    TH1D *NEGdndtheta1020;
    TH1D *NEGdndtheta2040;	
    TH1D *NEGdndtheta4080;	
    TH1D *NEGdndtheta80;	


    TH2D *NEGdndthetade;	


    TH1D *LEADdndtheta;
    TH1D *LEADdnde;	
    TH1D *LEADdndpt;
    TH1D *LEADdndpt2;	
    TH1D *LEADdndeta;
    TH1D *LEADdndphi;

    TH1D *LEADdndtheta0010;
    TH1D *LEADdndtheta1020;
    TH1D *LEADdndtheta2040;	
    TH1D *LEADdndtheta4080;	
    TH1D *LEADdndtheta80;	


    TH2D *LEADdndthetade;

	
    h = "CATdnde";
    hname = (h).c_str();	
    CATdnde = new TH1D(hname, "", nde, 0.0, Erange);	
	
    h = "CATdndtheta";
    hname = (h).c_str();	
    CATdndtheta = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	
	
    h = "CATdndpt";
    hname = (h).c_str();	
    CATdndpt = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "CATdndpt2";
    hname = (h).c_str();	
    CATdndpt2 = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "CATdndeta";
    hname = (h).c_str();	
    CATdndeta = new TH1D(hname, "", ndeta, 0.0, Etarange);	
	
    h = "CATdndtheta0010";
    hname = (h).c_str();	
    CATdndtheta0010 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "CATdndtheta1020";
    hname = (h).c_str();	
    CATdndtheta1020 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "CATdndtheta2040";
    hname = (h).c_str();	
    CATdndtheta2040 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "CATdndtheta4080";
    hname = (h).c_str();	
    CATdndtheta4080 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "CATdndtheta80";
    hname = (h).c_str();	
    CATdndtheta80 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);
	
    h = "CATdndphi";
    hname = (h).c_str();	
    CATdndphi = new TH1D(hname, "", ndphi, 0.0, Phirange);	


    h = "RADdnde";
    hname = (h).c_str();	
    RADdnde = new TH1D(hname, "", nde, 0.0, Erange);	
	
    h = "RADdndtheta";
    hname = (h).c_str();	
    RADdndtheta = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "RADdndtheta0010";
    hname = (h).c_str();	
    RADdndtheta0010 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "RADdndtheta1020";
    hname = (h).c_str();	
    RADdndtheta1020 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "RADdndtheta2040";
    hname = (h).c_str();	
    RADdndtheta2040 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "RADdndtheta4080";
    hname = (h).c_str();	
    RADdndtheta4080 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "RADdndtheta80";
    hname = (h).c_str();	
    RADdndtheta80 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);

	
    h = "RADdndpt";
    hname = (h).c_str();	
    RADdndpt = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "RADdndpt2";
    hname = (h).c_str();	
    RADdndpt2 = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "RADdndeta";
    hname = (h).c_str();	
    RADdndeta = new TH1D(hname, "", ndeta, 0.0, Etarange);	

    h = "RADdndphi";
    hname = (h).c_str();	
    RADdndphi = new TH1D(hname, "", ndphi, 0.0, Phirange);	




    h = "NEGdnde";
    hname = (h).c_str();	
    NEGdnde = new TH1D(hname, "", nde, 0.0, Erange);	
	
    h = "NEGdndtheta";
    hname = (h).c_str();	
    NEGdndtheta = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "NEGdndtheta0010";
    hname = (h).c_str();	
    NEGdndtheta0010 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "NEGdndtheta1020";
    hname = (h).c_str();	
    NEGdndtheta1020 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "NEGdndtheta2040";
    hname = (h).c_str();	
    NEGdndtheta2040 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "NEGdndtheta4080";
    hname = (h).c_str();	
    NEGdndtheta4080 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "NEGdndtheta80";
    hname = (h).c_str();	
    NEGdndtheta80 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);

	
    h = "NEGdndpt";
    hname = (h).c_str();	
    NEGdndpt = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "NEGdndpt2";
    hname = (h).c_str();	
    NEGdndpt2 = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "NEGdndeta";
    hname = (h).c_str();	
    NEGdndeta = new TH1D(hname, "", ndeta, 0.0, Etarange);	

    h = "NEGdndphi";
    hname = (h).c_str();	
    NEGdndphi = new TH1D(hname, "", ndphi, 0.0, Phirange);



    h = "LEADdnde";
    hname = (h).c_str();	
    LEADdnde = new TH1D(hname, "", nde, 0.0, Erange);	
	
    h = "LEADdndtheta";
    hname = (h).c_str();	
    LEADdndtheta = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "LEADdndtheta0010";
    hname = (h).c_str();	
    LEADdndtheta0010 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "LEADdndtheta1020";
    hname = (h).c_str();	
    LEADdndtheta1020 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "LEADdndtheta2040";
    hname = (h).c_str();	
    LEADdndtheta2040 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "LEADdndtheta4080";
    hname = (h).c_str();	
    LEADdndtheta4080 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);	

    h = "LEADdndtheta80";
    hname = (h).c_str();	
    LEADdndtheta80 = new TH1D(hname, "", ndtheta, 0.0, Thetarange);

	
    h = "LEADdndpt";
    hname = (h).c_str();	
    LEADdndpt = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "LEADdndpt2";
    hname = (h).c_str();	
    LEADdndpt2 = new TH1D(hname, "", ndpt, 0.0, Ptrange);	
	
    h = "LEADdndeta";
    hname = (h).c_str();	
    LEADdndeta = new TH1D(hname, "", ndeta, 0.0, Etarange);	

    h = "LEADdndphi";
    hname = (h).c_str();	
    LEADdndphi = new TH1D(hname, "", ndphi, 0.0, Phirange);













	

    h = "CATdndthetade";
    hname = (h).c_str();	
    CATdndthetade = new TH2D(hname, "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);


    h = "RADdndthetade";
    hname = (h).c_str();	
    RADdndthetade = new TH2D(hname, "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);


    h = "NEGdndthetade";
    hname = (h).c_str();	
    NEGdndthetade = new TH2D(hname, "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);


    h = "LEADdndthetade";
    hname = (h).c_str();	
    LEADdndthetade = new TH2D(hname, "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);

	
///////////////////////////////////////////////////////////////














//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

    TH1D *dndeleading[ntistep];
    TH1D *dndptleading[ntistep];
    TH1D *dndpt2leading[ntistep];			
    TH1D *dndthetaleading[ntistep];
    TH1D *dndphileading[ntistep];	
    TH1D *dndetaleading[ntistep];

    TH1D *dndpesudoptleading[ntistep];
    TH1D *dndpesudopt2leading[ntistep];

	
    TH2D *rdedrdz[ntistep];
    TH2D *dedrdz[ntistep];
    TH2D *rdndrdz[ntistep];
    TH2D *dndrdz[ntistep];

    TH2D *rdedrdzpos[ntistep];
    TH2D *dedrdzpos[ntistep];
    TH2D *rdndrdzpos[ntistep];
    TH2D *dndrdzpos[ntistep];

    TH2D *rdedrdzneg[ntistep];
    TH2D *dedrdzneg[ntistep];
    TH2D *rdndrdzneg[ntistep];
    TH2D *dndrdzneg[ntistep];

    TH2D *rdedrdzCAT[ntistep];
    TH2D *dedrdzCAT[ntistep];
    TH2D *rdndrdzCAT[ntistep];
    TH2D *dndrdzCAT[ntistep];

    TH2D *rdedrdzRAD[ntistep];
    TH2D *dedrdzRAD[ntistep];
    TH2D *rdndrdzRAD[ntistep];
    TH2D *dndrdzRAD[ntistep];

    TH2D *rdedxdz[ntistep];
    TH2D *dedxdz[ntistep];
    TH2D *rdndxdz[ntistep];
    TH2D *dndxdz[ntistep];

    TH2D *rdedxdzpos[ntistep];
    TH2D *dedxdzpos[ntistep];
    TH2D *rdndxdzpos[ntistep];
    TH2D *dndxdzpos[ntistep];

    TH2D *rdedxdzneg[ntistep];
    TH2D *dedxdzneg[ntistep];
    TH2D *rdndxdzneg[ntistep];
    TH2D *dndxdzneg[ntistep];		

    TH2D *rdedxdzCAT[ntistep];
    TH2D *dedxdzCAT[ntistep];
    TH2D *rdndxdzCAT[ntistep];
    TH2D *dndxdzCAT[ntistep];

    TH2D *rdedxdzRAD[ntistep];
    TH2D *dedxdzRAD[ntistep];
    TH2D *rdndxdzRAD[ntistep];
    TH2D *dndxdzRAD[ntistep];
	
    TH2D *dedetadphi[ntistep];
    TH2D *dedetadphipos[ntistep];
    TH2D *dedetadphineg[ntistep];
    TH2D *dedetadphiCAT[ntistep];
    TH2D *dedetadphiRAD[ntistep];	
    TH2D *dndetadphi[ntistep];
    TH2D *dndetadphipos[ntistep];
    TH2D *dndetadphineg[ntistep];
    TH2D *dndetadphiCAT[ntistep];
    TH2D *dndetadphiRAD[ntistep];


    TH2D *dndptde[ntistep];
    TH2D *dndpt2de[ntistep];			
    TH2D *dndthetade[ntistep];
    TH2D *dndphide[ntistep];	
    TH2D *dndetade[ntistep];			

    TH2D *dndptdepos[ntistep];
    TH2D *dndpt2depos[ntistep];			
    TH2D *dndthetadepos[ntistep];
    TH2D *dndphidepos[ntistep];	
    TH2D *dndetadepos[ntistep];

    TH2D *dndptdeneg[ntistep];
    TH2D *dndpt2deneg[ntistep];			
    TH2D *dndthetadeneg[ntistep];
    TH2D *dndphideneg[ntistep];	
    TH2D *dndetadeneg[ntistep];

    TH2D *dndptdeCAT[ntistep];
    TH2D *dndpt2deCAT[ntistep];			
    TH2D *dndthetadeCAT[ntistep];
    TH2D *dndphideCAT[ntistep];	
    TH2D *dndetadeCAT[ntistep];

    TH2D *dndptdeRAD[ntistep];
    TH2D *dndpt2deRAD[ntistep];			
    TH2D *dndthetadeRAD[ntistep];
    TH2D *dndphideRAD[ntistep];	
    TH2D *dndetadeRAD[ntistep];

    TH1D *dnde[ntistep];
    TH1D *dndpt[ntistep];
    TH1D *dndpt2[ntistep];
    TH1D *dndpesudopt[ntistep];
    TH1D *dndpesudopt2[ntistep];	
    TH1D *dndtheta[ntistep];
    TH1D *dndphi[ntistep];	
    TH1D *dndeta[ntistep];			

    TH1D *dndphi01[ntistep];	
    TH1D *dndphi12[ntistep];
    TH1D *dndphi23[ntistep];
    TH1D *dndphi34[ntistep];	
    TH1D *dndphi24[ntistep];
    TH1D *dndphi48[ntistep];
    TH1D *dndphi80[ntistep];	
	
    TH1D *dndepos[ntistep];
    TH1D *dndptpos[ntistep];
    TH1D *dndpt2pos[ntistep];
    TH1D *dndpesudoptpos[ntistep];
    TH1D *dndpesudopt2pos[ntistep];	
    TH1D *dndthetapos[ntistep];
    TH1D *dndphipos[ntistep];	
    TH1D *dndetapos[ntistep];

    TH1D *dndeneg[ntistep];
    TH1D *dndptneg[ntistep];
    TH1D *dndpt2neg[ntistep];
    TH1D *dndpesudoptneg[ntistep];
    TH1D *dndpesudopt2neg[ntistep];	
    TH1D *dndthetaneg[ntistep];
    TH1D *dndphineg[ntistep];	
    TH1D *dndetaneg[ntistep];	

    TH1D *dndeCAT[ntistep];
    TH1D *dndptCAT[ntistep];
    TH1D *dndpt2CAT[ntistep];
    TH1D *dndpesudoptCAT[ntistep];
    TH1D *dndpesudopt2CAT[ntistep];	
    TH1D *dndthetaCAT[ntistep];
    TH1D *dndphiCAT[ntistep];	
    TH1D *dndetaCAT[ntistep];

    TH1D *dndeRAD[ntistep];
    TH1D *dndptRAD[ntistep];
    TH1D *dndpt2RAD[ntistep];
    TH1D *dndpesudoptRAD[ntistep];
    TH1D *dndpesudopt2RAD[ntistep];	
    TH1D *dndthetaRAD[ntistep];
    TH1D *dndphiRAD[ntistep];	
    TH1D *dndetaRAD[ntistep];



	
    for (unsigned ik = 1; ik <= ntistep; ++ik)
    {	

    h = "dndeleading";    	
    dndeleading[ik] = new TH1D(Form("dndeleading%d",ik), "", nde, 0.0, Erange);

    h = "dndptleading";
    dndptleading[ik] = new TH1D(Form("dndptleading%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2leading";    	
    dndpt2leading[ik] = new TH1D(Form("dndpt2leading%d",ik), "", ndpt, 0.0, Ptrange);	
 
     h = "dndpesudoptleading";
    dndpesudoptleading[ik] = new TH1D(Form("dndpesudoptleading%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2leading";    	
    dndpesudopt2leading[ik] = new TH1D(Form("dndpesudopt2leading%d",ik), "", ndpt, 0.0, Ptrange);	
 
    h = "dndthetaleading";    	
    dndthetaleading[ik] = new TH1D(Form("dndthetaleading%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphileading";    	
    dndphileading[ik] = new TH1D(Form("dndphileading%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndetaleading";    	
    dndetaleading[ik] = new TH1D(Form("dndetaleading%d",ik), "", ndeta, -Etarange, Etarange);			
			
			
	
    h = "rdedrdz";    	
    rdedrdz[ik] = new TH2D(Form("rdedrdz%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "dedrdz";    	
    dedrdz[ik] = new TH2D(Form("dedrdz%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "rdndrdz";    	 
    rdndrdz[ik] = new TH2D(Form("rdndrdz%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndrdz";    	
    dndrdz[ik] = new TH2D(Form("dndrdz%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	

    h = "rdedrdzpos";    	
    rdedrdzpos[ik] = new TH2D(Form("rdedrdzpos%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "dedrdzpos";    	
    dedrdzpos[ik] = new TH2D(Form("dedrdzpos%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "rdndrdzpos";    	 
    rdndrdzpos[ik] = new TH2D(Form("rdndrdzpos%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndrdzpos";    	
    dndrdzpos[ik] = new TH2D(Form("dndrdzpos%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);	

    h = "rdedrdzneg";    	
    rdedrdzneg[ik] = new TH2D(Form("rdedrdzneg%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "dedrdzneg";    	
    dedrdzneg[ik] = new TH2D(Form("dedrdzneg%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "rdndrdzneg";    	 
    rdndrdzneg[ik] = new TH2D(Form("rdndrdzneg%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndrdzneg";    	
    dndrdzneg[ik] = new TH2D(Form("dndrdzneg%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
		

    h = "rdedrdzCAT";    	
    rdedrdzCAT[ik] = new TH2D(Form("rdedrdzCAT%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "dedrdzCAT";    	
    dedrdzCAT[ik] = new TH2D(Form("dedrdzCAT%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "rdndrdzCAT";    	 
    rdndrdzCAT[ik] = new TH2D(Form("rdndrdzCAT%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndrdzCAT";    	
    dndrdzCAT[ik] = new TH2D(Form("dndrdzCAT%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);	

    h = "rdedrdzRAD";    	
    rdedrdzRAD[ik] = new TH2D(Form("rdedrdzRAD%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "dedrdzRAD";    	
    dedrdzRAD[ik] = new TH2D(Form("dedrdzRAD%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);

    h = "rdndrdzRAD";    	 
    rdndrdzRAD[ik] = new TH2D(Form("rdndrdzRAD%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndrdzRAD";    	
    dndrdzRAD[ik] = new TH2D(Form("dndrdzRAD%d",ik), "", ndr, 0.0, Rrange, ndz, -Zrange, Zrange);






////////////////////////////////////////////////////////////////////////////

    h = "rdedxdz";    	
    rdedxdz[ik] = new TH2D(Form("rdedxdz%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "dedxdz";    	
    dedxdz[ik] = new TH2D(Form("dedxdz%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "rdndxdz";    	 
    rdndxdz[ik] = new TH2D(Form("rdndxdz%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndxdz";    	
    dndxdz[ik] = new TH2D(Form("dndxdz%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	

    h = "rdedxdzpos";    	
    rdedxdzpos[ik] = new TH2D(Form("rdedxdzpos%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "dedxdzpos";    	
    dedxdzpos[ik] = new TH2D(Form("dedxdzpos%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "rdndxdzpos";    	 
    rdndxdzpos[ik] = new TH2D(Form("rdndxdzpos%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndxdzpos";    	
    dndxdzpos[ik] = new TH2D(Form("dndxdzpos%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);	

    h = "rdedxdzneg";    	
    rdedxdzneg[ik] = new TH2D(Form("rdedxdzneg%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "dedxdzneg";    	
    dedxdzneg[ik] = new TH2D(Form("dedxdzneg%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "rdndxdzneg";    	 
    rdndxdzneg[ik] = new TH2D(Form("rdndxdzneg%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndxdzneg";    	
    dndxdzneg[ik] = new TH2D(Form("dndxdzneg%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);


    h = "rdedxdzCAT";    	
    rdedxdzCAT[ik] = new TH2D(Form("rdedxdzCAT%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "dedxdzCAT";    	
    dedxdzCAT[ik] = new TH2D(Form("dedxdzCAT%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "rdndxdzCAT";    	 
    rdndxdzCAT[ik] = new TH2D(Form("rdndxdzCAT%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndxdzCAT";    	
    dndxdzCAT[ik] = new TH2D(Form("dndxdzCAT%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);	

    h = "rdedxdzRAD";    	
    rdedxdzRAD[ik] = new TH2D(Form("rdedxdzRAD%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "dedxdzRAD";    	
    dedxdzRAD[ik] = new TH2D(Form("dedxdzRAD%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

    h = "rdndxdzRAD";    	 
    rdndxdzRAD[ik] = new TH2D(Form("rdndxdzRAD%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);
	
    h = "dndxdzRAD";    	
    dndxdzRAD[ik] = new TH2D(Form("dndxdzRAD%d",ik), "", ndr, -Rrange, Rrange, ndz, -Zrange, Zrange);

	
	 
    h = "dedetadphi";    	
    dedetadphi[ik] = new TH2D(Form("dedetadphi%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);

    h = "dndetadphi";    	
    dndetadphi[ik] = new TH2D(Form("dndetadphi%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);
	

    h = "dedetadphipos";    	
    dedetadphipos[ik] = new TH2D(Form("dedetadphipos%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);

    h = "dndetadphipos";    	
    dndetadphipos[ik] = new TH2D(Form("dndetadphipos%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange); 
	 
    h = "dedetadphineg";    	
    dedetadphineg[ik] = new TH2D(Form("dedetadphineg%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);

    h = "dndetadphineg";    
    dndetadphineg[ik] = new TH2D(Form("dndetadphineg%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange); 

    h = "dedetadphiCAT";    	
    dedetadphiCAT[ik] = new TH2D(Form("dedetadphiCAT%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);

    h = "dndetadphiCAT";    	
    dndetadphiCAT[ik] = new TH2D(Form("dndetadphiCAT%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange); 
	 
    h = "dedetadphiRAD";    	
    dedetadphiRAD[ik] = new TH2D(Form("dedetadphiRAD%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);

    h = "dndetadphiRAD";    
    dndetadphiRAD[ik] = new TH2D(Form("dndetadphiRAD%d",ik), "", ndphi, -Phirange, Phirange, ndeta, -Etarange, Etarange);
	
    h = "dndptde";    	
    dndptde[ik] = new TH2D(Form("dndptde%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);
	
    h = "dndpt2de";    	
    dndpt2de[ik] = new TH2D(Form("dndpt2de%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);	
 
    h = "dndthetade";    	
    dndthetade[ik] = new TH2D(Form("dndthetade%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);	 	 
	 
    h = "dndphide";    	
    dndphide[ik] = new TH2D(Form("dndphide%d",ik), "", ndphi, -Phirange, Phirange, nde, 0.0, Erange);

    h = "dndetade";    	
    dndetade[ik] = new TH2D(Form("dndetade%d",ik), "", ndeta, -Etarange, Etarange, nde, 0.0, Erange);
	

    h = "dndptdepos";    	
    dndptdepos[ik] = new TH2D(Form("dndptdepos%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);
	
    h = "dndpt2depos";    	
    dndpt2depos[ik] = new TH2D(Form("dndpt2depos%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);	
 
    h = "dndthetadepos";    	
    dndthetadepos[ik] = new TH2D(Form("dndthetadepos%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);	 	 
	 
    h = "dndphidepos";    	
    dndphidepos[ik] = new TH2D(Form("dndphidepos%d",ik), "", ndphi, -Phirange, Phirange, nde, 0.0, Erange);

    h = "dndetadepos";    	
    dndetadepos[ik] = new TH2D(Form("dndetadepos%d",ik), "", ndeta, Etarange, Etarange, nde, 0.0, Erange);

    h = "dndptdeneg";    	
    dndptdeneg[ik] = new TH2D(Form("dndptdeneg%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);
	
    h = "dndpt2deneg";    	
    dndpt2deneg[ik] = new TH2D(Form("dndpt2deneg%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);	
 
    h = "dndthetadeneg";    	
    dndthetadeneg[ik] = new TH2D(Form("dndthetadeneg%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);	 	 
	 
    h = "dndphideneg";    	
    dndphideneg[ik] = new TH2D(Form("dndphideneg%d",ik), "", ndphi, -Phirange, Phirange, nde, 0.0, Erange);

    h = "dndetadeneg";    	
    dndetadeneg[ik] = new TH2D(Form("dndetadeneg%d",ik), "", ndeta, Etarange, Etarange, nde, 0.0, Erange);
	

    h = "dndptdeCAT";    	
    dndptdeCAT[ik] = new TH2D(Form("dndptdeCAT%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);
	
    h = "dndpt2deCAT";    	
    dndpt2deCAT[ik] = new TH2D(Form("dndpt2deCAT%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);	
 
    h = "dndthetadeCAT";    	
    dndthetadeCAT[ik] = new TH2D(Form("dndthetadeCAT%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);	 	 
	 
    h = "dndphideCAT";    	
    dndphideCAT[ik] = new TH2D(Form("dndphideCAT%d",ik), "", ndphi, -Phirange, Phirange, nde, 0.0, Erange);

    h = "dndetadeCAT";    	
    dndetadeCAT[ik] = new TH2D(Form("dndetadeCAT%d",ik), "", ndeta, Etarange, Etarange, nde, 0.0, Erange);

    h = "dndptdeRAD";    	
    dndptdeRAD[ik] = new TH2D(Form("dndptdeRAD%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);
	
    h = "dndpt2deRAD";    	
    dndpt2deRAD[ik] = new TH2D(Form("dndpt2deRAD%d",ik), "", ndpt, 0.0, Ptrange, nde, 0.0, Erange);	
 
    h = "dndthetadeRAD";    	
    dndthetadeRAD[ik] = new TH2D(Form("dndthetadeRAD%d",ik), "", ndtheta, 0.0, Thetarange, 10*nde, 0.0, Erange);	 	 
	 
    h = "dndphideRAD";    	
    dndphideRAD[ik] = new TH2D(Form("dndphideRAD%d",ik), "", ndphi, -Phirange, Phirange, nde, 0.0, Erange);

    h = "dndetadeRAD";    	
    dndetadeRAD[ik] = new TH2D(Form("dndetadeRAD%d",ik), "", ndeta, Etarange, Etarange, nde, 0.0, Erange);


	
	
    h = "dnde";    	
    dnde[ik] = new TH1D(Form("dnde%d",ik), "", nde, 0.0, Erange);
	
    h = "dndpt";    	
    dndpt[ik] = new TH1D(Form("dndpt%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2";    	
    dndpt2[ik] = new TH1D(Form("dndpt2%d",ik), "", ndpt, 0.0, Ptrange);	

    h = "dndpesudopt";     	
    dndpesudopt[ik] = new TH1D(Form("dndpesudopt%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2";    	
    dndpesudopt2[ik] = new TH1D(Form("dndpesudopt2%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndtheta";    	
    dndtheta[ik] = new TH1D(Form("dndtheta%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphi";    	
    dndphi[ik] = new TH1D(Form("dndphi%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndphi01";    	
    dndphi01[ik] = new TH1D(Form("dndphi01%d",ik), "", ndphi, -Phirange, Phirange);	
	
    h = "dndphi12";    	
    dndphi12[ik] = new TH1D(Form("dndphi12%d",ik), "", ndphi, -Phirange, Phirange);	

    h = "dndphi23";    	
    dndphi23[ik] = new TH1D(Form("dndphi23%d",ik), "", ndphi, -Phirange, Phirange);	
	
    h = "dndphi34";    	
    dndphi34[ik] = new TH1D(Form("dndphi34%d",ik), "", ndphi, -Phirange, Phirange);	
	
    h = "dndphi24";    	
    dndphi24[ik] = new TH1D(Form("dndphi24%d",ik), "", ndphi, -Phirange, Phirange);	
	
    h = "dndphi48";    	
    dndphi48[ik] = new TH1D(Form("dndphi48%d",ik), "", ndphi, -Phirange, Phirange);	

    h = "dndphi80";    	
    dndphi80[ik] = new TH1D(Form("dndphi80%d",ik), "", ndphi, -Phirange, Phirange);



	
    h = "dndeta";    	
    dndeta[ik] = new TH1D(Form("dndeta%d",ik), "", ndeta, -Etarange, Etarange);	
	

    h = "dndepos";    	
    dndepos[ik] = new TH1D(Form("dndepos%d",ik), "", nde, 0.0, Erange);
	
    h = "dndptpos";    	
    dndptpos[ik] = new TH1D(Form("dndptpos%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2pos";    	
    dndpt2pos[ik] = new TH1D(Form("dndpt2pos%d",ik), "", ndpt, 0.0, Ptrange);	

    h = "dndpesudoptpos";    	
    dndpesudoptpos[ik] = new TH1D(Form("dndpesudoptpos%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2pos";    	
    dndpesudopt2pos[ik] = new TH1D(Form("dndpesudopt2pos%d",ik), "", ndpt, 0.0, Ptrange);
	
    h = "dndthetapos";    	
    dndthetapos[ik] = new TH1D(Form("dndthetapos%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphipos";    	
    dndphipos[ik] = new TH1D(Form("dndphipos%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndetapos";    	
    dndetapos[ik] = new TH1D(Form("dndetapos%d",ik), "", ndeta, -Etarange, Etarange);

    h = "dndeneg";    	
    dndeneg[ik] = new TH1D(Form("dndeneg%d",ik), "", nde, 0.0, Erange);
	
    h = "dndptneg";   	
    dndptneg[ik] = new TH1D(Form("dndptneg%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2neg";    	
    dndpt2neg[ik] = new TH1D(Form("dndpt2neg%d",ik), "", ndpt, 0.0, Ptrange);	

    h = "dndpesudoptneg";   	
    dndpesudoptneg[ik] = new TH1D(Form("dndpesudoptneg%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2neg";    	
    dndpesudopt2neg[ik] = new TH1D(Form("dndpesudopt2neg%d",ik), "", ndpt, 0.0, Ptrange);
	
    h = "dndthetaneg";    	
    dndthetaneg[ik] = new TH1D(Form("dndthetaneg%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphineg";    	
    dndphineg[ik] = new TH1D(Form("dndphineg%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndetaneg";   	
    dndetaneg[ik] = new TH1D(Form("dndetaneg%d",ik), "", ndeta, -Etarange, Etarange);


    h = "dndeCAT";    	
    dndeCAT[ik] = new TH1D(Form("dndeCAT%d",ik), "", nde, 0.0, Erange);
	
    h = "dndptCAT";    	
    dndptCAT[ik] = new TH1D(Form("dndptCAT%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2CAT";    	
    dndpt2CAT[ik] = new TH1D(Form("dndpt2CAT%d",ik), "", ndpt, 0.0, Ptrange);	

    h = "dndpesudoptCAT";    	
    dndpesudoptCAT[ik] = new TH1D(Form("dndpesudoptCAT%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2CAT";    	
    dndpesudopt2CAT[ik] = new TH1D(Form("dndpesudopt2CAT%d",ik), "", ndpt, 0.0, Ptrange);
	
    h = "dndthetaCAT";    	
    dndthetaCAT[ik] = new TH1D(Form("dndthetaCAT%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphiCAT";    	
    dndphiCAT[ik] = new TH1D(Form("dndphiCAT%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndetaCAT";    	
    dndetaCAT[ik] = new TH1D(Form("dndetaCAT%d",ik), "", ndeta, -Etarange, Etarange);

    h = "dndeRAD";    	
    dndeRAD[ik] = new TH1D(Form("dndeRAD%d",ik), "", nde, 0.0, Erange);
	
    h = "dndptRAD";   	
    dndptRAD[ik] = new TH1D(Form("dndptRAD%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpt2RAD";    	
    dndpt2RAD[ik] = new TH1D(Form("dndpt2RAD%d",ik), "", ndpt, 0.0, Ptrange);	

    h = "dndpesudoptRAD";   	
    dndpesudoptRAD[ik] = new TH1D(Form("dndpesudoptRAD%d",ik), "", ndpt, 0.0, Ptrange);	
	
    h = "dndpesudopt2RAD";    	
    dndpesudopt2RAD[ik] = new TH1D(Form("dndpesudopt2RAD%d",ik), "", ndpt, 0.0, Ptrange);
	
    h = "dndthetaRAD";    	
    dndthetaRAD[ik] = new TH1D(Form("dndthetaRAD%d",ik), "", ndtheta, 0.0, Thetarange);	 	 
	 
    h = "dndphiRAD";    	
    dndphiRAD[ik] = new TH1D(Form("dndphiRAD%d",ik), "", ndphi, -Phirange, Phirange);

    h = "dndetaRAD";   	
    dndetaRAD[ik] = new TH1D(Form("dndetaRAD%d",ik), "", ndeta, -Etarange, Etarange);


	 	
    }	

//............................................................................................Histograms for LBT transport end.





	
	//}


//..some setting before the event loop	
    int    ndrphoton;
    int    KATTgamma;	 
	int    numnucleon=10000;
	double numjettotal;
	double randomxy;
	double R1;
	double XXX;
	double YYY;
	   
	double numjet[20000]={0.0};
	double Xnucleon[20000]={0.0};
	double Ynucleon[20000]={0.0};	

//......Geometry profile

	
    double etag, phig, Eg, Ej0;

    double eta;
    double phi1;
    double tanphi;
    double phia;
    double etagamma;
    double phi1gamma;
    double tanphigamma;
    double phiagamma;
    double energygamma;
    double ptgamma;
    double Egamma;

    double P0gamma;
    double P1gamma;
    double P2gamma;
    double P3gamma;
	
    int numGamma = 0;

    double Xjgamma=0.0;
    double PTjet=0.0;
    double PTgamma=0.0;
    double Rjgamma=0.0;
    double nRjgamma=0.0;

    double numgammawithjet=0.0;	
    double numgammajet=0.0;	
	
//......Dijet
               double leadingJx;
               double leadingJy;
               double leadingJz;
               double leadingJM;
               double leadingJE;

               double leadingJpt;		  
               double leadingJeta;		  
		 		  
               double leadingJtheta;
               double leadingJphi; 

               double leadingJphi2;

               double leadingtanphi;
               double leadingJphia;

               double subleadingJx;
               double subleadingJy;
               double subleadingJz;
               double subleadingJM;
               double subleadingJE;

               double subleadingJpt;		  
               double subleadingJeta;		  
		 		  
               double subleadingJtheta;
               double subleadingJphi; 

               double subleadingJphi2;

               double subleadingtanphi;
               double subleadingJphia;						
		
               double deltaphi01;	
	
//..some setting before the event loop end




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//..Tloop
//LBT->temp0medium

//..Eloop


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




	
//..event loop

	int numEvent=0;

    for(int n=1; ; ++n)
	{		
	
    if(numEvent>=LBT->ncall) break;

	numEvent=numEvent+1;
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////	


        LBT->LBTclear();  //clear particle list in LBT


//.............gammajet information

        ndrphoton=0;

	    P0gamma=0.0;
	    P1gamma=0.0;
	    P2gamma=0.0;
	    P3gamma=0.0;

//......read in initial Gamma 4-momentum
	
        int ne;
		int IDgamma;
		double asd1, asd2;		
		double PGm[4]={0.0};			
				
	            ndrphoton=ndrphoton+1;		

                PGm[0]=LBT->P[0][1];
                PGm[1]=-LBT->P[1][1];
                PGm[2]=-LBT->P[2][1];
                PGm[3]=-LBT->P[3][1];
				
                P0gamma=PGm[0];
                P1gamma=PGm[1];
                P2gamma=PGm[2];
                P3gamma=PGm[3];

                etagamma=1.0/2.0*log((PGm[0]+PGm[3])/(PGm[0]-PGm[3]));					
				
                tanphigamma=abs(P2gamma/P1gamma);
                phiagamma=atan(tanphigamma);

                if(P1gamma>0.0 && P2gamma>0.0)
		        {
	            phi1gamma=phiagamma;
		        }
                if(P1gamma<0.0 && P2gamma>0.0)
		        {
	            phi1gamma=3.1415926-phiagamma;
		        }
		        if(P1gamma<0.0 && P2gamma<0.0)
		        {
                phi1gamma=3.1415926+phiagamma;		
		        }
		        if(P1gamma>0.0 && P2gamma<0.0)
                {
                phi1gamma=2*3.1415926-phiagamma;		
		        }

                energygamma=P0gamma;
                ptgamma=sqrt(pow(P1gamma,2)+pow(P2gamma,2));					

//				cout<<"    ptgamma   " << ptgamma <<"   etagamma  "<< etagamma <<endl;

				double JE0=energygamma;

//......fastjet vector initialize
        int ui = -123456;
        int npart=100000;
	
        vector<double> px(npart), py(npart), pz(npart), E(npart);
        vector<int> pid(npart);		
//......fastjet vector initialize end
		
        numGamma=numGamma+1;			

        h100Ngamma->Fill(ptgamma,ww10);
		
		double epsilon=0.000000001;

//......time evolution in LBT

		double ti=LBT->time0;	//initial time inside or outside class?
		
		//LBT->nj=1;
		
		LBT->np=LBT->nj;
		
		cout<<"nj"<<" "<<LBT->nj<<endl;
		
        int ntimestep=floor(LBT->timend/LBT->dt);	
		
		if(LBT->LBTswitch==0)
		{
		ntimestep=1;
		}	

		if(LBT->Singlestepswitch==1)
		{
		ntimestep=1;
		ti=LBT->radlength;
		}

		if(LBT->Singlestepswitch==1)
		{        
        cout<< "ti" << " " << ti << " " << "radlength" << " " << LBT->radlength <<endl; 
        }
		
//......test		
		LBT->ntest22=0;
		LBT->ntestrad=0;
//......test

		
//......time loop		
		for(int timestep=1;timestep<=ntimestep;++timestep)
        {

		if(LBT->LBTswitch==1)
		{		
        ti=ti+LBT->dt;			
        LBT->LinearBoltzmannTransport(n,ti);
		
		if(LBT->Reachtauend==1)
		{
//		timestep=ntimestep;	
		}	
		
        }
		
//......histogram filled in		
		
    if(switchroot>0)
    {



//.............................................................................................histogram for single interaction


        for(int i=1;i<=LBT->nj;i++) // recoil and radiated partons
        {

            //cout<<"i"<<" "<<i<<" "<<"np"<<" "<<LBT->np<<endl;
			
            double E0=LBT->ener;
		
		    double px=LBT->P[1][i];
		    double py=LBT->P[2][i];		
		    double pz=LBT->P[3][i];

            double energy=LBT->P[0][i];
            double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));
			
            double rp = sqrt(pow(LBT->V[1][i],2) + pow(LBT->V[3][i],2)); //r=x**2+z**2
            double rz = LBT->V[2][i]; //y->z
            double rE = LBT->P[0][i];	

            double phip,phi2,pxphi;

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));	

            //cout<<"theta1"<<"---------------------  "<<thetap<<endl;			

            //py=-100.0;
			//px=1.0;
			//pz=0.0;            
			
			pxphi=px+epsilon;			
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
			double phia=atan(tanphi);

			//cout<<"tanphi"<<" "<<tanphi<<endl;
			
			//cout<<"phia"<<" "<<phia<<endl;
			
/*
            if (px>=0 && py>=0)
            {
            phi2=phia;
            }
            if (px<0 && py>=0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>=0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }		

            //phip=phi2-3.1415926/2.0;
            phip=phi2;			
			//cout<<"phip"<<" "<<phip<<endl;

			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////	
			
			//cout<<"partonCAT"<<" "<<LBT->CAT[i]<<endl;

			//cout<<"partonE"<<" "<<LBT->P[0][i]<<endl;
			
            if(LBT->CAT[i] == 0)
			{
            LEADdnde -> Fill(energy, wwdnde);
			LEADdndtheta -> Fill(thetap, wwdndtheta);			


			if(LBT->P[0][i]<=1.0)
			{			
			LEADdndtheta0010 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>1.0 && LBT->P[0][i]<=2.0)
			{			
			LEADdndtheta1020 -> Fill(thetap, wwdndtheta);			
			}
			if(LBT->P[0][i]>2.0 && LBT->P[0][i]<=4.0)
			{
			LEADdndtheta2040 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>4.0 && LBT->P[0][i]<=8.0)
			{
			LEADdndtheta4080 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>8.0)
			{	
			LEADdndtheta80 -> Fill(thetap, wwdndtheta);									
			}
			
			LEADdndpt -> Fill(pesudopt, wwdndpt);
			LEADdndpt2 -> Fill(pesudopt,  2.0*pesudopt*wwdndpt);
			LEADdndphi -> Fill(phip, wwdndphi);
			LEADdndeta -> Fill(etap, wwdndeta);
			LEADdndthetade -> Fill(thetap, energy, wwdndthetade);
			}
			

        	//numEvent=numEvent+1; //...........................................numEvent++
			
        } //for(int i=1;i<=np;i++) // recoil and radiated partons


        for(int i=LBT->nj+1;i<=LBT->np;i++) // recoil and radiated partons
        {

            //cout<<"i"<<" "<<i<<" "<<"np"<<" "<<LBT->np<<endl;
			
            double E0=LBT->ener;
		
		    double px=LBT->P[1][i];
		    double py=LBT->P[2][i];		
		    double pz=LBT->P[3][i];

            double energy=LBT->P[0][i];
            double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V[1][i],2) + pow(LBT->V[3][i],2)); //r=x**2+z**2
            double rz = LBT->V[2][i]; //y->z
            double rE = LBT->P[0][i];	

            double phip,phi2,pxphi;

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));						

            //cout<<"theta2"<<"---------------------  "<<thetap<<endl;

            //py=-100.0;
			//px=1.0;
			//pz=0.0;            
			
			pxphi=px+epsilon;			
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
			double phia=atan(tanphi);

			//cout<<"tanphi"<<" "<<tanphi<<endl;
			
			//cout<<"phia"<<" "<<phia<<endl;
			
/*
            if (px>=0 && py>=0)
            {
            phi2=phia;
            }
            if (px<0 && py>=0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>=0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }		

            //phip=phi2-3.1415926/2.0;
            phip=phi2;			
			//cout<<"phip"<<" "<<phip<<endl;

			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////	
			
			//cout<<"partonCAT"<<" "<<LBT->CAT[i]<<endl;

			//cout<<"partonE"<<" "<<LBT->P[0][i]<<endl;
			
            if(LBT->CAT[i] == 2)
			{
            CATdnde -> Fill(energy, wwdnde);
			CATdndtheta -> Fill(thetap, wwdndtheta);			


			if(LBT->P[0][i]<=1.0)
			{			
			CATdndtheta0010 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>1.0 && LBT->P[0][i]<=2.0)
			{			
			CATdndtheta1020 -> Fill(thetap, wwdndtheta);			
			}
			if(LBT->P[0][i]>2.0 && LBT->P[0][i]<=4.0)
			{
			CATdndtheta2040 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>4.0 && LBT->P[0][i]<=8.0)
			{
			CATdndtheta4080 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>8.0)
			{	
			CATdndtheta80 -> Fill(thetap, wwdndtheta);									
			}
			
			CATdndpt -> Fill(pesudopt, wwdndpt);
			CATdndpt2 -> Fill(pesudopt,  2.0*pesudopt*wwdndpt);
			CATdndphi -> Fill(phip, wwdndphi);
			CATdndeta -> Fill(etap, wwdndeta);
			CATdndthetade -> Fill(thetap, energy, wwdndthetade);
			}
			
            if(LBT->CAT[i] == 0)
			{			
            RADdnde -> Fill(energy, wwdnde);
			RADdndtheta -> Fill(thetap, wwdndtheta);
			
			
			if(LBT->P[0][i]<=1.0)
			{
			RADdndtheta0010 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>1.0 && LBT->P[0][i]<=2.0)
			{			
			RADdndtheta1020 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>2.0 && LBT->P[0][i]<=4.0)
			{
			RADdndtheta2040 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>4.0 && LBT->P[0][i]<=8.0)
			{
			RADdndtheta4080 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>8.0)
			{	
			RADdndtheta80 -> Fill(thetap, wwdndtheta);			
			}
			RADdndpt -> Fill(pesudopt, wwdndpt);
			RADdndpt2 -> Fill(pesudopt,  2.0*pesudopt*wwdndpt);
			RADdndphi -> Fill(phip, wwdndphi);
			RADdndeta -> Fill(etap, wwdndeta);
			RADdndthetade -> Fill(thetap, energy, wwdndthetade);

        	//numEvent=numEvent+1; //...........................................numEvent++
			
            }

        	//numEvent=numEvent+1; //...........................................numEvent++
			
        } //for(int i=1;i<=np;i++) // recoil and radiated partons


        for(int i=LBT->nj+1;i<=LBT->np;i++) // recoil and radiated partons
        {

    	    if(LBT->P0[0][i] == 0) continue;

            //cout<<"i"<<" "<<i<<" "<<"np"<<" "<<LBT->np<<endl;
			
            double E0=LBT->ener;
		
		    double px=LBT->P0[1][i];
		    double py=LBT->P0[2][i];		
		    double pz=LBT->P0[3][i];

            double energy=LBT->P0[0][i];
            double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V0[1][i],2) + pow(LBT->V0[3][i],2)); //r=x**2+z**2
            double rz = LBT->V0[2][i]; //y->z
            double rE = LBT->P0[0][i];	

            double phip,phi2,pxphi;

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));						

            //cout<<"theta3"<<"---------------------  "<<thetap<<endl;

            //py=-100.0;
			//px=1.0;
			//pz=0.0;            
			
			pxphi=px+epsilon;			
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
			double phia=atan(tanphi);

			//cout<<"tanphi"<<" "<<tanphi<<endl;
			
			//cout<<"phia"<<" "<<phia<<endl;
			
/*
            if (px>=0 && py>=0)
            {
            phi2=phia;
            }
            if (px<0 && py>=0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>=0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }		

            //phip=phi2-3.1415926/2.0;
            phip=phi2;			
			//cout<<"phip"<<" "<<phip<<endl;

			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////	
			
			//cout<<"partonCAT"<<" "<<LBT->CAT[i]<<endl;

			//cout<<"partonE"<<" "<<LBT->P[0][i]<<endl;
			
            if(LBT->CAT0[i] == 0)
			{
            NEGdnde -> Fill(energy, wwdnde);
			NEGdndtheta -> Fill(thetap, wwdndtheta);			


			if(LBT->P[0][i]<=1.0)
			{			
			NEGdndtheta0010 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>1.0 && LBT->P[0][i]<=2.0)
			{			
			NEGdndtheta1020 -> Fill(thetap, wwdndtheta);			
			}
			if(LBT->P[0][i]>2.0 && LBT->P[0][i]<=4.0)
			{
			NEGdndtheta2040 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>4.0 && LBT->P[0][i]<=8.0)
			{
			NEGdndtheta4080 -> Fill(thetap, wwdndtheta);
			}
			if(LBT->P[0][i]>8.0)
			{	
			NEGdndtheta80 -> Fill(thetap, wwdndtheta);									
			}
			
			NEGdndpt -> Fill(pesudopt, wwdndpt);
			NEGdndpt2 -> Fill(pesudopt,  2.0*pesudopt*wwdndpt);
			NEGdndphi -> Fill(phip, wwdndphi);
			NEGdndeta -> Fill(etap, wwdndeta);
			NEGdndthetade -> Fill(thetap, energy, wwdndthetade);
			}

        	//numEvent=numEvent+1; //...........................................numEvent++
			
        } //for(int i=LBT->nj+1;i<=LBT->np;i++) // negative partons

		
//.............................................................................................histogram for single interaction end



///////////////////////////////////////////////////////////////////////////////////////////////////////////////





//.............................................................................................histogram for multiple interactions (leading parton)
        for(int i=1;i<=LBT->nj;i++) //shower partons
        {
            double E0=LBT->ener;
		
		    double px=LBT->P[1][i];
		    double py=LBT->P[2][i];		
		    double pz=LBT->P[3][i];

            double energy=LBT->P[0][i];
            double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V[1][i],2) + pow(LBT->V[3][i],2)); //r=x**2+z**2
            double rz = LBT->V[2][i]; //y->z
            double rE = LBT->P[0][i];	

            double phip,phi2,pxphi;

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));						

            //py=-100.0;
			//px=1.0;
			//pz=0.0;            
			
			pxphi=px+epsilon;			
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
			double phia=atan(tanphi);

			//cout<<"tanphi"<<" "<<tanphi<<endl;
			
			//cout<<"phia"<<" "<<phia<<endl;
			
/*
            if (px>=0 && py>=0)
            {
            phi2=phia;
            }
            if (px<0 && py>=0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>=0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }		

            //phip=phi2-3.1415926/2.0;
            phip=phi2;			
			//cout<<"phip"<<" "<<phip<<endl;

			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////	
			
			Sleadingpartonenergy->Fill(ti, energy*ww10);
			Sleadingpartonpt->Fill(ti, pt*ww10);			
			Sleadingelossf->Fill(ti, (E0-energy)/E0*ww10);
			Sleadingeloss->Fill(ti, (E0-energy)*ww10);
			Sleadingtheta->Fill(ti, thetap*ww10);
			Sleadingphi->Fill(ti, phip*ww10);
			Sleadingeta->Fill(ti, etap*ww10);	
            Sleadingpartonpesudopt->Fill(ti, pesudopt*ww10);
            Sleadingpartonpesudopt2->Fill(ti, pesudopt*pesudopt*ww10);			
            Sleadingpartonqhat->Fill(ti, pesudopt*pesudopt/ti*ww10);			
        } //for(int i=1;i<=nj;i++) //shower partons


//.............................................................................................histogram for multiple/single interactions (leading parton)




///////////////////////////////////////////////////////////////////////////////////////////////////////////////






//.............................................................................................histogram for multiple interactions

        for (unsigned ik = 1; ik <= ntistep; ++ik) // time step of the final histogram
        {			

        if (fabs(ti - ik * dtstep) < 0.000001) // check time
        {
		
        if(LBT->P[0][1]<2.0)
        {
		//exit(1);		
		}
		

	    //leading parton
        //for(int i=1;i<=1;i++) //shower partons		
        for(int i=1;i<=LBT->nj;i++) //.........................................................................shower partons
        {
            double E0=LBT->ener;
		
		    double px=LBT->P[1][i];
		    double py=LBT->P[2][i];		
		    double pz=LBT->P[3][i];

            double energy=LBT->P[0][i];
            double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V[1][i],2) + pow(LBT->V[3][i],2)); //r=x**2+z**2
            double rz = LBT->V[2][i]; //y->z
            double rE = LBT->P[0][i];	

            double phip,phi2,pxphi;

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));						

            //py=-100.0;
			//px=1.0;
			//pz=0.0;
			
			pxphi=px+epsilon;			
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
			double phia=atan(tanphi);

			//cout<<"tanphi"<<" "<<tanphi<<endl;
			
			//cout<<"phia"<<" "<<phia<<endl;
			
/*
            if (px>=0 && py>=0)
            {
            phi2=phia;
            }
            if (px<0 && py>=0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>=0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }		

            //phip=phi2-3.1415926/2.0;
            phip=phi2;			
			//cout<<"phip"<<" "<<phip<<endl;

			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////	
			
			leadingpartonenergy->Fill(ti-epsilon, energy*ww10);
			leadingpartonpt->Fill(ti-epsilon, pesudopt*ww10);			
			leadingelossf->Fill(ti-epsilon, (E0-energy)/E0*ww10);
			leadingeloss->Fill(ti-epsilon, (E0-energy)*ww10);
			leadingtheta->Fill(ti-epsilon, thetap*ww10);
			leadingphi->Fill(ti-epsilon, phip*ww10);
			leadingeta->Fill(ti-epsilon, etap*ww10);			


            dndedtleading0->Fill(energy, ti-epsilon, wwdnde);
            dndptdtleading0->Fill(pesudopt, ti-epsilon, wwdndpt);
            dndpt2dtleading0->Fill(pesudopt, ti-epsilon, 2.0*pesudopt*wwdndpt);			
			dndthetadtleading0->Fill(thetap, ti-epsilon, wwdndtheta);
			dndphidtleading0->Fill(phip, ti-epsilon, wwdndphi);	
			dndetadtleading0->Fill(etap, ti-epsilon, wwdndeta);

			
            dndeleading[ik]->Fill(energy, wwdnde);
            dndptleading[ik]->Fill(pt, wwdndpt);
            dndpt2leading[ik]->Fill(pt, 2.0*pt*wwdndpt);		

            dndpesudoptleading[ik]->Fill(pesudopt, wwdndpt);
            dndpesudopt2leading[ik]->Fill(pesudopt, 2.0*pesudopt*wwdndpt);
			
			dndthetaleading[ik]->Fill(thetap, wwdndtheta);
			dndphileading[ik]->Fill(phip, wwdndphi);	
			dndetaleading[ik]->Fill(etap, wwdndeta);
			
        } //for(int i=1;i<=nj;i++) //shower partons							
			
//////////////////////////////////////////...check process..........................................?
        double ip=0;
        for(int i=1;i<=LBT->np;i++)
        {
            if(LBT->P[0][i] == 0) ip+=1;
        }							
        double nnn=LBT->np-ip;	
//////////////////////////////////////////	

        //for(int i=1;i<=1;i++)
        for(int i=LBT->nj+1;i<=LBT->np;i++) //positive partons       //no leading parton
        {
    	    if(LBT->P[0][i] == 0) continue;

            double E0=LBT->ener;
		
		    double px=LBT->P[1][i];
		    double py=LBT->P[2][i];		
		    double pz=LBT->P[3][i];

            double energy=LBT->P[0][i];
            //double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));
			
			double pt=pesudopt;
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V[1][i],2) + pow(LBT->V[3][i],2)); //r=x**2+z**2
            double rz = LBT->V[2][i]; //y->z
            double rx = LBT->V[1][i]; //y->z			
            double rE = LBT->P[0][i];	

            double phip,phi2,pxphi;			

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));			

			pxphi=px-epsilon;
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
            double phia=atan(tanphi);

/*
            if (px>0 && py>0)
            {
            phi2=phia;
            }
            if (px<0 && py>0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/

	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }				
			
			
			
			
			
			
            //phip=phi2-3.1415926/2.0;
            phip=phi2;	



/*
			double tanphi = P[1][i]/P[2][i];			
            double phia = atan(tanphi); 
            double phip;		

            double thetap;			//?????????????
			
            if(P[2][i]>=0.0)
		    {
	        phip=phia;
		    }
		    if(P[2][i]<0.0 && P[1][i]>=0.0)
		    {
            phip=3.1415926+phia;		
		    }
		    if(P[2][i]<0.0 && P[1][i]<0.0)
            {
            phip=-3.1415926+phia;
		    }
*/














			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////				
		
            rdedrdz[ik]->Fill(rp, rz, rp * rE * wwdndrdz * 2.0);
            dedrdz[ik]->Fill(rp, rz, rE * wwdndrdz * 2.0);			
            dndrdz[ik]->Fill(rp, rz, wwdndrdz * 2.0);			
            rdndrdz[ik]->Fill(rp, rz, rp * wwdndrdz * 2.0);
			
            rdedrdzpos[ik]->Fill(rp, rz, rp * rE * wwdndrdz * 2.0);
            dedrdzpos[ik]->Fill(rp, rz, rE * wwdndrdz * 2.0);			
            dndrdzpos[ik]->Fill(rp, rz, wwdndrdz * 2.0);			
            rdndrdzpos[ik]->Fill(rp, rz, rp * wwdndrdz * 2.0);			

            if(LBT->CAT[i] == 2)
			{
            rdedrdzCAT[ik]->Fill(rp, rz, rp * rE * wwdndrdz * 2.0);
            dedrdzCAT[ik]->Fill(rp, rz, rE * wwdndrdz * 2.0);			
            dndrdzCAT[ik]->Fill(rp, rz, wwdndrdz * 2.0);			
            rdndrdzCAT[ik]->Fill(rp, rz, rp * wwdndrdz * 2.0);
            }
			
            if(LBT->CAT[i] == 0)
			{
            rdedrdzRAD[ik]->Fill(rp, rz, rp * rE * wwdndrdz * 2.0);
            dedrdzRAD[ik]->Fill(rp, rz, rE * wwdndrdz * 2.0);			
            dndrdzRAD[ik]->Fill(rp, rz, wwdndrdz * 2.0);			
            rdndrdzRAD[ik]->Fill(rp, rz, rp * wwdndrdz * 2.0);
			}

			
/////////////////////////////////////////////////////////////////			
            rdedxdz[ik]->Fill(rx, rz, rp * rE * wwdndrdz);
            dedxdz[ik]->Fill(rx, rz, rE * wwdndrdz);			
            dndxdz[ik]->Fill(rx, rz, wwdndrdz);			
            rdndxdz[ik]->Fill(rx, rz, rp * wwdndrdz);
			
            rdedxdzpos[ik]->Fill(rx, rz, rp * rE * wwdndrdz);
            dedxdzpos[ik]->Fill(rx, rz, rE * wwdndrdz);			
            dndxdzpos[ik]->Fill(rx, rz, wwdndrdz);			
            rdndxdzpos[ik]->Fill(rx, rz, rp * wwdndrdz);

            if(LBT->CAT[i] == 2)
			{
            rdedxdzCAT[ik]->Fill(rx, rz, rp * rE * wwdndrdz);
            dedxdzCAT[ik]->Fill(rx, rz, rE * wwdndrdz);			
            dndxdzCAT[ik]->Fill(rx, rz, wwdndrdz);			
            rdndxdzCAT[ik]->Fill(rx, rz, rp * wwdndrdz);
			}
			
            if(LBT->CAT[i] == 0)
			{
            rdedxdzRAD[ik]->Fill(rx, rz, rp * rE * wwdndrdz);
            dedxdzRAD[ik]->Fill(rx, rz, rE * wwdndrdz);			
            dndxdzRAD[ik]->Fill(rx, rz, wwdndrdz);			
            rdndxdzRAD[ik]->Fill(rx, rz, rp * wwdndrdz);
			}
			
/////////////////////////////////////////////////////////////////
			
			dedetadphi[ik]->Fill(phip, etap, energy * wwdndetadphi);
			dndetadphi[ik]->Fill(phip, etap, wwdndetadphi);			

			dedetadphipos[ik]->Fill(phip, etap, energy * wwdndetadphi);
			dndetadphipos[ik]->Fill(phip, etap, wwdndetadphi);

            if(LBT->CAT[i] == 2)
			{
			dedetadphiCAT[ik]->Fill(phip, etap, energy * wwdndetadphi);
			dndetadphiCAT[ik]->Fill(phip, etap, wwdndetadphi);
            }

            if(LBT->CAT[i] == 0)
			{
			dedetadphiRAD[ik]->Fill(phip, etap, energy * wwdndetadphi);
			dndetadphiRAD[ik]->Fill(phip, etap, wwdndetadphi);
			}

            dndptde[ik]->Fill(pesudopt, energy, wwdndptde);
            dndpt2de[ik]->Fill(pesudopt, energy, 2.0*pesudopt*wwdndptde);			
			dndthetade[ik]->Fill(thetap, energy, wwdndthetade);
			dndphide[ik]->Fill(phip, energy, wwdndphide);	
			dndetade[ik]->Fill(etap, energy, wwdndetade);			

            dndptdepos[ik]->Fill(pesudopt, energy, wwdndptde);
            dndpt2depos[ik]->Fill(pesudopt, energy, 2.0*pesudopt*wwdndptde);			
			dndthetadepos[ik]->Fill(thetap, energy, wwdndthetade);
			dndphidepos[ik]->Fill(phip, energy, wwdndphide);	
			dndetadepos[ik]->Fill(etap, energy, wwdndetade);			


            if(LBT->CAT[i] == 2)
			{			
            dndptdeCAT[ik]->Fill(pesudopt, energy, wwdndptde);
            dndpt2deCAT[ik]->Fill(pesudopt, energy, 2.0*pesudopt*wwdndptde);			
			dndthetadeCAT[ik]->Fill(thetap, energy, wwdndthetade);
			dndphideCAT[ik]->Fill(phip, energy, wwdndphide);	
			dndetadeCAT[ik]->Fill(etap, energy, wwdndetade);
            }
			
            if(LBT->CAT[i] == 0)
			{			
            dndptdeRAD[ik]->Fill(pesudopt, energy, wwdndptde);
            dndpt2deRAD[ik]->Fill(pesudopt, energy, 2.0*pesudopt*wwdndptde);			
			dndthetadeRAD[ik]->Fill(thetap, energy, wwdndthetade);
			dndphideRAD[ik]->Fill(phip, energy, wwdndphide);	
			dndetadeRAD[ik]->Fill(etap, energy, wwdndetade);
            }
			
            dnde[ik]->Fill(energy, wwdnde);
            dndpt[ik]->Fill(pesudopt, wwdndpt);
            dndpt2[ik]->Fill(pesudopt, 2.0*pt*wwdndpt);

            //dndpesudopt[ik]->Fill(pesudopt, wwdndpt);
            //dndpesudopt2[ik]->Fill(pesudopt, 2.0*pt*wwdndpt);
			
			dndtheta[ik]->Fill(thetap, wwdndtheta);
			dndphi[ik]->Fill(phip, wwdndphi);

            if(pt>0.0 && pt<=1.0)
			{
			dndphi01[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>1.0 && pt<=2.0)
			{
			dndphi12[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>2.0 && pt<=3.0)
			{
			dndphi23[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>3.0 && pt<=4.0)
			{
			dndphi34[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>2.0 && pt<=4.0)
			{
			dndphi24[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>4.0 && pt<=8.0)
			{
			dndphi48[ik]->Fill(phip, wwdndphi);			
			}

            if(pt>8.0)
			{
			dndphi80[ik]->Fill(phip, wwdndphi);			
			}
			
			dndeta[ik]->Fill(etap, wwdndeta);			

            dndepos[ik]->Fill(energy, wwdnde);
            dndptpos[ik]->Fill(pt, wwdndpt);
            dndpt2pos[ik]->Fill(pt, 2.0*pt*wwdndpt);
			dndthetapos[ik]->Fill(thetap, wwdndtheta);
			dndphipos[ik]->Fill(phip, wwdndphi);
			dndetapos[ik]->Fill(etap, wwdndeta);			

            if(LBT->CAT[i] == 2)
			{
            dndeCAT[ik]->Fill(energy, wwdnde);
            dndptCAT[ik]->Fill(pt, wwdndpt);
            dndpt2CAT[ik]->Fill(pt, 2.0*pt*wwdndpt);
			dndthetaCAT[ik]->Fill(thetap, wwdndtheta);
			dndphiCAT[ik]->Fill(phip, wwdndphi);
			dndetaCAT[ik]->Fill(etap, wwdndeta);
            }

            if(LBT->CAT[i] == 0)
			{			
            dndeRAD[ik]->Fill(energy, wwdnde);
            dndptRAD[ik]->Fill(pt, wwdndpt);
            dndpt2RAD[ik]->Fill(pt, 2.0*pt*wwdndpt);
			dndthetaRAD[ik]->Fill(thetap, wwdndtheta);
			dndphiRAD[ik]->Fill(phip, wwdndphi);
			dndetaRAD[ik]->Fill(etap, wwdndeta);
			}
			
        } //for(int i=2;i<=np;i++) //positive partons

//////////////////////////////////////////...check process		
        double ip0=0;
        for(int i=1;i<=LBT->np;i++)
        {
            if(LBT->P0[0][i] == 0) ip0+=1;
        }	
        double nnn0=LBT->np-ip0;
//////////////////////////////////////////

        //for(int i=0;i<=0;i++) //negative partons
		
        for(int i=LBT->nj+1;i<=LBT->np;i++) //negative partons
        {
    	    if(LBT->P0[0][i] == 0) continue;
			
			
            double E0=LBT->ener;
		
		    double px=LBT->P0[1][i];
		    double py=LBT->P0[2][i];		
		    double pz=LBT->P0[3][i];

            double energy=LBT->P0[0][i];
            //double pt=sqrt(pow(px,2)+pow(py,2));				
            double pesudopt=sqrt(pow(px,2)+pow(pz,2));

			double pt=pesudopt;
			
			double etap = 1.0/2.0*log((energy+pz)/(energy-pz));			
			
            double rp = sqrt(pow(LBT->V0[1][i],2) + pow(LBT->V0[3][i],2)); //r=x**2+z**2
            double rz = LBT->V0[2][i]; //y->z
            double rx = LBT->V0[1][i]; //y->z			
            double rE = LBT->P0[0][i];	

            double phip,phi2,pxphi;			

            double thetap=acos(py/sqrt(pow(px,2)+pow(py,2)+pow(pz,2)));						

			pxphi=px-epsilon;
			
            //double tanphi=abs(py/pxphi);
            double tanphi=(px/py);			
            double phia=atan(tanphi);

/*
            if (px>0 && py>0)
            {
            phi2=phia;
            }
            if (px<0 && py>0)
            {
            phi2=3.14159-phia;
            }
            if (px<0 && py<0) 
            {
            phi2=3.14159+phia;
            }
            if (px>0 && py<0)
            {
            phi2=2.0*3.14159-phia;
            }
*/
			
	        if(py>=0.0)
		    {
	        phi2=phia;
		    }
		    if(py<0.0 && px>=0.0)
		    {
            phi2=3.1415926+phia;		
		    }
		    if(py<0.0 && px<0.0)
            {
            phi2=-3.1415926+phia;
		    }			
	
            //phip=phi2-3.1415926/2.0;
            phip=phi2;
			
///////////////////////////////////////////////////
/*			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }
*/			
///////////////////////////////////////////////////					
			
            rdedrdz[ik]->Fill(rp, rz, -rp * energy * wwdndrdz * 2.0);							
            dedrdz[ik]->Fill(rp, rz, -energy * wwdndrdz * 2.0);
            dndrdz[ik]->Fill(rp, rz, -wwdndrdz * 2.0);
            rdndrdz[ik]->Fill(rp, rz, -rp * wwdndrdz * 2.0);
			
            rdedrdzneg[ik]->Fill(rp, rz, rp * rE * wwdndrdz * 2.0);
            dedrdzneg[ik]->Fill(rp, rz, rE * wwdndrdz * 2.0);
            dndrdzneg[ik]->Fill(rp, rz, wwdndrdz * 2.0);
            rdndrdzneg[ik]->Fill(rp, rz, rp * wwdndrdz * 2.0);


/////////////////////////////////////////////////////////////////			
            rdedxdz[ik]->Fill(rx, rz, -rp * rE * wwdndrdz);
            dedxdz[ik]->Fill(rx, rz, -rE * wwdndrdz);			
            dndxdz[ik]->Fill(rx, rz, -wwdndrdz);			
            rdndxdz[ik]->Fill(rx, rz, -rp * wwdndrdz);
			
            rdedxdzneg[ik]->Fill(rx, rz, rp * rE * wwdndrdz);
            dedxdzneg[ik]->Fill(rx, rz, rE * wwdndrdz);			
            dndxdzneg[ik]->Fill(rx, rz, wwdndrdz);			
            rdndxdzneg[ik]->Fill(rx, rz, rp * wwdndrdz);	
/////////////////////////////////////////////////////////////////

			
			dedetadphi[ik]->Fill(phip, etap, -energy * wwdndetadphi);
			dndetadphi[ik]->Fill(phip, etap, -wwdndetadphi);

			dedetadphineg[ik]->Fill(phip, etap, energy * wwdndetadphi);
			dndetadphineg[ik]->Fill(phip, etap, wwdndetadphi);


            dndptde[ik]->Fill(pt, energy, -wwdndptde);
            dndpt2de[ik]->Fill(pt, energy, -2.0*pt*wwdndptde);		
			dndthetade[ik]->Fill(thetap, energy, -wwdndthetade);
			dndphide[ik]->Fill(phip, energy, -wwdndphide);
			dndetade[ik]->Fill(etap, energy, -wwdndetade);

            dndptdeneg[ik]->Fill(pt, energy, wwdndptde);
            dndpt2deneg[ik]->Fill(pt, energy, 2.0*pt*wwdndptde);		
			dndthetadeneg[ik]->Fill(thetap, energy, wwdndthetade);
			dndphideneg[ik]->Fill(phip, energy, wwdndphide);
			dndetadeneg[ik]->Fill(etap, energy, wwdndetade);	
			
			
            dnde[ik]->Fill(energy, -wwdnde);
            dndpt[ik]->Fill(pesudopt, -wwdndpt);
            dndpt2[ik]->Fill(pesudopt, -2.0*pesudopt*wwdndpt);
			dndtheta[ik]->Fill(thetap, -wwdndtheta);
			dndphi[ik]->Fill(phip, -wwdndphi);
			

            if(pt>0.0 && pt<=1.0)
			{
			dndphi01[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>1.0 && pt<=2.0)
			{
			dndphi12[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>2.0 && pt<=3.0)
			{
			dndphi23[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>3.0 && pt<=4.0)
			{
			dndphi34[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>2.0 && pt<=4.0)
			{
			dndphi24[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>4.0 && pt<=8.0)
			{
			dndphi48[ik]->Fill(phip, -wwdndphi);			
			}

            if(pt>8.0)
			{
			dndphi80[ik]->Fill(phip, -wwdndphi);			
			}
		
			
			dndeta[ik]->Fill(etap, -wwdndeta);

            dndeneg[ik]->Fill(energy, wwdnde);
            dndptneg[ik]->Fill(pt, wwdndpt);
            dndpt2neg[ik]->Fill(pt, 2.0*pt*wwdndpt);
			dndthetaneg[ik]->Fill(thetap, wwdndtheta);
			dndphineg[ik]->Fill(phip, wwdndphi);
			dndetaneg[ik]->Fill(etap, wwdndeta);
			
        } // for(int i=2;i<=np;i++) //negative partons

        } // if (fabs(ti - ik * dtstep) < 0.000001) // check time
	
        } // for (unsigned ik = 1; ik <= ntistep; ++ik)				
		
	    }//if(switchroot>0)	

//......histogram filled end		

            int switchoutput=1;
            if(switchoutput>0)
            {		
		
            if(timestep==ntimestep)				
            {
				
                int ip=0;
                for(int i=1;i<=LBT->np;i++)
                {
                    if(LBT->P[0][i] == 0) ip+=1;
                }	
                int nnn=LBT->np-ip;
                numpositiveLBT<<n<<" "<<LBT->np<<endl;
                //numpositiveLBT<<n<<" "<<nnn<<endl;				
                for(int i=1;i<=LBT->np;i++)
                {
                    //if(LBT->P[0][i] == 0) continue;
                    positiveLBT<<n<<" "<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i] <<" "<<LBT->Vfrozen[1][i]<<" "<<LBT->Vfrozen[2][i]<<" "<<LBT->Vfrozen[3][i]<<" "<<LBT->Vfrozen[0][i]<<" "<<LBT->CAT[i]<<endl;
                }	
                    
                int ip0=0;
                for(int i=LBT->nj+1;i<=LBT->np;i++)
                {
                    if(LBT->P0[0][i] == 0) ip0+=1;
                }	
                int nnn0=LBT->np-ip0-LBT->nj;
                numnegativeLBT<<n<<" "<<LBT->np<<endl;
                //numnegativeLBT<<n<<" "<<nnn0<<endl;				
                for(int i=LBT->nj+1;i<=LBT->np;i++)
                {
                    //if(LBT->P0[0][i] == 0) continue;
                    negativeLBT<<n<<" "<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<LBT->Vfrozen0[1][i]<<" "<<LBT->Vfrozen0[2][i]<<" "<<LBT->Vfrozen0[3][i]<<" "<<LBT->Vfrozen0[0][i]<<" "<<LBT->CAT0[i]<<endl;
                }
            }//if(timestep==ntimestep)		
		
		    }//if(switchoutput>0)
		
		
            for (unsigned ik = 1; ik <= ntistep; ++ik)
            {		

            if (fabs(ti - ik * dtstep) < 0.000001) // check time
            {			
			
            int switchFastjet=1;
            if(switchFastjet>0)
            {		

                //if(timestep==ntimestep)				
                {
		        int njfast=-1;

                vector<fastjet::PseudoJet> input_particles;	    //@@@???

                for(int i=1;i<=LBT->np;i++) //positive
                {

                if(LBT->P[0][i] == 0) continue;
			
//		        njfast=i-1;

			    njfast=njfast+1;
		
	            E[njfast] = LBT->P[0][i];
	            px[njfast] = LBT->P[1][i];
	            py[njfast] = LBT->P[2][i];
	            pz[njfast] = LBT->P[3][i];

///////////////////////////////////////////////...test		
//              cout <<  px[k0] << " " << py[k0] << " " << pz[k0] << " " << E[k0] << " " << pid[k0] << endl;
//              cout <<  px[njfast] << " " << py[njfast] << " " << pz[njfast] << " " << E[njfast] << " " << endl;
///////////////////////////////////////////////...test

                fastjet::PseudoJet pj;		
                if ( E[njfast] < 0 ) 
                { 
                pj.reset_momentum(-px[njfast],-py[njfast],-pz[njfast],-E[njfast]);
                pj.set_user_index(-njfast); 
                }
                else 
                {
                pj.reset_momentum(px[njfast],py[njfast],pz[njfast],E[njfast]);
                pj.set_user_index(njfast);
                }

                input_particles.push_back(pj);			
		
                }//for(int i=1;i<=LBT->np;i++) //positive
						
                for(int i=LBT->nj+1;i<=LBT->np;i++) //negative
                {
			
                if(LBT->P0[0][i] == 0 || LBT->P0[0][i] == 10000) continue;			

		        njfast=njfast+1;
		
	            E[njfast] = -LBT->P0[0][i];
	            px[njfast] = -LBT->P0[1][i];
	            py[njfast] = -LBT->P0[2][i];
	            pz[njfast] = -LBT->P0[3][i];
			
///////////////////////////////////////////////		
//              cout <<  px[k0] << " " << py[k0] << " " << pz[k0] << " " << E[k0] << " " << pid[k0] << endl;
//              cout <<  px[njfast] << " " << py[njfast] << " " << pz[njfast] << " " << E[njfast] << " " << endl;
///////////////////////////////////////////////

                fastjet::PseudoJet pj;
                if ( E[njfast] < 0 ) 
                { 
                pj.reset_momentum(-px[njfast],-py[njfast],-pz[njfast],-E[njfast]);
                pj.set_user_index(-njfast); 
                }
                else 
                {
                pj.reset_momentum(px[njfast],py[njfast],pz[njfast],E[njfast]);
                pj.set_user_index(njfast);
                }

                //input_particles.push_back(pj);			
			
                }//for(int i=nj+1;i<=LBT->np;i++) //negative		


//..................................fastjet loading 


                // create a jet definition: 
                // a jet algorithm with a given radius parameter	  
//	            cout<<"R "<<R<<endl;
	  
                fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
                //fastjet::JetDefinition jet_def(fastjet::kt_algorithm, R);	  
                // create an instance of the negative energy recombiner, with a given flag ui
                NegativeEnergyRecombiner uir(ui);
                // tell jet_def to use this new recombiner
                jet_def.set_recombiner(&uir);

                // run the jet clustering with the above jet definition
                //----------------------------------------------------------
                fastjet::ClusterSequence clust_seq(input_particles, jet_def);

                // get the resulting jets ordered in pt
                //----------------------------------------------------------
                double ptmin = 0.0;
                vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
				
//...................................................................................................................jet grooming				
                //ClusterSequence cs(event, jet_def);
                vector<PseudoJet> NO1jet = SelectorNHardest(1)(clust_seq.inclusive_jets());				
//...................................................................................................................jet grooming					
				
                //if (j == 0)Ej0 = inclusive_jets[0].E(); // Initial jet energy
                Ej0 = energygamma; // Initial jet energy
                //cout << "Ej0=" << Ej0;
                //if (Ej0 ==0) cout << "jet energy 0! "<< i << " " << j << " " <<inclusive_jets.size()<<endl;
                //cout << Ej0 << endl;

                // tell the user what was done
                //  - the description of the algorithm used
                //  - extract the inclusive jets with pt > 5 GeV
                //    show the output as
                //      {index, rap, phi, pt}
                //----------------------------------------------------------
                // cout << "Ran " << jet_def.description() << endl;

//..................................fastjet loading end


//..................................fastjet loading

//...................................................................................................... fastjet Jet!



//......jet finding


        numgammajet=0.0;

        //for (unsigned int k = 0; k < inclusive_jets.size(); k++)		
        for (unsigned int k = 0; k < 1; k++) 
        {

//......jet with positive energy  @@@???  
	  
        if (inclusive_jets[k].user_index()>=0) 
        {

          double Jx=inclusive_jets[k].px();
          double Jy=inclusive_jets[k].py();
          double Jz=inclusive_jets[k].pz();
          double JM=sqrt(pow(Jx,2)+pow(Jy,2)+pow(Jz,2));
          double JE=inclusive_jets[k].E();

		  double Jpesudopt=sqrt(pow(Jx,2)+pow(Jz,2));
//...917???		  
		  
          Jpt=inclusive_jets[k].pt();
		  
          //Jeta=1.0/2.0*log((JE+Jz)/(JE-Jz));
          Jeta=inclusive_jets[k].eta();

		  
          double Jtheta=acos(inclusive_jets[k].py()/sqrt(pow(inclusive_jets[k].px(),2)+pow(inclusive_jets[k].py(),2)+pow(inclusive_jets[k].pz(),2)));

          double Jphi = inclusive_jets[k].phi_std()-phig;

///////////////////////////////////////////////
//     cout <<  "Jx" << " " << "Jy" << " " << "Jz" << " " << "JE" << " " << "Jpt" << " " << "Jeta" << " " << "Jtheta" << " " << "Jphi" << " " << endl;	
//     cout <<   Jx  << " " <<  Jy  << " " <<  Jz  << " " <<  JE  << " " <<  Jpt  << " " <<  Jeta  << " " <<  Jtheta  << " " <<  Jphi  << " " << endl;
///////////////////////////////////////////////

            if (Jphi<0) Jphi += 2*3.14159;
            if (Jphi>3.14159) Jphi = 2*3.14159-Jphi;

            double Jphi2;
			
//////////////////////////////////////

            double tanphi=abs(Jy/Jx);
            double Jphia=atan(tanphi);

            if (Jx>0 && Jy>0)
            {
            Jphi2=Jphia;
            }
            if (Jx<0 && Jy>0)
            {
            Jphi2=3.14159-Jphia;
            }
            if (Jx<0 && Jy<0) 
            {
            Jphi2=3.14159+Jphia;
            }
            if (Jx>0 && Jy<0) 
            {
            Jphi2=2.0*3.14159-Jphia;
            }

//////////////////////////////////////
			
            deltaphi=Jphi2-phi1gamma;			
			
            if(deltaphi<0)
            {
            deltaphi=deltaphi+2.0*3.1415926;
            }

            if(deltaphi>3.1415926)
            {
            deltaphi=2.0*3.1415926-deltaphi;
            }

            leadingJetenergy->Fill(ti-epsilon, JE*wwJ10);
            leadingJetpt->Fill(ti-epsilon, Jpt*wwJ10);
            leadingJetpesudopt->Fill(ti-epsilon, Jpesudopt*wwJ10);			
            leadingJetelossf->Fill(ti-epsilon, (JE0-JE)/JE0*wwJ10);
            leadingJeteloss->Fill(ti-epsilon, (JE0-JE)*wwJ10);
            leadingJettheta->Fill(ti-epsilon, Jtheta*wwJ10);
            leadingJetphi->Fill(ti-epsilon, deltaphi*wwJ10);
            leadingJeteta->Fill(ti-epsilon, Jeta*wwJ10);

            SleadingJetenergy->Fill(ti, JE*wwJ10);
            SleadingJetpt->Fill(ti, Jpt*wwJ10);
            SleadingJetpesudopt->Fill(ti, Jpesudopt*wwJ10);			
            SleadingJetelossf->Fill(ti, (JE0-JE)/JE0*wwJ10);
            SleadingJeteloss->Fill(ti, (JE0-JE)*wwJ10);
            SleadingJettheta->Fill(ti, Jtheta*wwJ10);
            SleadingJetphi->Fill(ti, deltaphi*wwJ10);
            SleadingJeteta->Fill(ti, Jeta*wwJ10);
			
///////////////////////////////////////////////	
//     cout <<  "Jx" << " " << "Jy" << " " << "Jz" << " " << "JE" << " " << "Jpt" << " " << "Jeta" << " " << "Jtheta" << " " << "deltaphi" << " " << endl;	
//     cout <<   Jx  << " " <<  Jy  << " " <<  Jz  << " " <<  JE  << " " <<  Jpt  << " " <<  Jeta  << " " <<  Jtheta  << " " <<  deltaphi  << " " << endl;
///////////////////////////////////////////////

//.........cms cut
//         if(deltaphi>Jdeltaphicut && abs(Jeta)<Jetacut && Jpt>Jptcutmin && Jpt<Jptcutmax) 
		   {

///////////////////////////////////////////////	
//           cout <<  "Jx" << " " << "Jy" << " " << "Jz" << " " << "JE" << " " << "Jpt" << " " << "Jeta" << " " << "Jtheta" << " " << "deltaphi" << " " << endl;	
//           cout <<   Jx  << " " <<  Jy  << " " <<  Jz  << " " <<  JE  << " " <<  Jpt  << " " <<  Jeta  << " " <<  Jtheta  << " " <<  deltaphi  << " " << endl;
///////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming




  // declare a SoftDropTagger to play with
  double beta = 0.0;
  double zcut = 0.1;
  double mu   = 1.0;
  contrib::SoftDropTagger soft_drop(beta, zcut, mu);

//cout<<"LBT GROOMING--------------------------------------------------start"<<endl;  
  
  // apply the soft drop to the jets
  vector<PseudoJet> sd_jets = soft_drop(NO1jet);

//cout<<"LBT GROOMING--------------------------------------------------end"<<endl;   
  
//.........................................................................................loop

    // create handful shortcuts for te original and soft-dropped jets
    const PseudoJet &jetorigin    = NO1jet[k];
    const PseudoJet &sd_jet = sd_jets[k];

	
	
    // print the jet kinematic properties:
    //cout << "Original jet: pt=" << jetorigin.pt() << ", y=" << jetorigin.rap() << ", phi=" << jetorigin.phi() << ", mass=" << jetorigin.m() << endl;		 
		 
    // if soft-drop failed, just skip all the information
    if (sd_jet==0){
      cout << "Soft-dropped jet: zero PseudoJet" << endl;
      continue;
    }

    //cout << "Soft-dropped jet: pt=" << sd_jet.pt() << ", y=" << sd_jet.rap() << ", phi=" << sd_jet.phi() << ", mass=" << sd_jet.m() << endl;

    // show the last splitting information
    //cout << "soft-drop value = " << sd_jet.structure_of<contrib::SoftDropTagger>().soft_drop()<< endl;
    //cout << "mass-drop value = " << sd_jet.structure_of<contrib::SoftDropTagger>().mu()<< endl;

    //cout << "mass value = " << sd_jet.structure_of<contrib::SoftDropTagger>().mass()/Jpt<< endl;	
	
    // get the 2 prongs in the tagged jet
    //cout << "Distance between the 2 prongs=" << sd_jet.structure_of<contrib::SoftDropTagger>().Rg() << endl;

    // show the energy loss
    //cout << "Energy loss = " << 1-sd_jet.pt()/jetorigin.pt() << endl;

    // show the z_drop_max
    // (Commented out because still in beta testing)
    //cout << "# z_drop_max = " << sd_jet.structure_of<contrib::SoftDropTagger>().z_drop_max() << endl;

    //cout << endl;

	        //double zg=sd_jet.structure_of<contrib::SoftDropTagger>().soft_drop();
			
	        double zg=sd_jet.structure_of<contrib::SoftDropTagger>().zg();			
			
			double massjpt=sd_jet.structure_of<contrib::SoftDropTagger>().mass()/Jpt;

			double Rg=sd_jet.structure_of<contrib::SoftDropTagger>().Rg();

			//double massjpt=sd_jet.structure_of<contrib::SoftDropTagger>().mass()/Jptinitial;			
			
			
			//cout<<"Jptinitial"<<" "<<Jptinitial<<" "<<"zg"<<" "<<zg<<" "<<"--------------------------------------"<<endl;


            cout<<"ti"<<" "<<ti<<" "<<"zg"<<" "<<zg<<" "<<"massjpt"<<" "<<massjpt<<endl;
			
			if(zg>0.1 && massjpt>0.0)
            {

            numgroomedjet[ik]=numgroomedjet[ik]+1;
            							
//			numEvent=numEvent+1;	
//          if(numEvent>LBT->ncall) break;

            //cout<<"zg"<<" "<<zg<<endl;				
            //cout<<"massjpt"<<" "<<massjpt<<endl;
			
            hzg->Fill(zg,wwzg);
            hmassjpt->Fill(massjpt,wwmassjpt);
			
            SDzg[ik]->Fill(zg, wwzg);
	        SDmassjpt[ik]->Fill(massjpt, wwmassjpt);
	        SDmass[ik]->Fill(sd_jet.structure_of<contrib::SoftDropTagger>().mass(), wwgroomedmass);
            SDRg[ik]->Fill(Rg, wwRg);
            SDtheta[ik]->Fill(Rg, wwSDtheta);
            SDlnztheta[ik]->Fill(log(1.0/Rg), log(zg*Rg), wwSDlnztheta);
            //SDthetaE[ik]->Fill(Rg, subjetE, SDlnztheta);			
			
			
            }
	
            cout<<ik<<" "<<numgroomedjet[ik]<<endl;

	
    //cout<<"--------------------------- "<< Sjet1px << endl;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...jet grooming

	   	   
            Xjgamma=Xjgamma+Jpt/ptgamma;
		    PTjet=PTjet+Jpt;
			if(k==0)//leadingjet
			{
            PTgamma=PTgamma+ptgamma;			
		    nRjgamma=nRjgamma+1;
		    }

            if(timestep==ntimestep)	
			{
//			numEvent=numEvent+1;
			}			
			
//          if(numEvent>LBT->ncall) break;

			double pxj;
			double pyj;
			double pzj;
			double p0j;

			double ptj;

            double Pphi2;
			
            double Ptanphi;
            double Pphia;				

            double deltaphip;
			
			double Ppl;
			
			double Ppltemp[4]={0.0};
			
//..........energy flow and super shape...........................................................................................................super rho!
            for(int ii=1;ii<=LBT->np;ii++)
            {

            if(LBT->P[0][ii] == 0) continue;
			
			pxj=LBT->P[1][ii];
			pyj=LBT->P[2][ii];
			pzj=LBT->P[3][ii];
			p0j=LBT->P[0][ii];
			
			Ppltemp[1]=pxj;
			Ppltemp[2]=pyj;
			Ppltemp[3]=pzj;			
			Ppltemp[0]=p0j;			
			
			LBT->rotate(Jx,Jy,Jz,Ppltemp,1);

			Ppl=Ppltemp[3];

			LBT->rotate(Jx,Jy,Jz,Ppltemp,-1);			
						
			ptj=sqrt(pow(pxj,2)+pow(pyj,2));
			
            Ptanphi=abs(pyj/pxj);
            Pphia=atan(Ptanphi);			

            if (pxj>0 && pyj>0)
            {
            Pphi2=Pphia;
            }
            if (pxj<0 && pyj>0)
            {
            Pphi2=3.14159-Pphia;
            }
            if (pxj<0 && pyj<0) 
            {
            Pphi2=3.14159+Pphia;
            }
            if (pxj>0 && pyj<0) 
            {
            Pphi2=2.0*3.14159-Pphia;
            }									
			
            double deltaphip=Pphi2-Jphi2;					

               double deltaphipS=deltaphip;
			   
			   if(deltaphipS<-1.0/2.0*3.1415926)
			   {
			   deltaphipS=deltaphipS+2.0*3.1415926;
			   }

			
            if(deltaphip<0)
            {
            deltaphip=deltaphip+2.0*3.1415926;
            }

            if(deltaphip>3.1415926)
            {
            deltaphip=2.0*3.1415926-deltaphip;
            }							
			
/////////////////////////////////////////////////////			
            double deltaphip0=Pphi2-3.1415926/2.0;
            if(deltaphip0<0)
            {
            deltaphip0=deltaphip0+2.0*3.1415926;
            }

            if(deltaphip0>3.1415926)
            {
            deltaphip0=2.0*3.1415926-deltaphip0;
            }		
/////////////////////////////////////////////////////			
			
			
			
//////////////////////////////////////////////////
            double pJeta=1.0/2.0*log((p0j+pzj)/(p0j-pzj));
            double deltaetap = abs(pJeta-Jeta);
            //double deltaetap = abs(pJeta-0.0);

              double rparton = sqrt(pow(deltaphip,2)+pow(deltaetap,2));			
            //double rparton = sqrt(pow(deltaphip0,2)+pow(deltaetap,2));

			
            double ptparton = sqrt(pow(pxj,2)+pow(pyj,2));

			   //if(ptparton<0.5) continue;
			
               rholeadingJetsuper[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);
			   
			   if(LBT->CAT[ii]==2)
			   {
               rholeadingJetsuperCAT[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJetCAT[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);
			   }	
			
			   //if(ptparton>=0.5 && ptparton<1.0)
			   if(ptparton>=0.0 && ptparton<1.0)			   
			   {
               rholeadingJetsuper0510[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet0510[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   			   
			   }		
			   if(ptparton>=1.0 && ptparton<2.0)
			   {
               rholeadingJetsuper1020[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet1020[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=2.0 && ptparton<3.0)
			   {
               rholeadingJetsuper2030[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet2030[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=2.0 && ptparton<4.0)
			   {
               rholeadingJetsuper2040[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet2040[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=3.0 && ptparton<4.0)
			   {
               rholeadingJetsuper3040[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet3040[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }			
			   if(ptparton>=4.0 && ptparton<8.0)
			   {
               rholeadingJetsuper4080[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet4080[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }
			   
			   if(ptparton>=8.0)
			   {
               rholeadingJetsuper80[ik]->Fill(rparton, ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet80[ik]->Fill(deltaphipS, pJeta-Jeta, ptparton*wwSdndetadphi);			   
			   }

			
//////////////////////////////////////////////////

			
			if(ptj<=1.0)
			{						
            dndphiJet01[ik]->Fill(deltaphip, wwJdndphi);			
            dEtdphi01[ik]->Fill(deltaphip, ptj*wwJdndphi);
            dpldphi01[ik]->Fill(deltaphip, Ppl*wwJdndphi);			
            }
			
			if(ptj>1.0 && ptj<=2.0)
			{
            dndphiJet12[ik]->Fill(deltaphip, wwJdndphi);			
            dEtdphi12[ik]->Fill(deltaphip, ptj*wwJdndphi);
            dpldphi12[ik]->Fill(deltaphip, Ppl*wwJdndphi);			
            }			
			
			if(ptj>2.0 && ptj<=4.0)
			{
            dndphiJet24[ik]->Fill(deltaphip, wwJdndphi);			
            dEtdphi24[ik]->Fill(deltaphip, ptj*wwJdndphi);
            dpldphi24[ik]->Fill(deltaphip, Ppl*wwJdndphi);			
            }			
			
			if(ptj>4.0)
			{						
            dndphiJet40[ik]->Fill(deltaphip, wwJdndphi);			
            dEtdphi40[ik]->Fill(deltaphip, ptj*wwJdndphi);
            dpldphi40[ik]->Fill(deltaphip, Ppl*wwJdndphi);			
            }			

			
            }//for(int ii=1;ii<=np;ii++) positive end

			
            for(int ii=1;ii<=LBT->np;ii++)
            {
            if(LBT->P0[0][ii] == 0) continue;
			
			pxj=LBT->P0[1][ii];
			pyj=LBT->P0[2][ii];
			pzj=LBT->P0[3][ii];
			p0j=LBT->P0[0][ii];
			
			Ppltemp[1]=pxj;
			Ppltemp[2]=pyj;
			Ppltemp[3]=pzj;			
			Ppltemp[0]=p0j;			
			
			LBT->rotate(Jx,Jy,Jz,Ppltemp,1);

			Ppl=Ppltemp[3];

			LBT->rotate(Jx,Jy,Jz,Ppltemp,-1);			
			
			ptj=sqrt(pow(pxj,2)+pow(pyj,2));						

            Ptanphi=abs(pyj/pxj);
            Pphia=atan(Ptanphi);			

            if (pxj>0 && pyj>0)
            {
            Pphi2=Pphia;
            }
            if (pxj<0 && pyj>0)
            {
            Pphi2=3.14159-Pphia;
            }
            if (pxj<0 && pyj<0) 
            {
            Pphi2=3.14159+Pphia;
            }
            if (pxj>0 && pyj<0) 
            {
            Pphi2=2.0*3.14159-Pphia;
            }
						
            deltaphip=Pphi2-Jphi2;	

               double deltaphipS=deltaphip;
			   
			   if(deltaphipS<-1.0/2.0*3.1415926)
			   {
			   deltaphipS=deltaphipS+2.0*3.1415926;
			   }
			
			
            if(deltaphip<0)
            {
            deltaphip=deltaphip+2.0*3.1415926;
            }

            if(deltaphip>3.1415926)
            {
            deltaphip=2.0*3.1415926-deltaphip;
            }				
			
/////////////////////////////////////////////////////			
            double deltaphip0=Pphi2-3.1415926/2.0;
            if(deltaphip0<0)
            {
            deltaphip0=deltaphip0+2.0*3.1415926;
            }

            if(deltaphip0>3.1415926)
            {
            deltaphip0=2.0*3.1415926-deltaphip0;
            }		
/////////////////////////////////////////////////////


			
//////////////////////////////////////////////////
            double pJeta=1.0/2.0*log((p0j+pzj)/(p0j-pzj));
            double deltaetap = abs(pJeta-Jeta);
            //double deltaetap = abs(pJeta-0.0);


            double rparton = sqrt(pow(deltaphip,2)+pow(deltaetap,2));				
            //double rparton = sqrt(pow(deltaphip0,2)+pow(deltaetap,2));			  
            double ptparton = sqrt(pow(pxj,2)+pow(pyj,2));
			
			//if(ptparton<0.5) continue;
	  	  
            rholeadingJetsuper[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);	


			   //if(ptparton<0.5) continue;
			
               rholeadingJetsuper[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);
			   
			   //if(LBT->CAT[ii]==2)
			   {
               rholeadingJetsuperCAT[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJetCAT[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   			   
			   }	
			
			   //if(ptparton>=0.5 && ptparton<1.0)
			   if(ptparton>=0.0 && ptparton<1.0)			   
			   {
               rholeadingJetsuper0510[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet0510[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   			   
			   }		
			   if(ptparton>=1.0 && ptparton<2.0)
			   {
               rholeadingJetsuper1020[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet1020[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=2.0 && ptparton<3.0)
			   {
               rholeadingJetsuper2030[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet2030[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=2.0 && ptparton<4.0)
			   {
               rholeadingJetsuper2040[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet2040[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }
			   if(ptparton>=3.0 && ptparton<4.0)
			   {
               rholeadingJetsuper3040[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet3040[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }			
			   if(ptparton>=4.0 && ptparton<8.0)
			   {
               rholeadingJetsuper4080[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet4080[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }
			   
			   if(ptparton>=8.0)
			   {
               rholeadingJetsuper80[ik]->Fill(rparton, -ptparton/inclusive_jets[k].pt()*wwdRsuper);

			   SPTetaphileadingJet80[ik]->Fill(deltaphipS, pJeta-Jeta, -ptparton*wwSdndetadphi);			   
			   }

			
//////////////////////////////////////////////////		
			
			if(ptj<=1.0)
			{
            dndphiJet01[ik]->Fill(deltaphip, -wwJdndphi);			
            dEtdphi01[ik]->Fill(deltaphip, -ptj*wwJdndphi);	
            dpldphi01[ik]->Fill(deltaphip, -Ppl*wwJdndphi);			
            }
			
			if(ptj>1.0 && ptj<=2.0)
			{
            dndphiJet12[ik]->Fill(deltaphip, -wwJdndphi);			
            dEtdphi12[ik]->Fill(deltaphip, -ptj*wwJdndphi);
            dpldphi12[ik]->Fill(deltaphip, -Ppl*wwJdndphi);			
            }			
			
			if(ptj>2.0 && ptj<=4.0)
			{
            dndphiJet24[ik]->Fill(deltaphip, -wwJdndphi);			
            dEtdphi24[ik]->Fill(deltaphip, -ptj*wwJdndphi);	
            dpldphi24[ik]->Fill(deltaphip, -Ppl*wwJdndphi);			
            }			
			
			if(ptj>4.0)
			{
            dndphiJet40[ik]->Fill(deltaphip, -wwJdndphi);			
            dEtdphi40[ik]->Fill(deltaphip, -ptj*wwJdndphi);	
            dpldphi40[ik]->Fill(deltaphip, -Ppl*wwJdndphi);			
            }			
			
            }//for(int ii=1;ii<=np;i++) negative end
			
//..........energy flow end	


//..........energy flow and super shape...........................................................................................................super rho! end



			
//////////////////////////////////////////////////////
            double xjetgamma = Jpt/ptgamma;
            h100->Fill(xjetgamma,ww100);
			
//.................................................................................803

            h100Xjgamma->Fill(ptgamma,xjetgamma*ww10);
            h100Njgamma->Fill(ptgamma,ww10);
			h100Rjgamma->Fill(ptgamma,ww10);
			
//.................................................................................803			
			
//////////////////////////////////////////////////////

//            cout<<"Jpt**********************"<<" "<<Jpt<<endl;

            dndeleadingJet[ik]->Fill(JE, wwJdnde);
            dndptleadingJet[ik]->Fill(Jpt, wwJdndpt);		
            dndpt2leadingJet[ik]->Fill(Jpt, 2.0*Jpt*wwJdndpt);
            dndpesudoptleadingJet[ik]->Fill(Jpesudopt, wwJdndpt);
            dndpesudopt2leadingJet[ik]->Fill(Jpesudopt, 2.0*Jpesudopt*wwJdndpt);			
            dndthetaleadingJet[ik]->Fill(Jtheta, wwJdndtheta);
            dndphileadingJet[ik]->Fill(deltaphi, wwJdndphi);	
            dndetaleadingJet[ik]->Fill(Jeta, wwJdndeta);

//...jet Structure

        vector<fastjet::PseudoJet> constituents = inclusive_jets[k].constituents();
			
        for (unsigned int k1=0; k1<constituents.size(); k1++){

        double phiparton = inclusive_jets[k].phi_std()-constituents[k1].phi_std();
        if (phiparton<0) phiparton += 2*3.14159;
        if (phiparton>3.14159) phiparton = 2*3.14159-phiparton;
        double etaparton = inclusive_jets[k].eta()-constituents[k1].eta();

        double rparton = sqrt(pow(phiparton,2)+pow(etaparton,2));			  
        double ptparton = constituents[k1].pt();

          if (constituents[k1].user_index()>=0)
		  {			  
          rholeadingJet[ik]->Fill(rparton, constituents[k1].pt()/inclusive_jets[k].pt()*wwdR);	
          }

          if (constituents[k1].user_index()<0)
		  {			  
          rholeadingJet[ik]->Fill(rparton, -constituents[k1].pt()/inclusive_jets[k].pt()*wwdR);	
          }

		  
//...no negative partons		

          if (constituents[k1].user_index()>=0)
		  {

          double cosab = (constituents[k1].pz()*inclusive_jets[k].pz()+constituents[k1].py()*inclusive_jets[k].py()+constituents[k1].px()*inclusive_jets[k].px())/constituents[k1].E()/sqrt(pow(inclusive_jets[k].px(),2)+pow(inclusive_jets[k].py(),2)+pow(inclusive_jets[k].pz(),2));

          dndxplleadingJet[ik]->Fill(fabs(constituents[k1].E()*cosab)/inclusive_jets[k].E(), wwdndxpl);		  		  
          dndxpl0leadingJet[ik]->Fill(fabs(constituents[k1].E()*cosab)/energygamma, wwdndxpl0);
//...917???
          }

        }// for (unsigned int k1=0; k1<constituents.size(); k1++)

		}//if(deltaphi<Jdeltaphicut || Jpt<Jptcut || abs(Jeta)>Jetacut)
		
        }//if (inclusive_jets[k].user_index()>=0)

        }//for (unsigned int k = 0; k < 1; k++)

//...................................................................................................... fastjet Jet! end
		   
		   
		   
                }//if(timestep==ntimestep)



	        }//if(switchFastjet>0)

            }//if (fabs(ti - ik * dtstep) < 0.000001) // check time
			
            }//for (unsigned ik = 1; ik <= ntistep; ++ik)
		
        }//time loop end


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
//		cout<<"numEvent"<<" "<<numEvent<<" "<<"numGamma"<<" "<<numGamma<<" "<<"n"<<" "<<n<<endl;

/*		
		if(numgammajet>1.0)
		{	
		numgammawithjet=numgammawithjet+1.0;
        }
		
        if(numEvent>=LBT->ncall) break;
*/
		
		
//......cms cut
//		if(deltaphi<Jdeltaphicut || abs(Jeta)>Jetacut || Jpt<Jptcutmin || Jpt>Jptcutmax) continue;
						
		//cout<<"P[0][1]"<<" "<<LBT->P[0][1]<<endl;
		
//......time evolution end

//......output control
        int print=n%LBT->nprint;
	    if(print==0) 
		{
        cout << "n" << "    " << "ntest22" << "    " << "ntestrad" << endl;      
        cout << n << "  " << LBT->ntest22 << " " << LBT->ntestrad << endl;	
		
        cout << "n" << "    " << "nj" << "    " << "np" << endl;      
        cout << n << "  " << LBT->nj << " " << LBT->np << endl;			
		}

		//exit(1);
		
	}//for(int n=1;n<=LBT->ncall;n++)
//..event loop end



/*
//.................................................................................803

    h100Xjgamma->Divide(h100Njgamma);
	
	h100Rjgamma->Divide(h100Ngamma);

//.................................................................................803
*/

//..save the root file
    LBTjet->Write();
    LBTjet->Close();

    for(unsigned int ik = 1; ik <= ntistep; ik++)
    {
	inforLBT<<"numgroomedjet"<<" "<<"ti"<<" "<<ik<<endl;
	inforLBT<<numgroomedjet[ik]<<endl;
    }

    inforLBT.close();

    positiveLBT.close();	
    negativeLBT.close();	

    numpositiveLBT.close();	
    numnegativeLBT.close();


	
/*
    struct tm *local_end;
    time_t time_end;
    time_end=time(NULL);
    local_end=localtime(&time_end);
    
    char buf2[80];
    strftime(buf2,80,"Current Time: %Y-%m-%d %H:%M:%S",local_end);
    cout << "the program ends at:" << endl;
    cout << buf2 << endl;
*/
	
}


