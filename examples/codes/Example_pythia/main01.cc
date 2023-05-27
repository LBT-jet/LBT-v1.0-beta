// main72.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It compares SlowJet, FJcore and FastJet, showing that they
// find the same jets.

#include "Pythia8/Pythia.h"

// The FastJet3.h header enables automatic initialisation of
// fastjet::PseudoJet objects from Pythia8 Particle and Vec4 objects,
// as well as advanced features such as access to (a copy of)
// the original Pythia 8 Particle directly from the PseudoJet,
// and fastjet selectors that make use of the Particle properties.
// See the extensive comments in the header file for further details
// and examples.

//#include "Pythia8/FastJet3.h"


/////////////////////////////////////////////////////////////...fastjet

/// run it with    : ./01-basic < data/single-event.dat

#include "fastjet/ClusterSequence.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <fstream>
#include <math.h>
#include <vector> 

//#include "runglauber_v2.0.C"

using namespace std;
using namespace fastjet;

/////////////////////////////////////////////////////////////...fastjet end


/////////////////////////////////////////////////////////////...pythia+root

#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <fstream>
#include <math.h>
#include <vector>

// ROOT, for histogramming.
#include "TH1.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"
 
// ROOT, for saving file.
#include "TFile.h"

#include "TSystem.h"

using namespace Pythia8;

/////////////////////////////////////////////////////////////...pythia+root end

int KATT[2000000]={0};
double P[4][2000000]={0.0};

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

int ndrphoton;

int KATTgamma;

double energy;
double pt;

double VX;
double VY;		
double tishower;
double VZ;

double miu;
double sigma;
double range1;		
	   
double xxx;
double yyy;
	   
double fff;
		  
double rank0;

long  NUM1;                  //random number seed (negative integer)
		
		


//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
//////////////////////////////////////////


		
/////////////////////////////////////////////////////////////...fastjet

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

/////////////////////////////////////////////////////////////...fastjet end

int main() {

       NUM1=-33;

  gSystem->Load("libMathMore");

  // Number of events, generated and listed ones (for jets).
  int nEvent = 50000;
//  int nListJets = 3;
  
  double ptcut=20.0; 

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Process selection.

/*  
    //inclusive jet
    pythia.readString("HardQCD:all = on");
  
    pythia.readString("PhaseSpace:pTHatMin = 100.");    
*/
	
/*    //Dijet
    pythia.readString("HardQCD:gg2gg = on");  
    pythia.readString("HardQCD:gg2qqbar = on");  
    pythia.readString("HardQCD:qg2qg = on");  
    pythia.readString("HardQCD:qq2qq = on");  
    pythia.readString("HardQCD:qqbar2gg = on");  
    pythia.readString("HardQCD:qqbar2qqbarNew = on");

    pythia.readString("PhaseSpace:pTHatMin = 3.");
  
    pythia.readString("HardQCD:nQuarkNew = 3");
*/
  
    //gamma-jet

    pythia.readString("HardQCD:all = off");
	
    pythia.readString("PromptPhoton:qg2qgamma = on");
    pythia.readString("PromptPhoton:qqbar2ggamma = on");
    pythia.readString("PromptPhoton:gg2ggamma = on");

    pythia.readString("PhaseSpace:pTHatMin = 10.0");
    
/*	
    //Heavy-flavor
//  pythia.readString("HardQCD:gg2ccbar = on");  
//  pythia.readString("HardQCD:qqbar2ccbar = on");  
    pythia.readString("HardQCD:hardccbar = on");  
//  pythia.readString("HardQCD:gg2bbbar = on");  
//  pythia.readString("HardQCD:qqbar2bbbar = on");  
    pythia.readString("HardQCD:hardbbbar = on")  
  
    pythia.readString("PhaseSpace:pTHatMin = 200.");  
  
*/
  
  // Process Level.  
  
    pythia.readString("PartonLevel:all = on");   
    pythia.readString("PartonLevel:MPI = on");   
    pythia.readString("PartonLevel:ISR = on");   
    pythia.readString("PartonLevel:FSR = on");   
    pythia.readString("PartonLevel:FSRinProcess = on");
    pythia.readString("PartonLevel:FSRinResonances = on");	
    pythia.readString("PartonLevel:earlyResDec = on");		

    pythia.readString("HadronLevel:all = off");   

    // PDF selection
    pythia.readString("PDFinProcess:nQuarkIn = 3");

	// Number of allowed quark flavours in g → q qbar branchings (phase space permitting)
    pythia.readString("TimeShower:nGluonToQuark = 3");	

    pythia.readString("Random:setSeed = on");	
    pythia.readString("Random:Seed = 0");	

/*  
  // No event record printout.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");

*/ 
 
  // LHC initialization.
//  pythia.readString("Beams:idA = 2212");  
//  pythia.readString("Beams:idB = 2212");  
//  pythia.readString("Beams:eCM = 2760.0");
    pythia.readString("Beams:eCM = 5020.0");
//  pythia.readString("Beams:eCM = 200.0");
    pythia.init();

	int Kpythiajet=0;
	int njetparton=0;

///////////////////////////////////////////////////////

   // save datafile (parton information: flavor and four momentum)

   ofstream partoninfor("pythiaparton.data");
   
   ofstream gammainfor("gamma.dat");
   ofstream numparton("numpt.dat");
   ofstream partonrecombine("parton.dat");   
   

   // Create file on which histogram(s) can be saved.
   TFile* outFile = new TFile("pythia.root", "RECREATE");

   // Book histogram.
   TH1D *h11 = new TH1D("h11","rapidity distribution", 100, -10, 10);
   TH1D *h12 = new TH1D("h12","azimuthal distribution", 100, -3.14, 3.14);

   TH1D *h13 = new TH1D("h13","energy distribution", 200, 0, 200);
   TH1D *h14 = new TH1D("h14","transverse momentum distribution", 200, 0, 200);
   
   TH1D *h21 = new TH1D("h21","gamma rapidity distribution", 100, -10, 10);
   TH1D *h22 = new TH1D("h22","gamma azimuthal distribution", 50, -3.14, 3.14);

   TH1D *h23 = new TH1D("h23","gamma energy distribution", 200, 0, 400);
   TH1D *h24 = new TH1D("h24","gamma transverse momentum distribution", 200, 0, 400);   


   double ww11=1/(20./100)/nEvent;
   double ww12=1/(2.*3.14/100)/nEvent;
   double ww13=1/(200./200)/nEvent;
   double ww14=1/(200./200)/nEvent;
   
   double ww21=1/(20./100)/nEvent;
   double ww22=1/(2.*3.14/50)/nEvent; 
   double ww23=1/(300./200)/nEvent;
   double ww24=1/(300./200)/nEvent;


   
   int numEvent=0;

//////////////////////////////////////////////////   
   partonrecombine<<setiosflags(ios::left)<<setprecision(10);   
//////////////////////////////////////////////////
     
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; ; ++iEvent) {
  
       if(numEvent>nEvent) break;

       //partoninfor datafile mark
//       partoninfor<<"#"<<endl;
//       partoninfor<<iEvent<<endl;

       // vector<fastjet::PseudoJet> input_particles;

       // clock_t befGen = clock();
       if (!pythia.next()) continue;
       // Pythia running now, inside the nEvent loop
    
      njetparton=0;

//      partoninfor<<"Tracks in this event:"<<pythia.event.size()<<endl;
//      pythia.info.list(partoninfor);
//      pythia.event.list(partoninfor);

		ndrphoton=0;

		ptgamma=0.0;

		P0gamma=0.0;
	    P1gamma=0.0;
	    P2gamma=0.0;
	    P3gamma=0.0;

        int IDgamma=0;
		
///////////////////////////////////////////////////////  
    for (int i = 0; i < pythia.event.size(); ++i) //first find the direct photon
    {
        
//  cout<< i << endl;

        if(pythia.event[i].isFinal()) //The prompt photon may not be a final particle		
		{
		
		if(event[i].id()==22)
		//if(event[i].id()==22 && event[i].status()==23)		
		{
		
	    ndrphoton=ndrphoton+1;

//////////////////////////////////////////////		
 //       cout<<"Gamma!"<<"  "<<ndrphoton<<endl;
//////////////////////////////////////////////
		
	    KATTgamma=event[i].id();

        double ptgammathis=sqrt(pow(event[i].px(),2)+pow(event[i].py(),2));
		
	    //if(P0gamma<event[i].e())
	    if(ptgamma<ptgammathis)	
		{

        IDgamma=i;
	
	    P0gamma=event[i].e();
	    P1gamma=event[i].px();
	    P2gamma=event[i].py();
	    P3gamma=event[i].pz();
		
        etagamma=1.0/2.0*log((event[i].e()+event[i].pz())/(event[i].e()-event[i].pz()));		
        tanphigamma=P1gamma/P2gamma;
        phiagamma=atan(tanphigamma);

        if(P2gamma>=0.0)
		{
	    phi1gamma=phiagamma;
		}
		if(P2gamma<0.0 && P1gamma>=0.0)
		{
        phi1gamma=3.1415926+phiagamma;		
		}
		if(P2gamma<0.0 && P1gamma<0.0)
        {
        phi1gamma=-3.1415926+phiagamma;		
		}

        energygamma=P0gamma;
        ptgamma=sqrt(pow(event[i].px(),2)+pow(event[i].py(),2));				
		
		} //if(P0gamma<event[i].e())
		
	    } //if(event[i].id()==22 && event[i].status()==23)
		
		} //if(pythia.event[i].isFinal())
		
		} //for (int i = 0; i < pythia.event.size(); ++i)

//////////////////////////////////////////////		
//		cout<<"ptgamma"<<" "<<ptgamma<<endl;
//////////////////////////////////////////////
		
		if(ndrphoton==0 || ptgamma<20.0 || ptgamma>40.0 || abs(etagamma)>0.67) continue;			
		
//////////////////////////////////////////////		
		cout<<"numEvent"<<" "<<numEvent<<" "<<iEvent<<endl;
//////////////////////////////////////////////
		
        numEvent=numEvent+1;
		
        gammainfor<<numEvent<<" "<<P1gamma<<" "<<P2gamma<<" "<<P3gamma<<" "<<P0gamma<<endl;		
		
        h21->Fill(etagamma,ww21);
        h22->Fill(phi1gamma,ww22);
        h23->Fill(energygamma,ww23);
        h24->Fill(ptgamma,ww24);

		
    for (int i = 0; i < pythia.event.size(); ++i)  //the final partons
    {
	    
        VZ=0.0;
		tishower=0.0;
		
		//if(abs(event[i].id())==22) continue;		

		//if(abs(event[i].id())==22) continue;

        if(i==IDgamma) continue;
		
        if(pythia.event[i].isFinal())
        //if(pythia.event[i].isFinal() && pythia.event[i].isCharged()) //The prompt photon may not be a final particle		
		{
		
		//if(abs(event[i].id())==1 || abs(event[i].id())==2 || abs(event[i].id())==3 || abs(event[i].id())==21)
		//{
		
    	if(Kpythiajet==1) // switch on the Pythia-Fastjet
	    {

	    }
        else
        {

        njetparton=njetparton+1;

	    KATT[njetparton]=event[i].id();

	    P[0][njetparton]=event[i].e();
	    P[1][njetparton]=event[i].px();
	    P[2][njetparton]=event[i].py();
	    P[3][njetparton]=event[i].pz();
		
//		double sss=pow(event[i].e(),2)-(pow(event[i].px(),2)+pow(event[i].py(),2)+pow(event[i].pz(),2));

//		double sss=pow(P[0][njetparton],2)-(pow(P[1][njetparton],2)+pow(P[2][njetparton],2)+pow(P[3][njetparton],2));
		
//		cout<<"sss = "<<sss<<endl;
		

        // Fill in histograms. End event loop.

        eta=1.0/2.0*log((event[i].e()+event[i].pz())/(event[i].e()-event[i].pz()));		

        tanphi=P[1][njetparton]/P[2][njetparton];
        phia=atan(tanphi);
	
        if(P[2][njetparton]>=0.0)
		{
	    phi1=phia;
		}
		if(P[2][njetparton]<0.0 && P[1][njetparton]>=0.0)
		{
        phi1=3.1415926+phia;		
		}
		if(P[2][njetparton]<0.0 && P[1][njetparton]<0.0)
        {
        phi1=-3.1415926+phia;
		}
		
        energy=P[0][njetparton];
        pt=sqrt(pow(event[i].px(),2)+pow(event[i].py(),2));
        
        h11->Fill(eta,ww11);
        h12->Fill(phi1-phi1gamma,ww12);
        h13->Fill(energy,ww13);
        h14->Fill(pt,ww14);

            miu=0.0;
            sigma=0.2;
            range1=1.0;		
		
		
s5:	        VX=ran33(&NUM1)*range1-0.5*range1;

	        fff=exp((-0.5)*pow(((VX-miu)/sigma),2));
		  
	        rank0=ran33(&NUM1);

	       if(rank0>fff)
		   {
		   goto s5;
           }	   

s6:	        VY=ran33(&NUM1)*range1-0.5*range1;

	        fff=exp((-0.5)*pow(((VY-miu)/sigma),2));

	        rank0=ran33(&NUM1);
				  
	       if(rank0>fff)
		   {
		   goto s6;
		   }		
		   
//////////////////////////////////////////////////////////
           VX=0.0;
           VY=0.0; 
//////////////////////////////////////////////////////////		   



















/////////////////////////////////////////////////////////////////////////////////////






//......formation time +

        int ishower=0;

        double p0[4]={0.0};        
        double p4[4]={0.0};
		
		double qt;
		double timeplus=0.0;
		double timeplus_0=0.0;
		double timeplus_0_1=0.0;
		
        int IDmom1, IDmom2;
		int timebreaker = 0;

        int IDmom;

        int IDmom0 = i;
		
		//cout<<"------------------------------------"<<" "<<i<<endl;

        while(timebreaker == 0){
		
		      int IDiii = IDmom0;

              if(abs(event[IDiii].status())==23 || abs(event[IDiii].status())==21 || abs(event[IDiii].status())==12)
			  {
			  //cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++end"<<endl;				  

              ishower=42;

              timebreaker=1;			
		      }

              IDmom1=event[IDiii].mother1();
              IDmom2=event[IDiii].mother2();

              if(IDmom1==IDmom2 && IDmom1==0){
		      
			  //cout<<"IDmom1==IDmom2 && IDmom1==0 timebreak"<<endl;

              ishower=40;
			  
              timebreaker=1;
		
		      }//if(IDmom1==IDmom2 && IDmom1==0)

              if(IDmom1==IDmom2 && IDmom1>0){
		
              IDmom0=IDmom1;
		
		      }//if(IDmom1==IDmom2 && IDmom1>0)

              if(IDmom1>0 && IDmom2==0){

              IDmom0=IDmom1;
			  
			  
//........................................			  

              double IDdaughter1=event[IDmom0].daughter1();
              double IDdaughter2=event[IDmom0].daughter2();
			  
			  if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){


              p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
              p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
              p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
              p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();
			  

              double x_split=event[IDiii].e()/p4[0];
           
//........................................


              //double x_split=event[IDiii].e()/event[IDmom0].e();
			  
			  if(x_split>1) x_split=1.0/x_split;


              p0[0]=event[IDiii].e();
              p0[1]=event[IDiii].px();
              p0[2]=event[IDiii].py();
              p0[3]=event[IDiii].pz(); 
  			  
              //p4[0]=event[IDmom0].e();
              //p4[1]=event[IDmom0].px();
              //p4[2]=event[IDmom0].py();
              //p4[3]=event[IDmom0].pz();
			  
			  double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
			  double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			  
			  

			  //cout<<"daughter"<<" "<<p0[1]<<" "<<p0[2]<<" "<<p0[3]<<" "<<p0[0]<<" "<<endl;
			  //cout<<"mother"<<" "<<p4[1]<<" "<<p4[2]<<" "<<p4[3]<<" "<<p4[0]<<" "<<endl;

              rotate(p4[1],p4[2],p4[3],p0,1);
			  
              qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
              
			  rotate(p4[1],p4[2],p4[3],p0,-1);

              double kt_daughter=qt;
			  
			  double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

              if(x_split<0.5){
			  
			  //if(kt_daughter > 0.0001 && Q2 > 1.0)
			  if(kt_daughter > 0.0001)
			  {
              //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
              timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
			  }
			  
			  }

              }
			  
		      }//if(IDmom1>0 && IDmom2==0)
		
              if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0){

                    if(event[IDmom1].e()>event[IDmom2].e())
					{
                    IDmom0=IDmom1;			
		            }
                    if(event[IDmom1].e()<=event[IDmom2].e())
					{
                    IDmom0=IDmom2;			
		            }

//........................................			  

              double IDdaughter1=event[IDmom0].daughter1();
              double IDdaughter2=event[IDmom0].daughter2();
			  
			  if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){


              p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
              p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
              p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
              p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();
			  

              double x_split=event[IDiii].e()/p4[0];
           
//........................................



                    //double x_split=event[IDiii].e()/event[IDmom0].e();
					
			        if(x_split>1) x_split=1.0/x_split;

                    p0[0]=event[IDiii].e();
                    p0[1]=event[IDiii].px();
                    p0[2]=event[IDiii].py();
                    p0[3]=event[IDiii].pz(); 
  			  
                    //p4[0]=event[IDmom0].e();
                    //p4[1]=event[IDmom0].px();
                    //p4[2]=event[IDmom0].py();
                    //p4[3]=event[IDmom0].pz();

                    rotate(p4[1],p4[2],p4[3],p0,1);
			  
                    qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
              
			        rotate(p4[1],p4[2],p4[3],p0,-1);

                    double kt_daughter=qt;
					
			        double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

                    if(x_split<0.5){
			  
			        //if(kt_daughter > 0.0001 && Q2 > 1.0)
			        if(kt_daughter > 0.0001)
			        {
                    //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                    timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
		            }
					
					}

              }
			  
		      }//if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0)

              if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0){
			
                    if(event[IDmom1].e()>event[IDmom2].e())
			        {
                    IDmom0=IDmom1;			
		            }
                    if(event[IDmom1].e()<=event[IDmom2].e())
			        {
                    IDmom0=IDmom2;			
		            }			


//........................................			  

              double IDdaughter1=event[IDmom0].daughter1();
              double IDdaughter2=event[IDmom0].daughter2();
			  
			  if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){


              p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
              p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
              p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
              p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();
			  

              double x_split=event[IDiii].e()/p4[0];
           
//........................................


                    //double x_split=event[IDiii].e()/event[IDmom0].e();

			        if(x_split>1) x_split=1.0/x_split;

                    p0[0]=event[IDiii].e();
                    p0[1]=event[IDiii].px();
                    p0[2]=event[IDiii].py();
                    p0[3]=event[IDiii].pz(); 
  			  
                    //p4[0]=event[IDmom0].e();
                    //p4[1]=event[IDmom0].px();
                    //p4[2]=event[IDmom0].py();
                    //p4[3]=event[IDmom0].pz();

                    rotate(p4[1],p4[2],p4[3],p0,1);
			  
                    qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
              
			        rotate(p4[1],p4[2],p4[3],p0,-1);

                    double kt_daughter=qt;
			        
				    double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

                    if(x_split<0.5){
			  
			        //if(kt_daughter > 0.0001 && Q2 > 1.0)
			        if(kt_daughter > 0.0001)					
			        {						
                    //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                    timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;					
					}

                    }

              }

		      }//if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0)		




              if(abs(event[IDmom0].status())==23 || abs(event[IDmom0].status())==21 || abs(event[IDmom0].status())==12){

              cout<<"event[IDmom0].status()"<<" "<<event[IDmom0].status()<<endl;

              ishower=abs(event[IDmom0].status());

              if(abs(event[IDmom0].status()) == 23)
			  {
			  //ishower=23;
			  //cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++23"<<endl;  
			  }

              double x_split=event[IDiii].e()/event[IDmom0].e();
			  
			  if(x_split>1) x_split=0;


              p0[0]=event[IDiii].e();
              p0[1]=event[IDiii].px();
              p0[2]=event[IDiii].py();
              p0[3]=event[IDiii].pz(); 
  			  
              p4[0]=event[IDmom0].e();
              p4[1]=event[IDmom0].px();
              p4[2]=event[IDmom0].py();
              p4[3]=event[IDmom0].pz();

              rotate(p4[1],p4[2],p4[3],p0,1);
			  
              qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
              
			  rotate(p4[1],p4[2],p4[3],p0,-1);

              double kt_daughter=qt;

			  double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
			  
			  if(kt_daughter > 0.0001 && Q2 > 1.0)
			  {
			  
			  timeplus_0=2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;

              }

			  timeplus_0_1=2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;

              timebreaker=1;			
		      }

        }//while(timebreaker == 0)
/*
mother1 = mother2 = 0: for lines 0 - 2, where line 0 represents the event as a whole, and 1 and 2 the two incoming beam particles;
mother1 = mother2 > 0: the particle is a "carbon copy" of its mother, but with changed momentum as a "recoil" effect, e.g. in a shower;
mother1 > 0, mother2 = 0: the "normal" mother case, where it is meaningful to speak of one single mother to several products, in a shower or decay;
mother1 < mother2, both > 0, for abs(status) = 81 - 86: primary hadrons produced from the fragmentation of a string spanning the range from mother1 to mother2, so that all partons in this range should be considered mothers; and analogously for abs(status) = 101 - 106, the formation of R-hadrons;
mother1 < mother2, both > 0, except case 4: particles with two truly different mothers, in particular the particles emerging from a hard 2 → n interaction.
mother2 < mother1, both > 0: particles with two truly different mothers, notably for the special case that two nearby partons are joined together into a status 73 or 74 new parton, in the g + q → q case the q is made first mother to simplify flavour tracing.
*/


        //if(abs(ishower) != 23) timeplus=0.0;


//......formation time -












////////////////////////////////////////////////////////////////////////////////////




        tishower=timeplus;
  
		
		//tishower=2.0*event[i].e()/pow(pt,2)*0.197;
		VZ=VZ+event[i].pz()/event[i].e()*tishower;
		
		VX=VX+event[i].px()/event[i].e()*tishower;		
		VY=VY+event[i].py()/event[i].e()*tishower;		
		
//      cout parton information

//        partoninfor<<KATT[njetparton]<<" "<<P[1][njetparton]<<" "<<P[2][njetparton]<<" "<<P[3][njetparton]<<" "<<P[0][njetparton]<<endl;


//        fout<<setiosflags(ios::left)<<setprecision(10);
//        //fout<<"x: "<<x<<"\t"<<"y: "<<y<<endl;
//        fout<<setw(12)<<x<<"\t"<<setw(12)<<y<<endl;

        partonrecombine<<setw(12)<<numEvent<<"  "<<setw(12)<<event[i].id()<<"  "<<setw(12)<<event[i].px()<<"  "<<setw(12)<<event[i].py()<<"  "<<setw(12)<<event[i].pz()<<"  "<<setw(12)<<event[i].e()<<"  "<<setw(12)<<ishower<<"  "<<setw(12)<<VX<<"  "<<setw(12)<<VY<<"  "<<setw(12)<<VZ<<"  "<<setw(12)<<tishower<<endl;
	
		
//        partonrecombine<<numEvent<<"  "<<event[i].id()<<"  "<<event[i].px()<<"  "<<event[i].py()<<"  "<<event[i].pz()<<"  "<<event[i].e()<<"  "<<VX<<"  "<<VY<<"  "<<VZ<<"  "<<tishower<<endl;
	
        //} //if(abs(event[i].id())==1 && abs(event[i].id())==2 && abs(event[i].id())==3 && abs(event[i].id())==21)
		
        } //if(pythia.event[i].isFinal())
		
	} //end Kpythiajet!=1 all partons are in

    } //for (int i = 0; i < pythia.event.size(); ++i)

	numparton<<numEvent<<"  "<<njetparton<<endl;	
	
	} //iEvent end
	
	cout<<"Pythiaend"<<endl;

	outFile->cd();
	h11->Write();
	h12->Write();
    h13->Write();
    h14->Write();
	
	h21->Write();
	h22->Write();
    h23->Write();
    h24->Write();	
	
//////////////////////////////////////////////////////

    gammainfor.close();
    numparton.close();
    partonrecombine.close();
	partoninfor.close();
   
//////////////////////////////////////////////////////



  // Done.
  cout << "job ends!" << endl;
  return 0;
}


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran33(long *idum)

{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) { 
	if (-(*idum) < 1) *idum=1; 
	else *idum = -(*idum);
	for (j=NTAB+7;j>=0;j--) { 
	  k=(*idum)/IQ1;
	  *idum=IA1*(*idum-k*IQ1)-k*IR1;
	  if (*idum < 0) *idum += IM1;
	  if (j < NTAB) iv[j] = *idum;
	}
	iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 
  iy=iv[j]-idum2; 
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}
















void rotate(double px,double py,double pz,double pr[4],int icc){
  //     input:  (px,py,pz), (wx,wy,wz), argument (i)
  //     output: new (wx,wy,wz)
  //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
  //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

	   
  double wx,wy,wz,E,pt,w,cosa,sina,cosb,sinb;
  double wx1,wy1,wz1;	   	   

  wx=pr[1];
  wy=pr[2];
  wz=pr[3];

  E=sqrt(px*px+py*py+pz*pz);
  pt=sqrt(px*px+py*py);

  w=sqrt(wx*wx+wy*wy+wz*wz);

//  if(pt==0)
  if(pt<1e-6)
	{
	  cosa=1;
	  sina=0;
	} 
  else
	{
	  cosa=px/pt;
	  sina=py/pt;
	}

  if(E>1e-6) {

      cosb=pz/E;
      sinb=pt/E;
    
      if(icc==1) {
    	  wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
    	  wy1=-wx*sina+wy*cosa;
    	  wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
      } else {
    	  wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
    	  wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
    	  wz1=-wx*sinb+wz*cosb;
      }
      wx=wx1;
      wy=wy1;
      wz=wz1;
  } else {
      cout << "warning: small E in rotation" << endl;
  }

  pr[1]=wx;
  pr[2]=wy;
  pr[3]=wz;      

//  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}

