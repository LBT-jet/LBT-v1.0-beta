//
// $Id: runglauber_v2.0.C 17 2014-08-11 19:46:37Z loizides $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TNamed.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TVector3.h>
#include <Math/SpecFuncMathMore.h>

using namespace std;

#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_
#endif

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA="Au",
		      const char *sysB="Au",
		      const Double_t signn=42,
		      const Double_t sigwidth=-1,
                      const Double_t mind=0.4,
                      const char *fname="glau_auau_ntuple.root");

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA="Au",           
                        const char *sysB="Au",           
                        const Double_t signn=42,           
    		        const Double_t sigwidth=-1,
                        const Double_t mind=0.4,
                        const Bool_t verbose=0,
                        const char *fname="glau_auau_nucleons.root");

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs  = 0.4,
                       const char *sysA     = "p",
                       const char *sysB     = "Pb",
                       const Double_t signn = 70,
                       const Double_t mind  = 0.4,
                       const char *fname    = "glau_ppb_smeared_ntuple.root");

//---------------------------------------------------------------------------------
class TGlauNucleon : public TNamed
{
private:
  Double32_t fX;            //Position of nucleon
  Double32_t fY;            //Position of nucleon
  Double32_t fZ;            //Position of nucleon
  Int_t      fType;         //0 = neutron, 1 = proton
  Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
  Int_t      fNColl;        //Number of binary collisions

public:
  TGlauNucleon() : fX(0), fY(0), fZ(0), fInNucleusA(0), fNColl(0) {}
  virtual   ~TGlauNucleon() {}
  
  void       Collide()                                  {fNColl++;}
  Double_t   Get2CWeight(Double_t x) const              {return 2.*(0.5*(1-x)+0.5*x*fNColl);}
  Int_t      GetNColl()              const              {return fNColl;}
  Int_t      GetType()               const              {return fType;}
  Double_t   GetX()                  const              {return fX;}
  Double_t   GetY()                  const              {return fY;}
  Double_t   GetZ()                  const              {return fZ;}
  Bool_t     IsInNucleusA()          const              {return fInNucleusA;}
  Bool_t     IsInNucleusB()          const              {return !fInNucleusA;}
  Bool_t     IsSpectator()           const              {return !fNColl;}
  Bool_t     IsWounded()             const              {return fNColl;}
  void       Reset()                                    {fNColl=0;}
  void       SetType(Bool_t b)                          {fType = b;}
  void       SetInNucleusA()                            {fInNucleusA=1;}
  void       SetInNucleusB()                            {fInNucleusA=0;}
  void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
  void       RotateXYZ(Double_t phi, Double_t theta);

  ClassDef(TGlauNucleon,2) // TGlauNucleon class
};

void TGlauNucleon::RotateXYZ(Double_t phi, Double_t theta)
{
  TVector3 v(fX,fY,fZ);
  TVector3 vr;
  vr.SetMagThetaPhi(1,theta,phi);
  v.RotateUz(vr);
  fX = v.X();
  fY = v.Y();
  fZ = v.Z();
}

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
private:
  Int_t      fN;                   //Number of nucleons
  Int_t      fZ;                   //Nuclear charge
  Double_t   fR;                   //Parameters of function
  Double_t   fA;                   //Parameters of function
  Double_t   fW;                   //Parameters of function
  Double_t   fBeta2;               //Beta2 (deformed nuclei) 
  Double_t   fBeta4;               //Beta4 (deformed nuclei) 
  Double_t   fMinDist;             //Minimum separation distance
  Int_t      fF;                   //Type of radial distribution
  Int_t      fTrials;              //Store trials needed to complete nucleus
  TF1*       fFunction;            //Probability density function rho(r)
  TF2*       fFunction2;           //Probability density function rho(r,theta) for deformed nuclei
  TObjArray* fNucleons;            //Array of nucleons
  Double_t   fPhiRot;              //Angle phi for nucleus
  Double_t   fThetaRot;            //Angle theta for nucleus
  Double_t   fHe3Arr[20000][3][3]; //Array of events, 3 nucleons, 3 coordinates
  Int_t      fHe3Counter;          //Event counter

  void       Lookup(const char* name);
  
public:
  TGlauNucleus(const char* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
  virtual ~TGlauNucleus();
  
  using      TObject::Draw;
  void       Draw(Double_t xs, Int_t col);
  Double_t   GetA()             const {return fA;}
  TF1*       GetFunction()      const {return fFunction;}
  Int_t      GetN()             const {return fN;}
  Double_t   GetR()             const {return fR;}
  TObjArray *GetNucleons()      const {return fNucleons;}
  Double_t   GetPhiRot()        const {return fPhiRot;}
  Double_t   GetThetaRot()      const {return fThetaRot;}
  Int_t      GetTrials()        const {return fTrials;}
  Double_t   GetW()             const {return fW;}
  Double_t   GetMinDist()       const {return fMinDist;}
  void       SetN(Int_t in)           {fN=in;}
  void       SetR(Double_t ir);
  void       SetA(Double_t ia);
  void       SetW(Double_t iw);
  void       SetMinDist(Double_t min) {fMinDist=min;}
  void       ThrowNucleons(Double_t xshift=0.);

  ClassDef(TGlauNucleus,2) // TGlauNucleus class
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
private:
  TGlauNucleus fANucleus;       //Nucleus A
  TGlauNucleus fBNucleus;       //Nucleus B
  Double_t     fXSect;          //Nucleon-nucleon cross section
  Double_t     fXSectOmega;     //StdDev of Nucleon-nucleon cross section
  Double_t     fXSectLambda;    //Jacobian from tot to inelastic (Strikman)
  Double_t     fXSectEvent;     //Event value of Nucleon-nucleon cross section
  TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
  TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
  TObjArray*   fNucleons;       //Array which joins Nucleus A & B
  Int_t        fAN;             //Number of nucleons in nucleus A
  Int_t        fBN;             //Number of nucleons in nucleus B
  TNtuple*     fNt;             //Ntuple for results (created, but not deleted)
  Double_t     fMeanX2;         //<x^2> of wounded nucleons
  Double_t     fMeanY2;         //<y^2> of wounded nucleons
  Double_t     fMeanXY;         //<xy> of wounded nucleons
  Double_t     fMeanXParts;     //<x> of wounded nucleons
  Double_t     fMeanYParts;     //<x> of wounded nucleons
  Double_t     fMeanXSystem;    //<x> of all nucleons
  Double_t     fMeanYSystem;    //<x> of all nucleons  
  Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
  Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
  Double_t     fB_MC;           //Impact parameter (b)
  Int_t        fEvents;         //Number of events with at least one collision
  Int_t        fTotalEvents;    //All events within selected impact parameter range
  Double_t     fBMin;           //Minimum impact parameter to be generated
  Double_t     fBMax;           //Maximum impact parameter to be generated
  Int_t        fMaxNpartFound;  //Largest value of Npart obtained
  Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
  Int_t        fNpartA;         //Number of wounded (participating) nucleons in Nucleus A
  Int_t        fNpartB;         //Number of wounded (participating) nucleons in Nucleus B
  Double_t     fSumW;           //Number of wounded (participating) nucleons in current event
  Double_t     fSumWA;          //Number of wounded (participating) nucleons in Nucleus A
  Double_t     fSumWB;          //Number of wounded (participating) nucleons in Nucleus B
  Int_t        fNcoll;          //Number of binary collisions in current event
  Double_t     fSx2;            //Variance of x of wounded nucleons
  Double_t     fSy2;            //Variance of y of wounded nucleons
  Double_t     fSxy;            //Covariance of x and y of wounded nucleons
  Double_t     fPsiN[10];       //Psi N
  Double_t     fEccN[10];       //Ecc N
  Double_t     f2Cx;            //Two-component x
  TF1         *fPTot;           //Cross section distribution

  Bool_t       CalcResults(Double_t bgen);
  Bool_t       CalcEvent(Double_t bgen);

public:
  TGlauberMC(const char* NA = "Pb", const char* NB = "Pb", Double_t xsect = 42, Double_t xsectsigma=0);
  virtual     ~TGlauberMC() {Reset();}
  
  void         Draw(Option_t* option);
  Double_t     GetB()               const {return fB_MC;}
  Double_t     GetBMin()            const {return fBMin;}
  Double_t     GetBMax()            const {return fBMax;}
  Double_t     GetEcc(Int_t i=2)    const {return fEccN[i];}
  Int_t        GetNcoll()           const {return fNcoll;}
  Int_t        GetNpart()           const {return fNpart;}
  Int_t        GetNpartA()          const {return fNpartA;}
  Int_t        GetNpartB()          const {return fNpartB;}
  Int_t        GetNpartFound()      const {return fMaxNpartFound;}
  TObjArray   *GetNucleons();
  TNtuple*     GetNtuple()          const {return fNt;}
  Double_t     GetPsi(Int_t i=2)    const {return fPsiN[i];}
  Double_t     GetSx2()             const {return fSx2;}    
  Double_t     GetSy2()             const {return fSy2;}    
  Double_t     GetSxy()             const {return fSxy;}    
  Double_t     GetTotXSect()        const;
  Double_t     GetTotXSectErr()     const;
  TF1*         GetXSectDist()       const {return fPTot;}
  Double_t     GetXSectEvent()      const {return fXSectEvent;}
  Bool_t       NextEvent(Double_t bgen=-1);
  void         Reset()                    {delete fNt; fNt=0; }
  void         Run(Int_t nevents,Double_t b=-1);
  void         SetBmin(Double_t bmin)     {fBMin = bmin;}
  void         SetBmax(Double_t bmax)     {fBMax = bmax;}
  void         Set2Cx(Double_t x)         {f2Cx = x;}
  void         SetMinDistance(Double_t d) {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}

  const TGlauNucleus* GetNucleusA() const {return &fANucleus;}
  const TGlauNucleus* GetNucleusB() const {return &fBNucleus;}

  static void  PrintVersion()             {cout << "TGlauberMC " << Version() << endl;}
  static const char *Version()            {return "v2.0";}

  ClassDef(TGlauberMC,2) // TGlauberMC class
};

//---------------------------------------------------------------------------------
void runAndSaveNtuple(const Int_t n,
                      const char *sysA,
                      const char *sysB,
                      const Double_t signn,
                      const Double_t sigwidth,
                      const Double_t mind,
                      const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  mcg->Run(n);
  TNtuple  *nt=mcg->GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if (nt) nt->Write();
  out.Close();
}

//---------------------------------------------------------------------------------
void runAndSaveNucleons(const Int_t n,                    
                        const char *sysA,           
                        const char *sysB,           
                        const Double_t signn,
                        const Double_t sigwidth,
                        const Double_t mind,
                        const Bool_t verbose,
                        const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn,sigwidth);
  mcg->SetMinDistance(mind);
  TFile *out=0;
  if (fname) 
    out=new TFile(fname,"recreate",fname,9);

  for (Int_t ievent=0; ievent<n;ievent++) {
    //get an event with at least one collision
    while (!mcg->NextEvent()) {}

    //access, save and (if wanted) print out nucleons
    TObjArray* nucleons=mcg->GetNucleons();
    if (!nucleons) continue;
    if (out)
      nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey);

    if (verbose) {
      cout<<endl<<endl<<"EVENT NO: "<<ievent<<endl;
      cout<<"B = "<<mcg->GetB()<<"  Npart = "<<mcg->GetNpart()<<endl<<endl;
      printf("Nucleus\t X\t Y\t Z\tNcoll\n");
      Int_t nNucls=nucleons->GetEntries();
      for (Int_t iNucl=0;iNucl<nNucls;iNucl++) {
	TGlauNucleon *nucl=(TGlauNucleon *)nucleons->At(iNucl);
	Char_t nucleus='A';
	if (nucl->IsInNucleusB()) nucleus='B';
	Double_t x=nucl->GetX();
	Double_t y=nucl->GetY();
	Double_t z=nucl->GetZ();
	Int_t ncoll=nucl->GetNColl();
	printf("   %c\t%2.2f\t%2.2f\t%2.2f\t%3d\n",nucleus,x,y,z,ncoll);
      }
    }
  }
  if (out) delete out;
}

//---------------------------------------------------------------------------------
void runAndSmearNtuple(const Int_t n,
                       const Double_t sigs,
                       const char *sysA,
                       const char *sysB,
                       const Double_t signn,
                       const Double_t mind,
                       const char *fname)
{
  // Run Glauber and store ntuple with smeared eccentricities in file.

  TGlauberMC *mcg = new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  TFile *out = TFile::Open(fname,"recreate",fname,9);
  if (!out)
    return;
  TNtuple *nt = new TNtuple("nt","nt",
			    "Npart:Ncoll:B:Psi2:Ecc2:Psi3:Ecc3:Psi2G:Ecc2G:Psi3G:Ecc3G:Sx2:Sy2:Sx2G:Sy2G");
  nt->SetDirectory(out);

  const Int_t NSAMP = 100;
  TF1 *rad = new TF1("rad","x*TMath::Exp(-x*x/(2.*[0]*[0]))",0.0,1.5);
  rad->SetParameter(0,sigs);

  for (Int_t ievent=0; ievent<n; ++ievent) {
    while (!mcg->NextEvent()) {}

    const TGlauNucleus *nucA   = mcg->GetNucleusA();
    const TObjArray *nucleonsA = nucA->GetNucleons();
    const Int_t AN             = nucA->GetN();
    const TGlauNucleus *nucB   = mcg-> GetNucleusB();
    const TObjArray *nucleonsB = nucB->GetNucleons();
    const Int_t BN             = nucB->GetN();
    Double_t sinphi[10] = {0};
    Double_t cosphi[10] = {0};
    Double_t rn[10]     = {0};
    Double_t ecc[10]    = {0};
    Double_t psi[10]    = {0};
    Double_t sx2g       = 0;
    Double_t sy2g       = 0;

    for (Int_t s=0; s<NSAMP; ++s) {
      Int_t ni = 0;
      Double_t xvals[1000] = {0};
      Double_t yvals[1000] = {0};
      for (Int_t i = 0; i<AN; ++i) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(nucleonsA->At(i));
	if (!nucleonA->IsWounded())
	  continue;
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	xvals[ni]   = nucleonA->GetX() + sr*TMath::Cos(sp);
	yvals[ni]   = nucleonA->GetY() + sr*TMath::Sin(sp);
	ni++;
      }
      for (Int_t i = 0; i<BN; ++i) {
	TGlauNucleon *nucleonB=(TGlauNucleon*)(nucleonsB->At(i));
	if (!nucleonB->IsWounded())
	  continue;
	Double_t sr = rad->GetRandom();
	Double_t sp = gRandom->Uniform(-TMath::Pi(), +TMath::Pi());
	xvals[ni]   = nucleonB->GetX() + sr*TMath::Cos(sp);
	yvals[ni]   = nucleonB->GetY() + sr*TMath::Sin(sp);
	ni++;
      }

      Double_t MeanXParts  = 0;
      Double_t MeanYParts  = 0;
      Double_t MeanXParts2 = 0;
      Double_t MeanYParts2 = 0;
      for (Int_t i = 0; i<ni; ++i) {
	MeanXParts  += xvals[i];
	MeanYParts  += yvals[i];
	MeanXParts2 += xvals[i]*xvals[i];
	MeanYParts2 += yvals[i]*yvals[i];
      }
      MeanXParts  /= ni;
      MeanYParts  /= ni;
      MeanXParts2 /= ni;
      MeanYParts2 /= ni;
      sx2g        += MeanXParts2-MeanXParts*MeanXParts;
      sy2g        += MeanYParts2-MeanYParts*MeanYParts;

      for (Int_t j = 1; j<9; ++j) {
	for (Int_t i = 0; i<ni; ++i) {
	  Double_t x   = xvals[i] - MeanXParts;
	  Double_t y   = yvals[i] - MeanYParts;
	  Double_t r   = TMath::Sqrt(x*x+y*y);
	  Double_t phi = TMath::ATan2(y,x);
	  cosphi[j] += TMath::Power(r,j)*TMath::Cos(j*phi);
	  sinphi[j] += TMath::Power(r,j)*TMath::Sin(j*phi);
	  rn[j]     += TMath::Power(r,j);
	}
      }
    }
    for (Int_t j = 1; j<9; ++j) {
      psi[j] = (TMath::ATan2(sinphi[j],cosphi[j]) + TMath::Pi())/j;
      ecc[j] = TMath::Sqrt(sinphi[j]*sinphi[j] + cosphi[j]*cosphi[j]) / rn[j];
    }

    Float_t v[15];
    v[0]  = mcg->GetNpart();
    v[1]  = mcg->GetNcoll();
    v[2]  = mcg->GetB();
    v[3]  = mcg->GetPsi(2);
    v[4]  = mcg->GetEcc(2);
    v[5]  = mcg->GetPsi(3);
    v[6]  = mcg->GetEcc(3);
    v[7]  = psi[2];
    v[8]  = ecc[2];
    v[9]  = psi[3];
    v[10] = ecc[3];
    v[11] = mcg->GetSx2();
    v[12] = mcg->GetSy2();
    v[13] = sx2g/NSAMP;
    v[14] = sy2g/NSAMP;
    nt->Fill(v);
  }

  out->Write();
  out->Close();
  delete out;
}

//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------
TGlauNucleus::TGlauNucleus(const char* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
  fN(iN),fR(iR),fA(ia),fW(iw),fBeta2(0),fBeta4(0),fMinDist(0.4),
  fF(0),fTrials(0),fFunction(ifunc),fFunction2(0),
  fNucleons(0), fPhiRot(0), fThetaRot(0), fHe3Counter(-1)
{
  if (fN==0) {
    cout << "Setting up nucleus " << iname << endl;
    Lookup(iname);
  }
}

TGlauNucleus::~TGlauNucleus()
{
  if (fNucleons) {
    delete fNucleons;
  }
  delete fFunction;
  delete fFunction2;
}

void TGlauNucleus::Draw(Double_t xs, Int_t col)
{
  Double_t r = 0.5*TMath::Sqrt(xs/TMath::Pi()/10.);
  TEllipse en;
  en.SetLineColor(col);
  en.SetLineStyle(1);
  en.SetLineWidth(1);

  for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
    TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
    if (!gn->IsSpectator()) {
      en.SetFillColor(0);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    } else {
      en.SetFillColor(col);
      en.SetFillStyle(1001);
      en.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
    }
  }
}

void TGlauNucleus::Lookup(const char* name)
{
  SetName(name);

  TString tmp(name);

  if      (TString(name) == "p")    {fN = 1;   fR = 0.6;      fA = 0;      fW =  0;        fF = 0; fZ=1;}
  else if (TString(name) == "d")    {fN = 2;   fR = 0.01;     fA = 0.5882; fW =  0;        fF = 1; fZ=1;}
  else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;     fA = 0.5882; fW =  0;        fF = 3; fZ=1;}
  else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;     fA = 0.5882; fW =  0;        fF = 4; fZ=1;}
  else if (TString(name) == "He3")  {fN = 3;   fR = 0.01;     fA = 0.5882; fW =  0;        fF = 6; fZ=1;}
  else if (TString(name) == "H3")   {fN = 3;   fR = 0.01;     fA = 0.5882; fW =  0;        fF = 6; fZ=2;}
  else if (TString(name) == "O")    {fN = 16;  fR = 2.608;    fA = 0.513;  fW = -0.051;    fF = 1; fZ=8;}
  else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;     fA = 0.580;  fW = -0.233;    fF = 1; fZ=14;}
  else if (TString(name) == "Si2")  {fN = 28;  fR = 3.34;     fA = 0.580;  fW =  0;        fF = 8; fZ=14; fBeta2=-0.478; fBeta4=0.239;}
  else if (TString(name) == "S")    {fN = 32;  fR = 2.54;     fA = 2.191;  fW =  0.16;     fF = 2; fZ=16;}
  else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766;    fA = 0.586;  fW = -0.161;    fF = 1; fZ=20;}
  else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309;    fA = 0.517;  fW = -0.1308;   fF = 1; fZ=28;}
  else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;      fA = 0.596;  fW =  0;        fF = 1; fZ=29;}
  else if (TString(name) == "Cu2")  {fN = 63;  fR = 4.2;      fA = 0.596;  fW =  0;        fF = 8; fZ=29; fBeta2=0.162; fBeta4=-0.006;}  
  else if (TString(name) == "CuHN") {fN = 63;  fR = 4.28;     fA = 0.5;    fW =  0;        fF = 1; fZ=29;} // from arXiv:0904.4080v1
  else if (TString(name) == "W")    {fN = 186; fR = 6.58;     fA = 0.480;  fW =  0;        fF = 1; fZ=74;}
  else if (TString(name) == "Au")   {fN = 197; fR = 6.38;     fA = 0.535;  fW =  0;        fF = 1; fZ=79;}
  else if (TString(name) == "Au2")  {fN = 197; fR = 6.38;     fA = 0.535;  fW =  0;        fF = 8; fZ=79; fBeta2=-0.131; fBeta4=-0.031;}
  else if (TString(name) == "AuHN") {fN = 197; fR = 6.42;     fA = 0.44;   fW =  0;        fF = 1; fZ=79;} // from arXiv:0904.4080v1
  else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;     fA = 0.546;  fW =  0;        fF = 1; fZ=82;}
  // Uranium description taken from Heinz & Kuhlman, nucl-th/0411054.  In this code, fR is defined as 6.8*0.91, fW=6.8*0.26
  else if (TString(name) == "U")    {fN = 238; fR = 6.188;    fA = 0.54;   fW =  1.77;     fF = 5; fZ=92;}  
  else if (TString(name) == "U2")   {fN = 238; fR = 6.67;     fA = 0.44;   fW =  0;        fF = 8; fZ=92; fBeta2=0.280; fBeta4=0.093;}

  else {
    cout << "Could not find nucleus " << name << endl;
    return;
  }

  switch (fF) {
    case 0: // Proton
      fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
      fFunction->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
      fFunction->SetParameters(fR,fA,fW);
      break;
    case 2: // 3pG
      fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
      fFunction->SetParameters(fR,fA,fW);
      break;
    case 3: // Hulthen
      fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
      fFunction->SetParameters(1/4.38,1/.85);
      break;
    case 4: // Hulthen HIJING
      fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
      fFunction->SetParameters(2/4.38,2/.85);
      break;
    case 5: // Ellipsoid (Uranium)
      fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,15);
      fFunction->SetParameters(fR,fA,0); // same as 3pF but setting W to zero
      break;
    case 6: // He3/H3
      fFunction = 0; // read in file instead
     break;
    case 7: // Deformed nuclei, box method
      fFunction = 0; // no function: only need beta parameters and use uniform box distribution
     break;
    case 8: // Deformed nuclei, TF2 method
      fFunction2 = new TF2("f77","x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",
			   0,20,0.0,TMath::Pi());
      fFunction2->SetNpx(120);
      fFunction2->SetNpy(120);
      fFunction2->SetParameter(0,fR);
      fFunction2->SetParameter(1,fA);
      fFunction2->SetParameter(2,fBeta2);
      fFunction2->SetParameter(3,fBeta4);
      break;
   default:
     cerr << "Could not find function type " << fF << endl;
  }
  return;
}

void TGlauNucleus::SetR(Double_t ir)
{
  fR = ir;
  switch (fF) {
    case 0: // Proton
      fFunction->SetParameter(0,fR);
      break;
    case 1: // 3pF
      fFunction->SetParameter(0,fR);
      break;
    case 2: // 3pG
      fFunction->SetParameter(0,fR);
      break;
  }
}

void TGlauNucleus::SetA(Double_t ia)
{
  fA = ia;
  switch (fF) {
    case 0: // Proton
      break;
    case 1: // 3pF
      fFunction->SetParameter(1,fA);
      break;
    case 2: // 3pG
      fFunction->SetParameter(1,fA);
      break;
  }
}

void TGlauNucleus::SetW(Double_t iw)
{
  fW = iw;
  switch (fF) {
    case 0: // Proton
      break;
    case 1: // 3pF
      fFunction->SetParameter(2,fW);
      break;
    case 2: // 3pG
      fFunction->SetParameter(2,fW);
      break;
   }
}

void TGlauNucleus::ThrowNucleons(Double_t xshift)
{
  if (fNucleons==0) {
    fNucleons=new TObjArray(fN);
    fNucleons->SetOwner();
    for (Int_t i=0;i<fN;i++) {
      TGlauNucleon *nucleon=new TGlauNucleon(); 
      nucleon->SetType(0);
      if (i<fZ) nucleon->SetType(1);
      fNucleons->Add(nucleon); 
    }
  } 
   
  fTrials = 0;
  // just in case we need to rotate whole nucleus
  fPhiRot = gRandom->Rndm()*2*TMath::Pi();
  Double_t cosThetaRot = 2*gRandom->Rndm()-1;
  fThetaRot = TMath::ACos(cosThetaRot);

  Bool_t hulthen = (TString(GetName())=="dh" || TString(GetName())=="dhh");
  Bool_t helium3 = (TString(GetName())=="He3") || (TString(GetName())=="H3");
  if (fN==2 && hulthen) { //special treatmeant for Hulten
    Double_t r = fFunction->GetRandom()/2;
    Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
    Double_t ctheta = 2*gRandom->Rndm() - 1 ;
    Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
     
    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    nucleon1->Reset();
    nucleon1->SetXYZ(r * stheta * TMath::Cos(phi),
		     r * stheta * TMath::Sin(phi),
		     r * ctheta);
    nucleon2->Reset();
    nucleon2->SetXYZ(-nucleon1->GetX(),
		     -nucleon1->GetY(),
		     -nucleon1->GetZ());
    fTrials = 1;
    return;
  } else if (helium3) { 
    if (fHe3Counter == -1) {
      // read in the ascii file into the array and step through the counter
      char filename[100] = "he3_plaintext.dat";
      if ((TString(GetName())=="H3")) {
	sprintf(filename,"h3_plaintext.dat");
      }
      cout << "READING IN THE " << filename << " FILE UPON INITIALIZATION" << endl;
      ifstream myfile;
      myfile.open(filename);
      if (!myfile) {
	cout << "ERROR:  no file for He3/H3 found with name = " << filename << endl;
	gSystem->Exit(123);
      }
      cout << "Reading file for He3/H3 found with name = " << filename << endl;
      Int_t inputcounter = 0;
      while (myfile) {
	if (inputcounter > 6000) break;
	Double_t foo;
	myfile >> fHe3Arr[inputcounter][0][0] >> fHe3Arr[inputcounter][0][1] >> fHe3Arr[inputcounter][0][2]
	       >> fHe3Arr[inputcounter][1][0] >> fHe3Arr[inputcounter][1][1] >> fHe3Arr[inputcounter][1][2]
	       >> fHe3Arr[inputcounter][2][0] >> fHe3Arr[inputcounter][2][1] >> fHe3Arr[inputcounter][2][2]
	       >> foo >> foo >> foo >> foo;
	inputcounter++;
      }  
      myfile.close();
      fHe3Counter=0;
    } // done reading in the file the first time

    if (fHe3Counter > 5999) 
      fHe3Counter = 0;
    TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
    TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
    TGlauNucleon *nucleon3=(TGlauNucleon*)(fNucleons->At(2));
    nucleon1->Reset();
    nucleon1->SetXYZ(fHe3Arr[fHe3Counter][0][0],
		     fHe3Arr[fHe3Counter][0][1],
		     fHe3Arr[fHe3Counter][0][2]);
    nucleon1->RotateXYZ(fPhiRot,fThetaRot);
    nucleon2->Reset();
    nucleon2->SetXYZ(fHe3Arr[fHe3Counter][1][0],
		     fHe3Arr[fHe3Counter][1][1],
		     fHe3Arr[fHe3Counter][1][2]);
    nucleon2->RotateXYZ(fPhiRot,fThetaRot);
    nucleon3->Reset();
    nucleon3->SetXYZ(fHe3Arr[fHe3Counter][2][0],
		     fHe3Arr[fHe3Counter][2][1],
		     fHe3Arr[fHe3Counter][2][2]);
    nucleon3->RotateXYZ(fPhiRot,fThetaRot);
    fHe3Counter++;
    fTrials = 1;
  } else {
    for (Int_t i = 0; i<fN; i++) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while (1) {
	fTrials++;
	Bool_t nucleon_inside = 0;
	Double_t x=999;
	Double_t y=999;
	Double_t z=999;
	
	if (fF==5||fF==7) { // the extended way, throw in a box and test the weight
	  while (!nucleon_inside) {
	    x = (fR*2)*(gRandom->Rndm() * 2 - 1);
	    y = (fR*2)*(gRandom->Rndm() * 2 - 1);
	    z = (fR*2)*(gRandom->Rndm() * 2 - 1);
	    Double_t r = TMath::Sqrt(x*x+y*y);
	    Double_t theta = TMath::ATan2(r,z);
	    Double_t R = TMath::Sqrt(x*x+y*y+z*z);

	    Double_t Rtheta = fR;
	    if (fF==5)
	      Rtheta= fR + fW*TMath::Cos(theta)*TMath::Cos(theta);
	    if (fF==7)
	      Rtheta = fR*(1+fBeta2*ROOT::Math::sph_legendre(2,0,theta)+fBeta4*ROOT::Math::sph_legendre(4,0,theta));
	      
	    Double_t prob = 1/(1+TMath::Exp((R-Rtheta)/fA));
	    if (gRandom->Rndm()<prob) nucleon_inside=1;
	  }
	} else if (fF==8) { // use TF2
	  Double_t r;
	  Double_t theta;
	  fFunction2->GetRandom2(r,theta);
	  Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
	  x = r * TMath::Sin(phi) * TMath::Sin(theta);
	  y = r * TMath::Cos(phi) * TMath::Sin(theta);
	  z = r *                   TMath::Cos(theta);
	} else  { 
	  // the normal way
	  Double_t r = fFunction->GetRandom();
	  Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
	  Double_t ctheta = 2*gRandom->Rndm() - 1 ;
	  Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
	  x = r * stheta * TMath::Cos(phi);
	  y = r * stheta * TMath::Sin(phi);      
	  z = r * ctheta;      
	}
	
	nucleon->SetXYZ(x,y,z);
	if (fF==5) 
	  nucleon->RotateXYZ(fPhiRot,fThetaRot); // Uranium etc.

	if (fMinDist<0) break;
	Bool_t test=1;
	for (Int_t j = 0; j<i; j++) {
	  TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
	  Double_t xo=other->GetX();
	  Double_t yo=other->GetY();
	  Double_t zo=other->GetZ();
	  Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
				      (y-yo)*(y-yo)+
				      (z-zo)*(z-zo));
	  if (dist<fMinDist) {
	    test=0;
	    break;
	  }
	}
	if (test) break; //found nucleuon outside of mindist
      }
    }
  }
      
  Double_t sumx=0;       
  Double_t sumy=0;       
  Double_t sumz=0;       
  if (1) { // set the centre-of-mass to be at zero (+xshift)
    for (Int_t i = 0; i<fN; i++) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
    }
    sumx = sumx/fN;  
    sumy = sumy/fN;  
    sumz = sumz/fN;  
  }
  for (Int_t i = 0; i<fN; i++) {
    TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
    nucleon->SetXYZ(nucleon->GetX()-sumx + xshift,
		    nucleon->GetY()-sumy,
		    nucleon->GetZ()-sumz);
  }
}

//---------------------------------------------------------------------------------
ClassImp(TGlauberMC)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(const char* NA, const char* NB, Double_t xsect, Double_t xsectsigma) :
  fANucleus(NA),fBNucleus(NB),
  fXSect(0),fXSectOmega(0),fXSectLambda(0),fXSectEvent(0),
  fNucleonsA(0),fNucleonsB(0),fNucleons(0),
  fAN(0),fBN(0),fNt(0),
  fMeanX2(0),fMeanY2(0),fMeanXY(0),fMeanXParts(0),
  fMeanYParts(0),fMeanXSystem(0),fMeanYSystem(0),  
  fMeanX_A(0),fMeanY_A(0),fMeanX_B(0),fMeanY_B(0),fB_MC(0),
  fEvents(0),fTotalEvents(0),fBMin(0),fBMax(0),fMaxNpartFound(0),
  fNpart(0),fNpartA(0),fNpartB(0),fSumW(0),fSumWA(0),fSumWB(0), 
  fNcoll(0),fSx2(0),fSy2(0),fSxy(0),f2Cx(0),fPTot(0)
{
  fBMin = 0;
  fBMax = 20;
  fXSect = xsect;
  if (xsectsigma>0) {
    fXSectOmega = xsectsigma;
    fXSectLambda = 1;
    fPTot = new TF1("fPTot","((x/[2])/(x/[2]+[0]))*exp(-(((x/[2])/[0]-1 )**2)/([1]*[1]))/[2]",0,300);
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    fPTot->SetNpx(1000);
    fXSectLambda = fXSect/fPTot->GetHistogram()->GetMean();
    cout << "final lambda=" << fXSectLambda << endl;
    fPTot->SetParameters(fXSect,fXSectOmega,fXSectLambda);
    cout << "final <sigma>=" << fPTot->GetHistogram()->GetMean() << endl;
  }

  TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
  SetName(name);
  SetTitle(title);
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
  // prepare event
  fANucleus.ThrowNucleons(-bgen/2.);
  fNucleonsA = fANucleus.GetNucleons();
  fAN = fANucleus.GetN();
  for (Int_t i = 0; i<fAN; i++) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    nucleonA->SetInNucleusA();
  }
  fBNucleus.ThrowNucleons(bgen/2.);
  fNucleonsB = fBNucleus.GetNucleons();
  fBN = fBNucleus.GetN();
  for (Int_t i = 0; i<fBN; i++) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    nucleonB->SetInNucleusB();
  }
  
  // "ball" diameter = distance at which two balls interact
  if (fPTot)
    fXSectEvent = fPTot->GetRandom();
  else 
    fXSectEvent = fXSect;

  // "ball" diameter = distance at which two balls interact
  Double_t d2 = (Double_t)fXSectEvent/(TMath::Pi()*10); // in fm^2

  // for each of the A nucleons in nucleus B
  for (Int_t i = 0; i<fBN; i++) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    for (Int_t j = 0 ; j < fAN ;j++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
      Double_t dx = nucleonB->GetX()-nucleonA->GetX();
      Double_t dy = nucleonB->GetY()-nucleonA->GetY();
      Double_t dij = dx*dx+dy*dy;
      if (dij < d2) {
	nucleonB->Collide();
	nucleonA->Collide();
      }
    }
  }
  return CalcResults(bgen);
}

Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  // calc results for the given event
  fNpart=0;
  fNpartA=0;
  fNpartB=0;
  fNcoll=0;
  fMeanX2=0;
  fMeanY2=0;
  fMeanXY=0;
  fMeanXParts=0;
  fMeanYParts=0;
  fMeanXSystem=0;
  fMeanYSystem=0;
  fMeanX_A=0;
  fMeanY_A=0;
  fMeanX_B=0;
  fMeanY_B=0;
  fSumW=0;
  fSumWA=0;
  fSumWB=0;

  Double_t sinphi[10];
  Double_t cosphi[10];
  Double_t rn[10];
  for (Int_t i = 0;i<10;i++) {
    fEccN[i] = 0;
    fPsiN[i] = 0;
    sinphi[i] = 0;
    cosphi[i] = 0;
    rn[i] = 0;
  }

  for (Int_t i = 0; i<fAN; i++) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    fMeanXSystem  += xA;
    fMeanYSystem  += yA;
    fMeanX_A  += xA;
    fMeanY_A  += yA;
    if (nucleonA->IsWounded()) {
      Double_t w = nucleonA->Get2CWeight(f2Cx);
      fSumW  += w;
      fSumWA += w;
      fNpart++;
      fNpartA++;
      fMeanXParts  += xA * w;
      fMeanYParts  += yA * w;
      fMeanX2 += xA * xA * w;
      fMeanY2 += yA * yA * w;
      fMeanXY += xA * yA * w;
    }
  }

  for (Int_t i = 0; i<fBN; i++) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    fMeanXSystem  += xB;
    fMeanYSystem  += yB;
    fMeanX_B  += xB;
    fMeanY_B  += yB;
    if (nucleonB->IsWounded()) {
      Double_t w = nucleonB->Get2CWeight(f2Cx);
      fNpart++;
      fNpartB++;
      fSumW += w;
      fSumWB += w;
      fMeanXParts  += xB * w;
      fMeanYParts  += yB * w;
      fMeanX2 += xB * xB * w;
      fMeanY2 += yB * yB * w;
      fMeanXY += xB * yB * w;
      fNcoll += nucleonB->GetNColl();
    }
  }
  if (fNpart>0) {
    fMeanXParts /= fSumW;
    fMeanYParts /= fSumW;
    fMeanX2 /= fSumW;
    fMeanY2 /= fSumW;
    fMeanXY /= fSumW;
  } else {
    fMeanXParts = 0;
    fMeanYParts = 0;
    fMeanX2 = 0;
    fMeanY2 = 0;
    fMeanXY = 0;
  }
   
  if (fAN+fBN>0) {
    fMeanXSystem /= (fAN + fBN);
    fMeanYSystem /= (fAN + fBN);
  } else {
    fMeanXSystem = 0;
    fMeanYSystem = 0;
  }
  if (fAN>0) {
    fMeanX_A /= fAN;
    fMeanY_A /= fAN;
  } else {
    fMeanX_A = 0;
    fMeanY_A = 0;
  }

  if (fBN>0) {
    fMeanX_B /= fBN;
    fMeanY_B /= fBN;
  } else {
    fMeanX_B = 0;
    fMeanY_B = 0;
  }
  
  fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
  fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
  fSxy=fMeanXY-fMeanXParts*fMeanYParts;

  if (fNpart>0) {
    // do full moments relative to meanX and meanY
    for (Int_t n = 1;n<6;n++) {
      for (Int_t ia = 0;ia<fAN;ia++) {
	TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(ia));
	if (nucleonA->IsWounded()) {
	  Double_t xA=nucleonA->GetX() - fMeanXParts;
	  Double_t yA=nucleonA->GetY() - fMeanYParts;
	  Double_t r = TMath::Sqrt(xA*xA+yA*yA);
	  Double_t phi = TMath::ATan2(yA,xA);
	  cosphi[n] += TMath::Power(r,n)*TMath::Cos(n*phi);
	  sinphi[n] += TMath::Power(r,n)*TMath::Sin(n*phi);
	  rn[n] += TMath::Power(r,n);
	}
      }
      for (Int_t ib = 0;ib<fBN;ib++) {
	TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(ib));
	if (nucleonB->IsWounded()) {
	  Double_t xB=nucleonB->GetX() - fMeanXParts;
	  Double_t yB=nucleonB->GetY() - fMeanYParts;
	  Double_t r = TMath::Sqrt(xB*xB+yB*yB);
	  Double_t phi = TMath::ATan2(yB,xB);
	  cosphi[n] += TMath::Power(r,n)*TMath::Cos(n*phi);
	  sinphi[n] += TMath::Power(r,n)*TMath::Sin(n*phi);
	  rn[n] += TMath::Power(r,n);
	}
      }
      cosphi[n] /= fNpart;
      sinphi[n] /= fNpart;
      rn[n] /= fNpart;
      fPsiN[n] = (TMath::ATan2(sinphi[n],cosphi[n]) + TMath::Pi())/n;
      fEccN[n] = TMath::Sqrt(sinphi[n]*sinphi[n]+cosphi[n]*cosphi[n])/rn[n];
    }
  }
  fB_MC = bgen;
  fTotalEvents++;
  if (fNpart>0) fEvents++;
  if (fNpart==0) return kFALSE;
  if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;
  return kTRUE;
}

void TGlauberMC::Draw(Option_t* /*option*/)
{
  TEllipse e;
  e.SetFillColor(0);
  e.SetFillStyle(0);
  e.SetLineColor(1);
  e.SetLineStyle(2);
  e.SetLineWidth(1);
  e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
  e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
  fANucleus.Draw(fXSect, 3);
  fBNucleus.Draw(fXSect, 7);
  
  Double_t sy2 = GetSy2();
  Double_t sx2 = GetSx2();
  Double_t d;
  Double_t phase = 0;
  if (sy2<sx2) {
    d = sx2;
    sx2 = sy2;
    sy2 = d;
    phase = 3.14/2.;
  }

  Double_t x1 = (0.5*(sy2-sx2)+TMath::Sqrt(pow((sy2-sx2),2.)-4*pow(GetSxy(),2)));
  Double_t ang = TMath::ATan2(-GetSxy(),x1)+phase;
  TLine l;
  l.DrawLine(-10*TMath::Cos(ang),-10*TMath::Sin(ang),10*TMath::Cos(ang),10*TMath::Sin(ang));
}

Double_t TGlauberMC::GetTotXSect() const
{
  return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
  return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
    TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
  if (!fNucleonsA || !fNucleonsB) return 0;
  if (fNucleons) return fNucleons;

  fNucleonsA->SetOwner(0);
  fNucleonsB->SetOwner(0);
  TObjArray *allnucleons=new TObjArray(fAN+fBN);
  allnucleons->SetOwner();
  for (Int_t i = 0; i<fAN; i++) {
    allnucleons->Add(fNucleonsA->At(i));
  }
  for (Int_t i = 0; i<fBN; i++) {
    allnucleons->Add(fNucleonsB->At(i));
  }
  fNucleons = allnucleons;
  return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  if (bgen<0) 
    bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);

  return CalcEvent(bgen);
}

void TGlauberMC::Run(Int_t nevents, Double_t b)
{
  TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
  TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
  if (fNt == 0) {
    fNt = new TNtuple(name,title,
		      "Npart:Ncoll:B"
                      ":MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB"
                      ":NpartA:NpartB:PhiA:ThetaA:PhiB:ThetaB:Psi1:Ecc1:Psi2:Ecc2:Psi3:Ecc3:Psi4:Ecc4:Psi5:Ecc5");
    fNt->SetDirectory(0);
  }

  for (Int_t i = 0;i<nevents;i++) {
    while (!NextEvent(b)) {}

    Float_t v[33];
    v[0]  = GetNpart();
    v[1]  = GetNcoll();
    v[2]  = fB_MC;
    v[3]  = fMeanXParts;
    v[4]  = fMeanYParts;
    v[5]  = fMeanX2;
    v[6]  = fMeanY2;
    v[7]  = fMeanXY;
    v[8]  = fSx2;
    v[9]  = fSy2;
    v[10] = fSxy;
    v[11] = fMeanXSystem;
    v[12] = fMeanYSystem;
    v[13] = fMeanX_A;
    v[14] = fMeanY_A;
    v[15] = fMeanX_B;
    v[16] = fMeanY_B;
    v[17] = fNpartA;
    v[18] = fNpartB;
    v[19] = fANucleus.GetPhiRot();
    v[20] = fANucleus.GetThetaRot();
    v[21] = fBNucleus.GetPhiRot();
    v[22] = fBNucleus.GetThetaRot();
    v[23] = fPsiN[1];
    v[24] = fEccN[1];
    v[25] = fPsiN[2];
    v[26] = fEccN[2];
    v[27] = fPsiN[3];
    v[28] = fEccN[3];
    v[29] = fPsiN[4];
    v[30] = fEccN[4];
    v[31] = fPsiN[5];
    v[32] = fEccN[5];

    fNt->Fill(v);

    if (!(i%50)) 
      cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;
  }
  cout << endl << "Done!" << endl;
}
#endif
