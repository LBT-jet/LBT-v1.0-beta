//....Generate proton + proton events within PYTHIA8, 
//....then construct a leading jet for each events using fastjet
#include <iostream>    //....C++ headers
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include <cstdio>    //....C headers
#include <cstdlib>
#include <cmath>
#include <ctime>

//#include "fastjet/ClusterSequence.hh"    //....FASTJET headers
#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support

#include "TH1.h"    //....ROOT headers
#include "TH2.h"
#include "TFile.h"

using namespace std;
using namespace fastjet;

//....Classes definition
//....FASTJET class dealing with negative particles
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
      } else { pab.set_user_index(0);}

    }

  private:
    const int _ui;  
};

int main(int argc, char* argv[]) {
  //....time counting begins
  struct tm *local_start;
  time_t time_start;
  time_start=time(NULL);
  local_start=localtime(&time_start);

  char buf1[80];
  strftime(buf1,80,"Current Time: %Y-%m-%d %H:%M:%S",local_start);
  cout << "the program starts at:" <<endl;
  cout << buf1 << endl;

  //....Generator. Initialization. 
  char charName[256];

  sprintf(charName, "./pp.txt");
  ifstream ppdata(charName);
  cout << charName << endl;

  sprintf(charName, "./positive.txt");
  ifstream posdata(charName);
  cout << charName << endl;

  sprintf(charName, "./negative.txt");
  ifstream negdata(charName);
  cout << charName << endl;

  if (!ppdata || !posdata  || !negdata) {
    cout << "no input data!" << endl;
    exit(1);
  }

  const int nE = 10;
  //const int nR = 9, ny = 7;
  //double R[nR] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  //double y[ny] = {0.0, 0.3, 0.8, 1.2, 1.6, 2.1, 2.8};
  const int nR = 1, ny = 2;
  double R[nR] = {0.4};
  double y[ny] = {0, 2.8};
  sprintf(charName, "dNdpT.root");
  TFile *rootFile = new TFile(charName, "Recreate");
  TH1D *hdNdpT_pp[nR][ny];
  TH1D *hdNdpT_AA[nR][ny];
  TH1D *hpTloss[nR][ny];

  const int nptbins = 1400;
  const double ptmin = 0.0, ptmax = 1400.0;

  for (int i = 0; i < nR; ++i){
    for (int j = 0; j < ny; ++j){
      //.. dNdpT
      sprintf(charName, "dNdpT_pp_R%d_y%d", i, j);
      hdNdpT_pp[i][j] = new TH1D(charName, "", nptbins, ptmin, ptmax);

      sprintf(charName, "dNdpT_AA_R%d_y%d", i, j);
      hdNdpT_AA[i][j] = new TH1D(charName, "", nptbins, ptmin, ptmax);

      sprintf(charName, "pTloss_R%d_y%d", i, j);
      hpTloss[i][j] = new TH1D(charName, "", nptbins, ptmin, ptmax);
    }
  }


  double v2 = 0, v3 = 0, Psi2 = 0, Psi3 = 0;
//  cout << Psi2 << "  " << Psi3 << "  " << v2 << "  " << v3 << endl;


  for (int i = 0; i < nE; i++) {
    cout << "event: " << i << endl;
    vector<fastjet::PseudoJet> input_particles_pp;
    int id_pp;
    double px_pp, py_pp, pz_pp, e_pp, timeplus_pp;
    string buf_pp;

    ppdata >> buf_pp;
    while (buf_pp != "#") {
      id_pp = atoi(buf_pp.c_str());
      ppdata >> e_pp >> px_pp >> py_pp >> pz_pp >> timeplus_pp;
      ppdata >> buf_pp;

      fastjet::PseudoJet pj_pp;
      pj_pp.reset_momentum(px_pp, py_pp, pz_pp, e_pp);
      pj_pp.set_user_index(fabs(id_pp));
      input_particles_pp.push_back(pj_pp); 
    }
    int iE_pp, np_pp;
    ppdata >> iE_pp >> np_pp;

//    cout << "pp data ready" << endl;

    vector<fastjet::PseudoJet> input_particles;
    int id, icat;
    double px, py, pz, e, xx, xy, xz, xt;
    string buf;

    posdata >> buf;
    while (buf != "#") {
      id = atoi(buf.c_str());
      posdata >> px >> py >> pz >> e >> icat; 
      //posdata >> px >> py >> pz >> e >> xx >> xy >> xz >> xt >> icat; 
      posdata >> buf;

      fastjet::PseudoJet pj;
      pj.reset_momentum(px, py, pz, e);
      pj.set_user_index(fabs(id));
      input_particles.push_back(pj); 
    }
    int iE, np;
    posdata >> iE >> np;
//    cout << "positive data ready" << endl;

    negdata >> buf;
    while (buf != "#") {
      id = atoi(buf.c_str());
      negdata >> px >> py >> pz >> e; 
      negdata >> buf;

      fastjet::PseudoJet pj;
      pj.reset_momentum(-px, -py, -pz, -e);
      pj.set_user_index(-fabs(id));
      input_particles.push_back(pj); 
    }
    negdata >> iE >> np;
//    cout << "negative data ready" << endl;

    for (int iR = 0; iR < nR; iR++) {
//      cout << iR << endl;
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R[iR]);
      NegativeEnergyRecombiner uir(-123456);
      jet_def.set_recombiner(&uir);

      //            fastjet::ClusterSequence clust_seq(input_particles, jet_def);

      double maxrap = 5.0;
      int n_repeat = 1; // default is 1
      double ghost_area = 0.01; // this is the default
      fastjet::GhostedAreaSpec area_spec(maxrap, n_repeat, ghost_area);

      fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

      //.... pp fastjet
      // The only change is the usage of a ClusterSequenceArea rather thana ClusterSequence
//      cout << "pp UES" << endl;
      fastjet::ClusterSequenceArea clust_seq_pp(input_particles_pp, jet_def, area_def);

      //....Get the resulting jets ordered in pt
      vector<fastjet::PseudoJet> jets_pp = sorted_by_pt(clust_seq_pp.inclusive_jets(0.0));

      if (jets_pp.size() == 0 ) continue;

      //.. UE subtraction: precedure

      vector<double> jpt_pp;
      vector<int> triggerjet_pp;
      for(int ij = 0; ij < jets_pp.size() ; ++ij ) {
        jpt_pp.push_back(jets_pp[ij].mt());
        triggerjet_pp.push_back(0);
      }

      for(int ij = 0; ij < jets_pp.size() ; ++ij ) {
        vector<PseudoJet> csts = jets_pp[ij].constituents();
        bool trigger1 = false;  //.. for at least one tower mt > 3.0
        bool trigger2 = false;  //.. mtmax / mtavg > 4.0
        for (int ic = 0; ic < csts.size(); ++ic) {
          if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
            if (csts[ic].mt() > 3.0) trigger1 = true;
          }
        }

        if (trigger1) {
          int ncsts = 0;
          double avgmt = 0.0;
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              ncsts += 1;
              avgmt += csts[ic].mt();
            }
          }
          avgmt /= ncsts;

          double maxmt;
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              maxmt = csts[ic].mt();
            }
          }

          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              if (csts[ic].mt() > maxmt) maxmt = csts[ic].mt();
            }
          }

          if (maxmt / avgmt > 4) trigger2 = true;
          if (trigger2)  {
            triggerjet_pp[ij] = 1;
          }
        }
      }

      double areatotal_pp = 0.0;
      double uemt_pp = 0.0;
      for(int ij = 0; ij < jets_pp.size() ; ++ij ) {
        if (triggerjet_pp[ij]) {
          areatotal_pp += jets_pp[ij].area();
        }
        else {
          vector<PseudoJet> csts = jets_pp[ij].constituents();
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              if (csts[ic].user_index() < 0) uemt_pp -= csts[ic].mt();
              else uemt_pp += csts[ic].mt();
            }
          }
        }
      }

      uemt_pp /= ((2 * 3.14) * 2 * maxrap - areatotal_pp);


      for(int ij = 0; ij < jets_pp.size() ; ++ij ) {
        if (triggerjet_pp[ij]) {
          double jetphi = jets_pp[ij].phi_std();
          double dphiJet2 = jetphi - 0;
          if(dphiJet2 > M_PI) dphiJet2 = -2*M_PI + dphiJet2;
          if(dphiJet2 < -M_PI) dphiJet2 = +2*M_PI + dphiJet2;
          jpt_pp[ij] -= uemt_pp * jets_pp[ij].area();
        }
      }
//      cout << "pp UES ready at event " << i << " " << "at jet " << jets_pp.size() << " R:" << iR  << endl;

      //.... AA fastjet
      // The only change is the usage of a ClusterSequenceArea rather thana ClusterSequence
      fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def, area_def);

      //....Get the resulting jets ordered in pt
      vector<fastjet::PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(0.0));

      if (jets.size() == 0 ) continue;

      //.. UE subtraction: precedure

      vector<double> jpt;
      vector<int> triggerjet;
      for(int ij = 0; ij < jets.size() ; ++ij ) {
        jpt.push_back(jets[ij].mt());
        triggerjet.push_back(0);
      }

      for(int ij = 0; ij < jets.size() ; ++ij ) {
        vector<PseudoJet> csts = jets[ij].constituents();
        bool trigger1 = false;  //.. for at least one tower mt > 3.0
        bool trigger2 = false;  //.. mtmax / mtavg > 4.0
        for (int ic = 0; ic < csts.size(); ++ic) {
          if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
            if (csts[ic].mt() > 3.0) trigger1 = true;
          }
        }

        if (trigger1) {
          int ncsts = 0;
          double avgmt = 0.0;
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              ncsts += 1;
              avgmt += csts[ic].mt();
            }
          }
          avgmt /= ncsts;

          double maxmt;
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              maxmt = csts[ic].mt();
            }
          }

          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              if (csts[ic].mt() > maxmt) maxmt = csts[ic].mt();
            }
          }

          if (maxmt / avgmt > 4) trigger2 = true;
          if (trigger2)  {
            triggerjet[ij] = 1;
          }
        }
      }

      double areatotal = 0.0;
      double uemt = 0.0;
      for(int ij = 0; ij < jets.size() ; ++ij ) {
        if (triggerjet[ij]) {
          areatotal += jets[ij].area();
        }
        else {
          vector<PseudoJet> csts = jets[ij].constituents();
          for (int ic = 0; ic < csts.size(); ++ic) {
            if (csts[ic].rap() >= -maxrap && csts[ic].rap() < maxrap) {
              if (csts[ic].user_index() < 0) uemt -= csts[ic].mt();
              else uemt += csts[ic].mt();
            }
          }
        }
      }

      uemt /= ((2 * 3.14) * 2 * maxrap - areatotal);


      for(int ij = 0; ij < jets.size() ; ++ij ) {
        if (triggerjet[ij]) {
          double jetphi = jets[ij].phi_std();
          double dphiJet2 = jetphi - Psi2;
          if(dphiJet2 > M_PI) dphiJet2 = -2*M_PI + dphiJet2;
          if(dphiJet2 < -M_PI) dphiJet2 = +2*M_PI + dphiJet2;
          jpt[ij] -= uemt * jets[ij].area() * (1 + 2*v2*cos(2*dphiJet2));
        }
      }

//      cout << "AA UES ready at event " << i << " " << "at jet " << jets.size() << " R:" << iR  << endl;

      //.. pp inclusive jet
      for(int ij = 0; ij < jets_pp.size(); ++ij ) {
        //.. pp jets
        double jetrap_pp = fabs(jets_pp[ij].rap());
        double jetmt_pp = jpt_pp[ij];
//        cout << "pp jet:" << ij << " " << jetrap_pp << "  " << jetmt_pp << endl;

        //.. pp dNdpT

        for (int iy = 0; iy < ny; ++iy) {
          if (iy == ny - 1 ) {
            if (jetrap_pp < y[iy]) {
              hdNdpT_pp[iR][iy]->Fill(jetmt_pp, 1.0 / 2.0 / y[iy]);
            } 
          }
          else {
            if (jetrap_pp >= y[iy] && jetrap_pp < y[iy + 1]) {
              hdNdpT_pp[iR][iy]->Fill(jetmt_pp, 1.0 / 2.0 / (y[iy + 1] - y[iy]));
            }
          }
        }
      }

      //.. AA inclusive jet
      for(int ij = 0; ij < jets.size(); ++ij ) {
        //.. AA jets
        double jetrap = fabs(jets[ij].rap());
        double jetmt = jpt[ij];
//        cout << "AA jet:" << ij << " " << jetrap << "  " << jetmt << endl;

        //.. AA dNdpT

        for (int iy = 0; iy < ny; ++iy) {
          if (iy == ny - 1 ) {
            if (jetrap < y[iy]) {
              hdNdpT_AA[iR][iy]->Fill(jetmt, 1.0 / 2.0 / y[iy]);
            } 
          }
          else {
            if (jetrap >= y[iy] && jetrap < y[iy + 1]) {
              hdNdpT_AA[iR][iy]->Fill(jetmt, 1.0 / 2.0 / (y[iy + 1] - y[iy]));
            }
          }
        }
      }

      // .. leading jet pT loss
      if (jets_pp.size() > 0 && jets.size() > 0) {
        double jetrap_pp = fabs(jets_pp[0].rap());
        double jetmt_pp = jpt_pp[0];
        double jetrap = fabs(jets[0].rap());
        double jetmt = jpt[0];
        for (int iy = 0; iy < ny; ++iy) {
          if (iy == ny - 1 ) {
            if (jetrap_pp < y[iy]) {
              hpTloss[iR][iy]->Fill(jetmt_pp, jetmt_pp - jetmt);
            } 
          }
          else {
            if (jetrap_pp >= y[iy] && jetrap_pp < y[iy + 1]) {
              hpTloss[iR][iy]->Fill(jetmt_pp, jetmt_pp - jetmt);
            }
          }
        }
      }
    }
  }

  for (int iR = 0; iR < nR; ++iR) {
    for (int iy = 0; iy < ny; ++iy) {
      hdNdpT_pp[iR][iy]->Scale(1.0/nE);
      hdNdpT_AA[iR][iy]->Scale(1.0/nE);
      hpTloss[iR][iy]->Scale(1.0/nE);
    }
  }

  rootFile->Write();
  rootFile->Close();

  //....time counting ends
  struct tm *local_end;
  time_t time_end;
  time_end=time(NULL);
  local_end=localtime(&time_end);

  char buf2[80];
  strftime(buf2,80,"Current Time: %Y-%m-%d %H:%M:%S",local_end);
  cout << "the program ends at:" << endl;
  cout << buf2 << endl;

  int cost,nh,nm,ns;
  cost=difftime(time_end,time_start);

  nh=cost/3600;
  nm=(cost%3600)/60;
  ns=(cost%3600)%60;

  cout << "the program costs:" << endl;
  cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;

  return 0;
}
