#include"LBT.h"

//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
//////////////////////////////////////////

int main(int argc, char **argv) {	  				
  //.... The program starts
  struct tm *local_start;
  time_t time_start;
  time_start = time(NULL);
  local_start = localtime(&time_start);

  char buf1[80];
  strftime(buf1, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_start);
  cout << "the program starts at:" << endl;
  cout << buf1 << endl;
  //.... Time counts


  char charName[1024];

  //.... Output files in pp
  sprintf(charName, "pp.txt");
  ofstream ppdata(charName);
  sprintf(charName, "cs.txt");
  ofstream cs(charName);


  //.... Output files in AA
  sprintf(charName, "positive.txt");
  ofstream positiveLBT(charName);
  sprintf(charName, "negative.txt");
  ofstream negativeLBT(charName);


  //.... PYTHIA 8 Generator. Initialization. 
  Pythia pythia;
  Event& event = pythia.event;
  sprintf(charName, "m01.cmnd");
  pythia.readFile(charName);
  int nEvent = pythia.mode("Main:numberOfEvents");
  pythia.init();


  //.... LBTclass. Initialization
  LBTclass *LBT = new LBTclass();
  LBT->LBTinitialize();

  //.... LBT input files
  char geo_route[1024];
  sprintf(geo_route,"%s/geometry.dat",LBT->bulk3D_route);   

  ifstream f_geo(geo_route);
  cout << geo_route << endl;
  if (!f_geo.is_open())  {
    cout << "No geometry!" << endl;
    exit(1);
  }

    //.... Geometry profile
    const int numnucleon = 10000;
    double numjet[numnucleon] = {0.0};
    double Xnucleon[numnucleon] = {0.0};
    double Ynucleon[numnucleon] = {0.0};
    double numjettotal = 0, randomxy = 0;

//    cout << "geometry!" << endl;
    for (int i = 0; i < numnucleon; ++i) {
      f_geo >> numjet[i] >> Xnucleon[i] >> Ynucleon[i];
//      cout << numjet[i] << "   " << Xnucleon[i] << "   " << Ynucleon[i] << endl;
      numjettotal += numjet[i];
    }

  double sigma = 0.0;
  double sigmaErr = 0.0;
  int numEvent=0;


  //.... Event loop
  for(int n=1; ; ++n) {		
    if (!pythia.next()) continue;
    if (numEvent == nEvent) break;

    vector<int> id;
    vector<double> e, px, py, pz, timePlus;

    for (int i = 0; i < event.size(); ++i) {

      if (event[i].isFinal() && fabs(event[i].eta()) < 5.0) {
        if (fabs(event[i].id()) == 1 || fabs(event[i].id()) == 2 ||  fabs(event[i].id()) == 3 || event[i].id() == 21) {

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
                        if (abs(event[IDiii].status())==23 || abs(event[IDiii].status())==21 || abs(event[IDiii].status())==12) {
                            ishower=42;
                            timebreaker=1;			
                        }

                        IDmom1=event[IDiii].mother1();
                        IDmom2=event[IDiii].mother2();

                        if (IDmom1==IDmom2 && IDmom1==0) {
                            //cout<<"IDmom1==IDmom2 && IDmom1==0 timebreak"<<endl;
                            ishower=40;
                            timebreaker=1;
                        }//if(IDmom1==IDmom2 && IDmom1==0)

                        if(IDmom1==IDmom2 && IDmom1>0) {
                            IDmom0=IDmom1;
                        }//if(IDmom1==IDmom2 && IDmom1>0)

                        if(IDmom1>0 && IDmom2==0) {
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

                                //double x_split=event[IDiii].e()/event[IDmom0].e();

                                if (x_split>1) x_split=1.0/x_split;

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
                                    if(kt_daughter > 0.0001) {
                                        //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                        timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    }
                                }
                            }
                        }//if(IDmom1>0 && IDmom2==0)

                        if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0) {
                            if(event[IDmom1].e()>event[IDmom2].e()) {
                                IDmom0=IDmom1;			
                            }
                            if(event[IDmom1].e()<=event[IDmom2].e()) {
                                IDmom0=IDmom2;			
                            }

                            double IDdaughter1=event[IDmom0].daughter1();
                            double IDdaughter2=event[IDmom0].daughter2();

                            if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){

                                p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
                                p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
                                p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
                                p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();

                                double x_split=event[IDiii].e()/p4[0];

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

                                if(x_split<0.5) {
                                    //if(kt_daughter > 0.0001 && Q2 > 1.0)
                                    if(kt_daughter > 0.0001) {
                                        //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                        timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                    }
                                }
                            }
                        }//if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0)

                        if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0) {
                            if(event[IDmom1].e()>event[IDmom2].e()) {
                                IDmom0=IDmom1;			
                            }
                            if(event[IDmom1].e()<=event[IDmom2].e()) {
                                IDmom0=IDmom2;			
                            }			

                            double IDdaughter1=event[IDmom0].daughter1();
                            double IDdaughter2=event[IDmom0].daughter2();

                            if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){

                                p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
                                p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
                                p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
                                p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();

                                double x_split=event[IDiii].e()/p4[0];

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

                                if(x_split<0.5) {
                                    //if(kt_daughter > 0.0001 && Q2 > 1.0)
                                    if(kt_daughter > 0.0001) {						
                                        //timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                        timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;					
                                    }
                                }
                            }
                        }//if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0)		


                        if(abs(event[IDmom0].status())==23 || abs(event[IDmom0].status())==21 || abs(event[IDmom0].status())==12) {

                            // cout<<"event[IDmom0].status()"<<" "<<event[IDmom0].status()<<endl;

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
                    
          id.push_back(event[i].id());
          e.push_back(event[i].e());
          px.push_back(event[i].px());
          py.push_back(event[i].py());
          pz.push_back(event[i].pz());
          timePlus.push_back(timeplus);
          ppdata << event[i].id() << "  " << event[i].e() << "  " << event[i].px() << "  " << event[i].py() << "  " << event[i].pz() << " " << timeplus << endl;
        }
      }
    }
    ppdata << "#  " << numEvent << "  " << id.size() << endl;

    //pythia.stat();
    sigma += pythia.info.sigmaGen();
    sigmaErr += pythia.info.sigmaErr();


    if(numEvent == LBT->ncall) break;


    //.... Clear particle list in LBT
    LBT->LBTclear();  


    double R1 = 0, XXX = 0, YYY = 0;
    randomxy = LBT->ran0(&LBT->NUM1);
    for (int i = 0; i < numnucleon; ++i) {
      R1 += numjet[i] / numjettotal;
      if (randomxy < R1) {
        XXX = Xnucleon[i];
        YYY = Ynucleon[i];
        break;
      }
    }


    //.... Jet shower partons
    int nj0 = 0;
    for ( int i = 0; i < id.size(); i++) {
      int fjetshower = id[i];
      double p0jetshower = e[i];
      double pxjetshower = px[i];
      double pyjetshower = py[i];
      double pzjetshower = pz[i];
      double timeplus = timePlus[i];
      p0jetshower = sqrt(pow(pxjetshower, 2) + pow(pyjetshower, 2) + pow(pzjetshower, 2));

      nj0 += 1;
      LBT->KATT1[nj0] = fjetshower;
      LBT->P[1][nj0] = pxjetshower;
      LBT->P[2][nj0] = pyjetshower;
      LBT->P[3][nj0] = pzjetshower;
      LBT->P[0][nj0] = p0jetshower;

      //.... Formation time of the jet parton
//      LBT->tiform[nj0] = 2.0 * p0jetshower / (pow(pxjetshower, 2) + pow(pyjetshower, 2));
      LBT->tiform[nj0] = timeplus;
      LBT->V[1][nj0] = XXX;
      LBT->V[2][nj0] = YYY;
      LBT->V[3][nj0] = 0.0;

      //.... V[0][i] is the distance it travels before the next collision
      for (int i = 0; i <= 3; ++i) {
        LBT->Prad[i][nj0] = LBT->P[i][nj0];
      }
      LBT->tirad[nj0] = LBT->tiform[nj0];
      LBT->tiscatter[nj0] = LBT->tiform[nj0];

      // test
//      cout << "shower " << nj0 << endl;
//      cout << LBT->KATT1[nj0] << endl;
//      for (int j = 0; j < 4; j++) {
//        cout << LBT->P[j][nj0] << endl;
//      }
//      for (int j = 1; j < 4; j++) {
//        cout << LBT->V[j][nj0] << endl;
//      }
      // test
    }

    //.... test		
    LBT->ntest22=0;
    LBT->ntestrad=0;
    //.... test		

    LBT->nj=nj0;
    LBT->np=LBT->nj;


    //.... time evolution in LBT
    double ti=LBT->time0;	//.... initial time inside or outside class?
    int ntimestep=floor(LBT->timend/LBT->dt);	
    for(int timestep=1;timestep<=ntimestep;++timestep) {
      if(LBT->LBTswitch==1) {		
        ti=ti+LBT->dt;
//        cout << "n, ti" << endl;
//        cout << n << " " << ti << endl;
        LBT->LinearBoltzmannTransport(n,ti);
      }
    } //.... time loop end


    int ip = 0;
    for(int i=1;i<=LBT->np;i++)
    {
      if(LBT->P[0][i] == 0) ip+=1;
    }	
    int nnn=LBT->np-ip;
    for(int i=1;i<=LBT->np;i++)
    {
      if(LBT->P[0][i] == 0) continue;
      positiveLBT<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i]<<" "<<LBT->CAT[i]<<endl;
      //positiveLBT<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i]<<" "<<LBT->Vfrozen[1][i]<<" "<<LBT->Vfrozen[2][i]<<" "<<LBT->Vfrozen[3][i]<<" "<<LBT->Vfrozen[0][i]<<" "<<LBT->CAT[i]<<endl;
    }	
    positiveLBT<< "# " <<numEvent<<" "<<nnn<<endl;

    int ip0=0;
    for(int i=LBT->nj+1;i<=LBT->np;i++)
    {
      if(LBT->P0[0][i] == 0) ip0+=1;
    }	
    int nnn0=LBT->np-ip0-LBT->nj;
    for(int i=LBT->nj+1;i<=LBT->np;i++)
    {
      if(LBT->P0[0][i] == 0) continue;
      negativeLBT<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<endl;
      //negativeLBT<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<LBT->Vfrozen0[1][i]<<" "<<LBT->Vfrozen0[2][i]<<" "<<LBT->Vfrozen0[3][i]<<" "<<LBT->Vfrozen0[0][i]<<endl;
    }
    negativeLBT<<"# " <<numEvent<<" "<<nnn0<<endl;	
    //.... time evolution end

    numEvent++;

    //.... Output control
    int print=numEvent%LBT->nprint;
    if(print==0) {
      cout << "n" << "    " << "ntest22" << "    " << "ntestrad" << endl;      
      cout << n << "  " << LBT->ntest22 << " " << LBT->ntestrad << endl;	

      cout << "n" << "    " << "np" << "    " << "LBT->nj" << endl;      
      cout << n << "  " << LBT->np << " " << LBT->nj << endl;			

      //exit(1);
    }

  }
  //.... Event loop end


  cs << sigma/nEvent << "  " << sigmaErr/nEvent << endl;


  cs.close();
  ppdata.close();
  positiveLBT.close();	
  negativeLBT.close();	


  //.... Time counts
  struct tm *local_end;
  time_t time_end;
  time_end = time(NULL);
  local_end = localtime(&time_end);

  char buf2[80];
  strftime(buf2, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_end);
  cout << "the program ends at:" << endl;
  cout << buf2 << endl;

  int cost, nh, nm, ns;
  cost = difftime(time_end, time_start);

  nh = cost / 3600;
  nm = (cost % 3600) / 60;
  ns = (cost % 3600) % 60;

  cout << "the program costs:" << endl;
  cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;
  //.... The program ends
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


