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
//  int jobid = atoi(argv[1]);


  //.... Output files in AA
  sprintf(charName, "positive.txt");
  ofstream positiveLBT(charName);
  sprintf(charName, "negative.txt");
  ofstream negativeLBT(charName);


  //.... PYTHIA 8 Generator. Initialization. 


  //.... LBTclass. Initialization
  LBTclass *LBT = new LBTclass();
//  LBT->jobid = atoi(argv[1]);
  LBT->LBTinitialize();

  //.... LBT input files
  int numEvent=0;


  //.... Event loop
  for(int n=1; ; ++n) {		

    if(numEvent == LBT->ncall) break;


    //.... Clear particle list in LBT
    LBT->LBTclear();  


    //.... Jet shower partons
    int nj0 = 0;
      int fjetshower = 1;
      double p0jetshower = 30;
      double pxjetshower = 0;
      double pyjetshower = 30;
      double pzjetshower = 0;
      double timeplus = 0;
      p0jetshower = sqrt(pow(pxjetshower, 2) + pow(pyjetshower, 2) + pow(pzjetshower, 2));

      nj0 += 1;
      LBT->KATT1[nj0] = fjetshower;
      LBT->P[1][nj0] = pxjetshower;
      LBT->P[2][nj0] = pyjetshower;
      LBT->P[3][nj0] = pzjetshower;
      LBT->P[0][nj0] = p0jetshower;

      //.... Formation time of the jet parton
//      LBT->tiform[nj0] = 2.0 * p0jetshower / (pow(pxjetshower, 2) + pow(pyjetshower, 2));
      LBT->tiform[nj0] = 0;
      LBT->V[1][nj0] = 0;
      LBT->V[2][nj0] = 0;
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
        //cout << "n, ti" << endl;
        //cout << n << " " << ti << endl;
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
      positiveLBT<<LBT->KATT1[i]<<" "<<LBT->P[1][i]<<" "<<LBT->P[2][i]<<" "<<LBT->P[3][i]<<" "<<LBT->P[0][i]<<" "<<LBT->Vfrozen[1][i]<<" "<<LBT->Vfrozen[2][i]<<" "<<LBT->Vfrozen[3][i]<<" "<<LBT->Vfrozen[0][i]<<" "<<LBT->CAT[i]<<endl;
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
      negativeLBT<<LBT->KATT10[i]<<" "<<-LBT->P0[1][i]<<" "<<-LBT->P0[2][i]<<" "<< -LBT->P0[3][i]<<" "<<-LBT->P0[0][i]<<" "<<LBT->Vfrozen0[1][i]<<" "<<LBT->Vfrozen0[2][i]<<" "<<LBT->Vfrozen0[3][i]<<" "<<LBT->Vfrozen0[0][i]<<endl;
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


