#include"LBT.h"



//const double    CA=3.0; 
//const double    sctr=0.197;	     // 1 GeV^-1 = 0.197 fm
//const double    pi=3.1415926; 




//const double  pi=3.1415926;  //conflict with fastjet  ???
//const double    CA=3.0; 
//const double    sctr=0.197;	     // 1 GeV^-1 = 0.197 fm
//const double    pi=3.1415926; 

//...input with current machine time
//...random number seed (any negative integer)
	  
//   long  NUM1 = -33;









//..............................................................subroutine
//..............................................................







void LBTclass::LinearBoltzmannTransport(int &n, double &ti){
		
  int CT;                        //collision type 1 2 3 13 4 5 6 7 8
  int KATTC0;                    //flavor code 1/d 2/u 3/s -1/dbar -2/ubar -3/sbar 21/gluon 
  int KATT2;
  int KATT3;
  int KATT30;
		
  double RTE;                    //scattering rate (energy temperature)
  double E;                      //parton energy
  double T;                      //local temperature

  double T1;                     //index for difference method
  double T2;
  double E1;
  double E2;
  int iT1;
  int iT2;
  int iE1;
  int iE2;	
		
  int nrad;
  int idlead;
  int idlead1;
  int idlead2;
  double Xtau;                   //....................main & titau	
  double Vx;
  double Vy;		
  double Veta;

  double tcar;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar;		
  double ycar;
  double zcar;					

  double tcar0;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar0;		
  double ycar0;
  double zcar0;

  
  double rans;
		
  int KATTx;
  int KATTy;

  double Ejp;
  double Elab;
		
  double qt;
		
  double px0;
  double py0;		
  double pz0;

  double Ncoll22;                    //...................average number of elastic scattering for particle i in the current step  
  int Nscatter;                   //...................number of elastic scattering for particle i in the current step

  int np0=np;                     //...................number of particles in the current step (np0)
                                  //...................number of particles in the beginning of current step (np)

  int Krad;								  
  int free=0;
  int free0=0;
  double ed=0.0;
  double fraction=0.0;
  double VX=0.0;
  double VY=0.0;
  double VZ=0.0;
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////...for heavy quark  
  int nnpp=np;                    //...................number of particles in the current np loop  
  
  double vc0b[4]={0.0};             //flow velocity     
  double pMag,vMag,flowFactor;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////...for heavy quark  
  
  
  double kt2;
  
  double vp0[4]={0.0};
  double p0temp[4]={0.0}; 

  double p0temp1=0.0; 
  double p0temp2=0.0; 
		
  double pcx[4]={0.0};
  double pcy[4]={0.0};		

  double pcx1[4]={0.0};
  double pcy1[4]={0.0};
  
  double Reactionrate;

  double probNcoll22;
  double probNrad;
  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////...for heavy quark  
  double dt_lrf,maxFncHQ,Tdiff,lim_low,lim_high,lim_int;
  double probCol,probRad,probTot;
  int hydro_ctl;

  double KPfactor,KTfactor,Kfactor;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////...for heavy quark

  double dtlastrad; 

  double npstep;  
  
  int nptest; 
  
  //...........................................................................................................NCUT
  ncut=0;
  ncut0=0;
  //...........................................................................................................NCUTEND		  



  Reachtauend=1;




////////////////////////////////////////////////////////////////////////////...2019
//switchtwcoll = 0;
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////...2019
//NRcut = 0;                   //scattering rank
//NScut = 0;                   //splitting rank

//NR0cut = 0;                  //scattering rank for negative
//NS0cut = 0;                  //splitting rank for negative

//dMD = 0.0;                     //distance between mother and daughter parton
//dMDcut = 0.0;                  //distance of decoherence

//collenercut = 0.0;
//collthetacut = 0.0;

//radenercut = 0.0;
//radthetacut = 0.0;

//switchphasecut = 0;
////////////////////////////////////////////////////////////////////////////



















  
  //...collision of the jet parton (i<=nj), recoiled partons, radiated gluons with the background thermal partons
  //...thermal partons
  //...np number of partons in current step
  //...nj number of jet partons 



  
    for(int i=1;i<=np;i++)
	{
	
    icl22=0;
    qt=0;
    idlead=i;
	npstep=np0;
	
//...........test	
//			 cout<<"P[0][i]"<<" "<<P[0][i]<<" "<<P[0][1]<<endl;	
//...........test		
	
	
	
		  
//......propagation of each positive particle and its negative counterpart in the current time step (for all particles)

    if(switchCoLBT_Hydro==0)
	{
		
	      if(P[0][i]==0.0) //anti-overflow NOT NECESSARY
		  {
		  V[1][i]=0.0;
		  V[2][i]=0.0;
		  V[3][i]=0.0;
		  
		  CAT[i]=1;
		  }

	      if(P0[0][i]==0.0) //anti-overflow NOT NECESSARY negative particles
		  {
		  V0[1][i]=0.0;
		  V0[2][i]=0.0;
		  V0[3][i]=0.0;
		  
		  CAT0[i]=1;
		  }
//..........................................................................anti-overflow $$$

		  
//        if(CAT0[i] == 0)
          if(CAT0[i] != 1)	//negative particles  
          {
//...........................................t-z coordinates	
          if(tauswitch==0)
		  {
          V0[1][i]=V0[1][i]+dt*P0[1][i]/P0[0][i];
          V0[2][i]=V0[2][i]+dt*P0[2][i]/P0[0][i];
		  V0[3][i]=V0[3][i]+dt*P0[3][i]/P0[0][i];

		  tcar0=ti;
		  zcar0=V0[3][i];
		  xcar0=V0[1][i];
		  ycar0=V0[2][i];

          double xxx00,yyy00,eta00,tau00,ed00,fraction00,VX00,VY00,VZ00;

		  free0=0;

          xxx00=V0[1][i];
		  yyy00=V0[2][i];
		  eta00=1.0/2.0*log((ti+V0[3][i])/(ti-V0[3][i]));
		  tau00=sqrt(pow(ti,2)-pow(V0[3][i],2));

		  if(switchmedium==0)
		  {
		  if(tau00<tauhydro0 || tau00>tauhydrofile || fabs(xxx00)>10.8 || fabs(yyy00)>10.8 || fabs(eta00)>4.8)     //$$$
		  {
		  free0=1;
		  }
		  }

		  if(free0==0)
          {
		  		  
		  bulklinear(tau00,xxx00,yyy00,eta00,ed00,temp00,fraction00,VX00,VY00,VZ00);
		  
		  vc0[1]=VX00;
		  vc0[2]=VY00;
		  vc0[3]=VZ00;
		  qhat00=DebyeMass2(Kqhat0,alphas,temp00);
		  
          if(fraction00!=0.0)
		  {
          Vfrozen0[0][i]=ti;
          Vfrozen0[1][i]=V0[1][i];
          Vfrozen0[2][i]=V0[2][i];
          Vfrozen0[3][i]=V0[3][i];
          }
		  
		  if(fraction00==0.0)
		  {
		  free0=1;
		  }

          }//if(free0==0)
		  
		  }//if(tauswitch==0) //negative propagation

//...........................................tau-eta coordinates		  
          if(tauswitch==1) //negative propagation
		  {
		  
		  vp0[0]=0.0;
		  vp0[1]=V0[1][i];
		  vp0[2]=V0[2][i];
		  vp0[3]=V0[3][i];

		  pc0[1]=P0[1][i];
		  pc0[2]=P0[2][i];
		  pc0[3]=P0[3][i];
		  pc0[0]=P0[0][i];
		  
		  titau(ti,vc0,vp0,pc0,Vx,Vy,Veta,Xtau);

		  V0[1][i]=V0[1][i]+dt*Vx;
		  V0[2][i]=V0[2][i]+dt*Vy;
		  V0[3][i]=V0[3][i]+dt*Veta;

		  //time and position in Cartesian coordinate			 
		  tcar0=ti*cosh(V0[3][i]);
		  zcar0=ti*sinh(V0[3][i]);
		  xcar0=V0[1][i];
		  ycar0=V0[2][i];

          //...........................................Hydro part		
          double xxx00,yyy00,eta00,tau00,ed00,fraction00,VX00,VY00,VZ00;
		  
		  free0=0;

          xxx00=V0[1][i];
		  yyy00=V0[2][i];
		  eta00=V0[3][i];
		  tau00=ti;

		  //???$$$
		  if(switchmedium==0)
		  {		  
		  if(tau00<tauhydro0 || tau00>tauhydrofile || fabs(xxx00)>10.8 || fabs(yyy00)>10.8 || fabs(eta00)>4.8)
		  {
		  free0=1;
		  }
		  }

          if(free0==0)
          {

		  bulklinear(tau00,xxx00,yyy00,eta00,ed00,temp00,fraction00,VX00,VY00,VZ00);	  
		  
		  vc0[1]=VX00;
		  vc0[2]=VY00;
		  vc0[3]=VZ00;		  
		  qhat00=DebyeMass2(Kqhat0,alphas,temp00);		  
		  
          if(fraction00!=0.0)
		  {
          Vfrozen0[0][i]=tcar0;
          Vfrozen0[1][i]=xcar0;
          Vfrozen0[2][i]=ycar0;
          Vfrozen0[3][i]=zcar0;
          }
		  
		  if(fraction00==0.0)
		  {
		  free0=1;
		  }
		  
          }//if(free==0)		  
		  
		  }//if(tauswitch==1) //negative propagation

		  
          }//if(CAT0[i] != 1) //negative particles 


    }//if(switchCoLBT_Hydro==0)


	      if(CAT[i]!=1)      //positive particles
		  {
		  
		  if(switchCoLBT_Hydro==0)    //$$$ ???
		  {
		  		  
//...........................................t-z coordinates	
          if(tauswitch==0)
		  {		  
          V[1][i]=V[1][i]+dt*P[1][i]/P[0][i];
          V[2][i]=V[2][i]+dt*P[2][i]/P[0][i];
          V[3][i]=V[3][i]+dt*P[3][i]/P[0][i];
		  
		  tcar=ti;
		  zcar=V[3][i];
		  xcar=V[1][i];
		  ycar=V[2][i];

          //...........................................Hydro part		
          double xxx,yyy,eta,tau;

		  free=0;

          xxx=V[1][i];
		  yyy=V[2][i];
		  eta=1.0/2.0*log((ti+V[3][i])/(ti-V[3][i]));
		  tau=sqrt(pow(ti,2)-pow(V[3][i],2));


////////////////////////////////////////////////////////////////////////////...2019
          int IDmom = Mother[i];
		  
          dMD = sqrt(pow((V[1][i]-V[1][IDmom]),2)+pow((V[2][i]-V[2][IDmom]),2)+pow((V[3][i]-V[3][IDmom]),2));
		  
		  //cout<<"ti"<<" "<<ti<<" "<<"i"<<" "<<i<<" "<<"imom"<<" "<<IDmom<<" "<<"NR[i]"<<" "<<NR[i]<<endl;
		  
		  
		  
		  //if(i==1) cout<<"2019--------------------"<<" "<<ti<<endl;
		  

		  //if(dMD<dMDcut)
		  //{
		  //free=1;
		  //}
		  
		  //if(i>nj && ti<tiform[i])
		  //{
		  //free=1;
		  //}	

		  if(i>nj && NR[i]>=NRcut)
		  {
		  free=1;
		  }	

		  //if(NS[i]>=NScut) //another position
		  //{
		  //free=1;
		  //}	
		  
////////////////////////////////////////////////////////////////////////////...2019





		  
		  if(i<=nj && ti<tiform[i])
		  {
		  free=1;
		  }

		  if(Kfishbone==1 && i>nj)
		  {
		  free=1;
		  }		  
	  
		  //???$$$
		  if(switchmedium==0)
		  {
			  
		  if(tau<tauhydro0 || tau>tauhydrofile || fabs(xxx)>10.8 || fabs(yyy)>10.8 || fabs(eta)>4.8)
		  {
		  free=1;	  		  
		  }
		  
//...............................................................805	  
		  if(tau<=tauhydrofile)
		  {
		  Reachtauend=0;	  		  
		  }		  
//...............................................................805
		  
          }		  
		  
          if(free==0)
          {

		  bulklinear(tau,xxx,yyy,eta,ed,temp0,fraction,VX,VY,VZ);
		  
		  vc0[1]=VX;
		  vc0[2]=VY;
		  vc0[3]=VZ;
		  qhat0=DebyeMass2(Kqhat0,alphas,temp0);
		  	  
          if(fraction!=0.0)
		  {
          Vfrozen[0][i]=ti;
          Vfrozen[1][i]=V[1][i];
          Vfrozen[2][i]=V[2][i];
          Vfrozen[3][i]=V[3][i];
          }
		  
		  if(fraction==0.0)
		  {
		  free=1;
		  }
  
		  pc0[1]=P[1][i];
          pc0[2]=P[2][i];
          pc0[3]=P[3][i];
          pc0[0]=P[0][i];
		  
//...............................................................................................NEW ??????????????????????????????????????????????         
          vc0b[1]=VX;
          vc0b[2]=VY;
          vc0b[3]=VZ;		  
		  
          vMag=sqrt(vc0b[1]*vc0b[1]+vc0b[2]*vc0b[2]+vc0b[3]*vc0b[3]);
          flowFactor=(1.0-(pc0[1]*vc0b[1]+pc0[2]*vc0b[2]+pc0[3]*vc0b[3])/pc0[0])/sqrt(1.0-vMag*vMag);			  
//...............................................................................................NEW ??????????????????????????????????????????????
		  
		  //@@@???	  
		  trans(vc0,pc0);
			  
		  if(pc0[0]<sqrt(qhat0))
		  {
		  free=1;
		  }
		  
		  if(i>nj && pc0[0]<Ecut)
		  {
		  free=1;
          }
		  
		  transback(vc0,pc0);

          }//if(free==0)
			 		  
		  }//if(tauswitch==0)	  
		  
		  
//...........................................tau-eta coordinates		
	      if(tauswitch==1)
		  {

		  vp[0]=0.0;
		  vp[1]=V[1][i];
		  vp[2]=V[2][i];
		  vp[3]=V[3][i];

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];

		  titau(ti,vc0,vp,pc0,Vx,Vy,Veta,Xtau);

		  V[1][i]=V[1][i]+dt*Vx;
		  V[2][i]=V[2][i]+dt*Vy;
		  V[3][i]=V[3][i]+dt*Veta;

		  //time and position in Cartesian coordinate			 
		  tcar=ti*cosh(V[3][i]);
		  zcar=ti*sinh(V[3][i]);
		  xcar=V[1][i];
		  ycar=V[2][i];

          //...........................................Hydro part		
          double xxx,yyy,eta,tau;

		  free=0;

          xxx=V[1][i];
		  yyy=V[2][i];
		  eta=V[3][i];
		  tau=ti;

		  if(i<=nj && tcar<tiform[i])
		  {
		  free=1;
		  }

		  //???$$$
		  if(switchmedium==0)
		  {				  
		  if(tau<tauhydro0 || tau>tauhydrofile || fabs(xxx)>10.8 || fabs(yyy)>10.8 || fabs(eta)>4.8)
		  {
		  free=1;
		  }
          } 
		  		  
          if(free==0)
          {

		  bulklinear(tau,xxx,yyy,eta,ed,temp0,fraction,VX,VY,VZ);	  
		  
		  vc0[1]=VX;
		  vc0[2]=VY;
		  vc0[3]=VZ;
		  qhat0=DebyeMass2(Kqhat0,alphas,temp0);

          vp[0]=0.0;
		  vp[1]=V[1][i];
		  vp[2]=V[2][i];
		  vp[3]=V[3][i];

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];  		  

		  titau(ti,vc0,vp,pc0,Vx,Vy,Veta,Xtau);

          if(fraction!=0.0)
		  {
          Vfrozen[0][i]=tcar;
          Vfrozen[1][i]=xcar;
          Vfrozen[2][i]=ycar;
          Vfrozen[3][i]=zcar;
          }

		  if(fraction==0.0)
		  {
		  free=1;
		  }

		  pc0[1]=P[1][i];
          pc0[2]=P[2][i];
          pc0[3]=P[3][i];
          pc0[0]=P[0][i];

//...............................................................................................NEW ??????????????????????????????????????????????????  
          vc0b[1]=VX;
          vc0b[2]=VY;
          vc0b[3]=VZ;

          vMag=sqrt(vc0b[1]*vc0b[1]+vc0b[2]*vc0b[2]+vc0b[3]*vc0b[3]);
          flowFactor=(1.0-(pc0[1]*vc0b[1]+pc0[2]*vc0b[2]+pc0[3]*vc0b[3])/pc0[0])/sqrt(1.0-vMag*vMag);		  
//...............................................................................................NEW ??????????????????????????????????????????????????
		  
		  //@@@???			  
		  trans(vc0,pc0);
			  
		  if(pc0[0]<sqrt(qhat0))
		  {
		  free=1;
		  }
		  
		  if(i>nj && pc0[0]<Ecut)
		  {
		  free=1;
          }
		  
		  transback(vc0,pc0);

          }//if(free==0)

		  }//if(tauswitch==1)

          }//if(switchCoLBT_Hydro==0)



//......propagation of each positive particle and its negative counterpart in the current time step	end	

//......CAT[i]=1 free streaming particle

//......free streaming when energy  < Ecut
//......free streaming when energy0 < Ecut0......in the next step


//................................CoLBT-hydro

          if(switchCoLBT_Hydro==1)     //$$$ ???
          {
		  free=0;
          temp0=tempHydro[i];
          vc0[1]=vxHydro[i];
          vc0[2]=vyHydro[i];
          vc0[3]=vzHydro[i];
          fraction=fractionHydro[i];
		  
		  VX=vc0[1];
		  VY=vc0[2];
		  VZ=vc0[3];		  

		  pc0[1]=P[1][i];
          pc0[2]=P[2][i];
          pc0[3]=P[3][i];
          pc0[0]=P[0][i];
		  
		  trans(vc0,pc0);
			  
		  if(pc0[0]<sqrt(qhat0))
		  {
		  free=1;
		  }
		  
		  if(i>nj && pc0[0]<Ecut)
		  {
		  free=1;
          }
		  
		  transback(vc0,pc0);		  
		  		  
		  if(temp0==0.0)
		  {
		  free=1;
		  }
		  free0=1;		  

          vc0b[1]=VX;
          vc0b[2]=VY;
          vc0b[3]=VZ;
		  
          vMag=sqrt(vc0b[1]*vc0b[1]+vc0b[2]*vc0b[2]+vc0b[3]*vc0b[3]);
          flowFactor=(1.0-(pc0[1]*vc0b[1]+pc0[2]*vc0b[2]+pc0[3]*vc0b[3])/pc0[0])/sqrt(1.0-vMag*vMag);		  
		  
          }

//................................CoLBT-hydro end


//............................................................................................................................................................end of propagation and medium information



//............................................................................................................................................................interaction rate (scattering rate and radiation rate)

          if(free==0)
		  {

          dt_lrf=dt*flowFactor;				  
		  
		  E=P[0][i];
		  T=temp0;
		  KATTC0=KATT1[i];	  		  
		  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new		  
          double qhatTP,RTE1,RTE2;
          double PLen;
		  
          // modification: interplotate rate in l.r.f. 
       	  pc0[1]=P[1][i];
       	  pc0[2]=P[2][i];
       	  pc0[3]=P[3][i];
       	  pc0[0]=P[0][i];		  		  
		  
          trans(vc0,pc0);
          E=pc0[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
          PLen=sqrt(pc0[1]*pc0[1]+pc0[2]*pc0[2]+pc0[3]*pc0[3]);
          transback(vc0,pc0);
	
		  T=temp0;
		  KATTC0=KATT1[i];

//........get scattering number in the current step
		  
		  lam(KATTC0,RTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2); //modified: use P instead	  
		 
          if(tauswitch==0)
		  {
		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];
  
          Ncoll22=fraction*dt_lrf*alphas/0.3*RTE/0.197;		  		  	  		  
          }
//...........................................tau-eta coordinates
          if(tauswitch==1)
		  {
		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];

          Ncoll22=fraction*dt_lrf*alphas/0.3*RTE*Xtau/0.1970;			  	  
		  }


//.......get scattering number in the current step end
		  
		  probNcoll22=0.0;
		  probNrad=0.0;

//.......get radiation number in the current step		  
		  		  
          if(Kradiation==1)
		  {
// get qhat from table????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
                  if(KATTC0==21) {
                      RTE1=(qhatG[iT2][iE1]-qhatG[iT1][iE1])*(T-T1)/(T2-T1)+qhatG[iT1][iE1];
                      RTE2=(qhatG[iT2][iE2]-qhatG[iT1][iE2])*(T-T1)/(T2-T1)+qhatG[iT1][iE2];
                  } else if(KATTC0==4||KATTC0==-4) {
                      RTE1=(qhatHQ[iT2][iE1]-qhatHQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE1];
                      RTE2=(qhatHQ[iT2][iE2]-qhatHQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE2];
                  } else {
                      RTE1=(qhatLQ[iT2][iE1]-qhatLQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE1];
                      RTE2=(qhatLQ[iT2][iE2]-qhatLQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE2];
                  }

                  qhatTP=(RTE2-RTE1)*(PLen-E1)/(E2-E1)+RTE1;

                  //qhatTP=qhatTP*Kfactor;
				  //qhatTP=qhatTP;
				  
				  //...................................alphas/0.3
				  qhatTP=pow(alphas/0.3,2)*qhatTP;				  
				  
                  qhat_over_T3=qhatTP;  // what is read in is qhat/T^3 of quark
                  if(KATTC0==21) 
				  {
				  D2piT=8.0*pi/(qhat_over_T3/CA*CF);
                  } else {
				  D2piT=8.0*pi/qhat_over_T3;
				  }				  
				  
                  qhat=qhatTP*pow(T,3); // for light quark and gluon, need additional table to make it correct

// get qhat from table???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????? end

//........get qhat from a analytical formula 		  
          if(abs(KATT1[i])==1 || abs(KATT1[i])==2 || abs(KATT1[i])==3)
          {
          isp=1;
		  //qhat=4.0/3.0*42*1.2/3.1415926*pow(alphas,2)*pow(T,3)*log(5.8*E*T/(4.0*qhat0));
          n_sp1 +=1;
          }
          else if(KATT1[i]==21)
          {
          isp=2;
		  //qhat=3.0*42*1.2/3.1415926*pow(alphas,2)*pow(T,3)*log(5.8*E*T/(4.0*qhat0));
          n_sp2 +=1;
          }

//........time interval since last radiation (ti-tirad)

//........time step in the local rest frame ???

          Tint_lrf[i]=Tint_lrf[i]+dt_lrf;
		  
          //Tdiff=Tint_lrf[i];

		  Tdiff=abs(tcar-tirad[i])*flowFactor;
		  dtlastrad=abs(tcar-tirad[i])*flowFactor;
	  

		  Elab=pc0[0]; // p4-the initial 4-momentum of the jet parton in the lab frame
          trans(vc0,pc0);
          Ejp=pc0[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
          transback(vc0,pc0);

		  
          if(tauswitch==0)
		  {	  
          if(KATT1[i]==21)
          {			  
          radng[i]=alphas/0.3*1.0/2.0*nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ);
          }
          if(KATT1[i]!=21)
          {			  
          radng[i]=alphas/0.3*nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ);
          }      		  
          }
		  if(tauswitch==1)
		  {
          radng[i]=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*Xtau;
          }
		  
          lim_low=sqrt(6.0*pi*alphas)*temp0/E;
          if(abs(KATT1[i]==4)) lim_high=1.0;
          else lim_high=1.0-lim_low;
          lim_int=lim_high-lim_low;
          if(lim_int>0.0) probRad=1.0-exp(-radng[i]);
          else probRad=0.0;
		  
          if(radng[i]<0.00000000001) radng[i]=0.0;

		  
		  }//if(Kradiation==1)


//.......get radiation number in the current step end		  	
  
          probNcoll22=1.0-exp(-Ncoll22);
		  
          probNrad=1.0-exp(-radng[i]);

		  //probNrad=probRad;		  
		  
          if(Kradiation==1)
          {
          Reactionrate=probNcoll22*(1.0-probNrad)+probNrad;
          }
          else
          {
		  probNrad=0.0;
          Reactionrate=probNcoll22;
          }

//............................................................................................................................................................interaction rate (scattering rate and radiation rate) end




//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////...single interaction




          int switchcoll = 1;
          int switchrad = 1;

           if(Singlestepswitch == 1)
		   {
           if(Force22 == 1)
		   {

//...Force interaction
           switchcoll = 2;           
           switchrad = 2;
//...Force interaction

           }

           if(Force23 == 1)
		   {

//...Force interaction
           switchcoll = 2;           
           switchrad = 2;
//...Force interaction

           }
		   
           if(Force2n == 1)
		   {

//...Force interaction
           switchcoll = 2;           
           switchrad = 3;
//...Force interaction

           }
           }		   
		  
//........interaction begin 
		  
          if(ran0(&NUM1)<Reactionrate || switchcoll > 1) // !Yes, collision!		  		  
		  //if(V[0][i]<0.0) // !Yes, collision!
		  {
			  Nscatter=KPoisson(Ncoll22);
              
			  if(Nscatter<=1)
			  {
              Nscatter=1;
			  }
			  
			  Nscatter=1;
			  
			  
              for(int nsca=1;nsca<=Nscatter;nsca++)
			  {			
			
			  np0=np0+1;
			  pc0[1]=P[1][i];
			  pc0[2]=P[2][i];
			  pc0[3]=P[3][i];
			  pc0[0]=P[0][i];

			  //...........collision between a jet parton or a recoiled parton and a thermal parton from the medium 
			  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new			  
			  
			  flavor(CT,KATTC0,KATT2,KATT3,RTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2);
              if(CT==11||CT==12) { // for heavy quark scattering
			  collHQ22(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
              } else { // for light parton scattering
			  colljet22(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
              }			  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new

              icl22 = 1; //???
			  
			  ntest22=ntest22+1;
			 
			  KATT1[i]=KATTC0;
			  KATT1[np0]=KATT2;
			  KATT10[np0]=KATT3;

//???			  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new			  
//			  if(pc0[0]<pc2[0] && abs(KATTC0)!=4 && KATTC0==KATT2) { //disable switch for heavy quark, only allow switch for identical particles			  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new


//.........................................................................................................................66

                 if(ExchangeIDswitch == 1)
				 {	 
			     if(pc0[0]<pc2[0])
				 {
				  //if(Singlestepswitch==0 && Force22!=1 && Force23==0)
				  {
				  for(int k=0;k<=3;k++)
					{                 
					  p0temp[k]=pc2[k];
					  pc2[k]=pc0[k];
					  pc0[k]=p0temp[k];
					}
			      KATT1[i]=KATT2;
     			  KATT1[np0]=KATTC0;
				  }
				 }
				 }//if(Singlestepswitch == 0)
//.........................................................................................................................66
//
			    for(int j=0;j<=3;j++)
				{
//			    if(i<=nj) // always make sure the jet parton carries the largest momentum (maybe now all the parton includes recoiled for there could be multiple scattering in a single time step)
				
				P[j][i]=pc0[j];
				P[j][np0]=pc2[j];
				V[j][np0]=V[j][i];
				
//		        if(i<=nj || nsca!=1) // then for each recoiled parton there is a negative parton
//				{
				P0[j][np0]=pc3[j];
				V0[j][np0]=V[j][i];
//				}
				}

////////////////////////////////////////////////////////////////////////////...2019
                Mother[np0] = i;
				tiform[np0] = tcar;
				
				NR[i] = NR[i]+1;
////////////////////////////////////////////////////////////////////////////...2019


                CAT[np0]=2;

//.............................................................. collision between "negative" partons

//...background thermal parton scattering with the same thermal parton the associated parton has scattered.
//...negative parton can only get the rank i>nj and it can only act on the first collision in the current time step. 
//...There is no physical reason. only because the recoiled parton and its negative partner appear in pairs. 

//...Scattering between P0[i] and the newest negative parton  

                if(CAT0[i] == 0)
                {

                if(switchtwcoll==0) free0=1;
                
                if(free0 == 0)
                {
		        if(i>nj && nsca==1)
				{

                temp00=temp0;

		        qhat00=DebyeMass2(Kqhat0,alphas,temp00);

                pc00[1]=P0[1][i];
                pc00[2]=P0[2][i];
                pc00[3]=P0[3][i];
                pc00[0]=P0[0][i];
				
				pc30[1]=pc3[1];
				pc30[2]=pc3[2];
				pc30[3]=pc3[3];
				pc30[0]=pc3[0];

	            E=pc00[0];
	            T=temp00;
	            KATTC0=KATT10[i];
				KATT30=KATT3;
                
				twflavor(CT,KATTC0,KATT30,E,T);		 			
					   				   
                twcoll(CT,qhat00,vc0,pc00,pc30);		   
				   				   				   
	            KATT1[np0]=KATT2;
	            KATT10[i]=KATTC0;
	            KATT10[np0]=KATT30;

                for(int j=0;j<=3;j++)
				{
                 P0[j][i]=pc00[j];
                 P0[j][np0]=pc30[j];
                 V0[j][np0]=V[j][i];
                }
               	V0[0][np0]=-log(1.0-ran0(&NUM1));

////////////////////////////////////////////////////////////////////////////...2019
                //Mother[np0] = i;
				//tiform[np0] = tcar;
				
				NR0[i] = NR0[i]+1;
////////////////////////////////////////////////////////////////////////////...2019	


				}//if(i>nj && nsca==1)
                }//if(free0 == 0)
                }//if(CAT0[i] == 0)

//.............................................................. "negative" end
											
             
			  }//for(int nsca=1;nsca<=Nscatter;nsca++)			  			  



//............test

//              cout<<"colljet22 P[0][i]"<<" "<<P[0][i]<<" "<<"colljet22 P[0][np0]"<<" "<<P[0][np0]<<endl;

//............test

			  
//...............................................................................................the point where the elastic scattering end	




//work need to be done		  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			  			  
			    if(ran0(&NUM1)<probNrad/Reactionrate && Kradiation==1 && qt!=0) Krad=1;
                else Krad=0;			  			  	  			  
			  
			    if(Kradiation==1 && Singlestepswitch==1)
				{
                if(Force23==1 || Force2n==1)
                {
				Krad=1;						
				}
				}


			  
			    if(Krad==1) // Switch on the radiation processes
				{			  					  

						    for(int j=0;j<=3;j++) // prepare for the colljet23
							{
							  pc01[j]=pc4[j];							  
							  pb[j]=pc4[j];
							}
							
							
						  //colljet23(temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,tcar,tiscatter[i],tirad[i],dtlastrad,Elab,Ejp);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new							  
		                  collHQ23(KATT1[i],temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,Tdiff,Ejp,maxFncHQ,lim_low,lim_int);						  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new							  

			              ntestrad=ntestrad+1;					  
						  
						  if(icl23!=1)
						  {
							
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new	
                          int ctGluon=1;							
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new								
							
                              ng_coll23+=1; 
							  nrad=KPoisson(radng[i]);   
							  
                              ng_nrad+=nrad;
				  
							  np0=np0+1;
                              int npg=np0;  //....radiated gluons' id starts

							  KATT1[np0]=21;				  
				  
							  for(int j=0;j<=3;j++)
								{
                                  //re-assignment
								  P[j][i]=pc01[j];            //jet parton after colljet23
								  P[j][np0-1]=pc2[j];
								  P0[j][np0-1]=pc3[j];

								  P[j][np0]=pc4[j];           //radiated gluon from colljet23
                   				  V[j][np0]=V[j][i];
								  P0[j][np0]=0.0;
								  V0[j][np0]=0.0;
								  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new		  
                                  Vfrozen[j][np0]=V[j][i];
                                  Vfrozen0[j][np0]=0.0;
                                  Tfrozen[np0]=temp0;
                                  Tfrozen0[np0]=0.0;								  
                                  if(j!=0)
								  {
                                  vcfrozen[j][np0]=vc0b[j];
                                  vcfrozen0[j][np0]=0.0;
                                  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new								  
								}

								
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new									
                          Tfrozen[np0]=temp0;
                          Tfrozen0[np0]=0.0;

                          // pass initial pT and weight information
                          P[5][np0]=P[5][i];
                          P0[5][np0]=P0[5][i];
                          WT[np0]=WT[i];
                          WT0[np0]=WT0[i];

                          Tint_lrf[i]=0.0;  //reset radiation infomation for heavy quark

                          eGluon=eGluon+pc4[0];
                          nGluon=nGluon+1.0;

               	          V[0][np0]=-log(1.0-ran0(&NUM1));								
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new									

			
//                             for(int j=0; j<=3; ++j)
//                             {
//              				     Prad[j][i]=P[j][i];
//              				     Prad[j][np0-1]=P[j][np0-1];
//              				     Prad[j][np0]=P[j][np0];
//                             }							 													 
							 

               	              V[0][np0]=-log(1.0-ran0(&NUM1));

			                  // rotate to the frame in which jet moves along z-axis
			                  rotate(pb[1],pb[2],pb[3],pc4,1);

			                  // calculate qhat in lab frame since last radiation   
			                  kt2=pow(pc4[1],2)+pow(pc4[2],2);

			                  rotate(pb[1],pb[2],pb[3],pc4,-1);
								
							  tiform[np0]=tcar+Elab*(pc4[0]/Elab)*(1-pc4[0]/Elab)/pow(kt2,2);
							  

////////////////////////////////////////////////////////////////////////////...2019
                Mother[np0] = i;
				//tiform[np0] = tcar;
				
				NS[i] = NS[i]+1;
////////////////////////////////////////////////////////////////////////////...2019


							  //...colljet23 done				  				  



//............test

              //cout<<"colljet23 P[0][i]"<<" "<<P[0][i]<<" "<<"colljet23 P[0][np0]"<<" "<<P[0][np0]<<" "<<"colljet23 P[0][np0-1]"<<" "<<P[0][np0-1]<<endl;

			  //cout<<np0<<endl;
			  
			  nptest=np0;
			  
//............test

                if(Singlestepswitch==1)
				{
                if(Force2n==-1)
                {
				nrad=1;
                }
                }				

							  
							  if(nrad>1)
								{
//                                    idlead1=i;
//                                    idlead2=npstep+1;																		
									
								Doll33:

//............test
              //cout<<"Doll33 start"<<endl;

              //cout<<"pc2"<<" "<<pc2[0]<<" "<<"pc01"<<" "<<pc01[0]<<" "<<"pc4"<<" "<<pc4[0]<<" "<<"pb"<<" "<<pb[0]<<endl;			  
//............test
                                  
//............test								  
								  for(int j=0;j<=3;j++) pc2[j]=pc4[j];
//............test
								
								  //radiation(qhat0,vc0,pc01,pc2,pc4,pb,iclrad,tcar,tiscatter[i],tirad[i],dtlastrad,Elab,Ejp);	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new
			                      radiationHQ(KATT1[i],qhat0,vc0,pc2,pc01,pc4,pb,iclrad,Tdiff,Ejp,maxFncHQ,temp0,lim_low,lim_int);  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...new                   


//............test
              //cout<<"Doll33 end"<<endl;

              //cout<<"pc2"<<" "<<pc2[0]<<" "<<"pc01"<<" "<<pc01[0]<<" "<<"pc4"<<" "<<pc4[0]<<" "<<"pb"<<" "<<pb[0]<<endl;			  
//............test



			                      ntestrad=ntestrad+1;								  
			  
								  if(iclrad!=1)
								  {
							      np0=np0+1;
							      ng_radiation+=1;
                                      
							      KATT1[np0]=21;
								  
//                                cout << "add one gluon" << endl; //???
                                  ctGluon++;
									
								  for(int j=0;j<=3;j++)
								  {
                                  //P[j][idlead1]=pc01[j];
                                  //P[j][idlead2]=pc2[j];

                                  P[j][i]=pc01[j];
                                  P[j][np0-1]=pc2[j];								  
								  
                                  P[j][np0]=pc4[j];
                   				  V[j][np0]=V[j][i];
								  P0[j][np0]=0.0;
								  V0[j][np0]=0.0;
  
                                  Vfrozen[j][np0]=V[j][i];
                                  Vfrozen0[j][np0]=0.0;
                                  if(j!=0) 
							      {                                 
                                  vcfrozen[j][np0]=vc0b[j];
                                  vcfrozen0[j][np0]=0.0;
                                  }
										  
								  }

//                                  for(int j=0; j<=3; ++j)
//                                  {
//                       		      Prad[j][idlead1]=P[j][idlead1];
//                       			  Prad[j][idlead2]=P[j][idlead2];
//                       		      Prad[j][np0]=P[j][np0];
//                                  }
								  

                                  Tfrozen[np0]=temp0;
                                  Tfrozen0[np0]=0.0;
 
 
 
///////////////////////////////////////////////////////////////////////////////////// 
                                  // pass initial pT and weight information
                                  P[5][np0]=P[5][i];
                                  P0[5][np0]=P0[5][i];
                                  WT[np0]=WT[i];
                                  WT0[np0]=WT0[i];

								  
                                  eGluon=eGluon+pc4[0];
                                  nGluon=nGluon+1.0;
/////////////////////////////////////////////////////////////////////////////////////
								  
								  
								  
                      	          V[0][np0]=-log(1.0-ran0(&NUM1));


			                      rotate(pb[1],pb[2],pb[3],pc2,1);
  
			                      kt2=pow(pc2[1],2)+pow(pc2[2],2);

			                      rotate(pb[1],pb[2],pb[3],pc2,-1);
								
							      tiform[np0]=tcar+Elab*(pc2[0]/Elab)*(1-pc2[0]/Elab)/pow(kt2,2);										


////////////////////////////////////////////////////////////////////////////...2019
                Mother[np0] = i;
				//tiform[np0] = tcar;
				
				NS[i] = NS[i]+1;
////////////////////////////////////////////////////////////////////////////...2019									  
								  
								  
								  
								  
								  
//............test

              //cout<<np0<<endl;

              //cout<<"radiation P[0][nptest]"<<" "<<P[0][nptest]<<endl;

              //cout<<"radiation P[0][i]"<<" "<<P[0][i]<<" "<<"radiation P[0][np0]"<<" "<<P[0][np0]<<" "<<"radiation P[0][np0-1]"<<" "<<P[0][np0-1]<<endl;

//............test									  
								  
								  
								  
								  
								  
								  
								  
								  
								  
							      if(nrad>2)
								  {				  
								  nrad=nrad-1;
										  
								  goto Doll33;
								  }
								  
								  //133................................................................................................end
								  
							      }//if(icl!=1)radiation
			   
								}//if(nrad>1)

								
								
//............test

              //cout<<"radiation P[0][i]"<<" "<<P[0][i]<<endl;

//............test								
								
								
								
								
								
								
								
                                //....tirad information 
                                //....the jet parton and the recoiled parton
                                tirad[i]=tcar;
                                tirad[npg-1]=tcar;
                                //....radiated gluons in radiation()
	                            for(unsigned ig=npg; ig<=np0; ++ig)
                                {
                                tirad[ig]=tcar;
                                }

							} //if(icl!=1)colljet23			  
				  				  
				}  //if(Krad==1)


//...........collision end				  
//...........determine the leading partons and set radng[i]
			  
              //....tiscatter information 
			
              //tiscatter[i]=tcar;
              V[0][i]=-log(1.0-ran0(&NUM1));

//........................................................................................................				
             for(unsigned ip=npstep+1; ip<=np0; ++ip)
             {
             //tiscatter[ip]=tcar;
             V[0][ip]=-log(1.0-ran0(&NUM1));
             V0[0][ip]=-log(1.0-ran0(&NUM1));

             Vfrozen[0][ip]=ti;
             Vfrozen[1][ip]=V[1][ip];
             Vfrozen[2][ip]=V[2][ip];
             Vfrozen[3][ip]=V[3][ip];				
				
             Vfrozen0[0][ip]=ti;
             Vfrozen0[1][ip]=V0[1][ip];
             Vfrozen0[2][ip]=V0[2][ip];
             Vfrozen0[3][ip]=V0[3][ip];
             }

            //....find the jet parton with maximal energy
			//....find the ID of the parton [idlead] with maximal energy produced in this step [i], compare with parton i, exchange their identity if P[0][idlead]>P[0][i]
             idlead=npstep+1;
             p0temp1=P[0][idlead];
             for(unsigned ip=npstep+1; ip<=np0-1; ++ip)
             {
                 if(p0temp1 < P[0][ip+1])
                 {
                     p0temp1=P[0][ip+1];
                     idlead=ip+1;
                 }
             }

			 
//...........test			 
			 
			 //cout<<"P[0][i]"<<" "<<P[0][i]<<endl;
			 //cout<<"P[0][idlead]"<<" "<<P[0][idlead]<<endl;
			 
//...........test			 
			 
			  if(P[0][idlead]>P[0][i])   //i<=nj
                {	

	              KATTx=KATT1[i];
				  KATTy=KATT1[idlead];
				  KATT1[i]=KATTy;
		          KATT1[idlead]=KATTx;					

				  for(int j=0;j<=3;j++)
					{
					  pcx[j]=P[j][i];
					  pcy[j]=P[j][idlead];
					  P[j][i]=pcy[j];
					  P[j][idlead]=pcx[j];
					}
                }
								
			  //...........sorting block i npstep~np0(positive) npstep+1(negative)

			  //........................................................................................................
			  
			  //...........branch of the recoiled parton
			  for(int m=npstep+1;m<=np0;m++)
				{
				  if(CAT[i]==2)
					{			
					  CAT[m]=2;             				 
					}
				}

			}
		  //...!Yes, collision!			
        }
        //....for if(free==0)
        }
	  //..........if CAT[i]=0 end


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	  
	}
    //......for np end

  //........time step end, np: the number of particles at this point  		
  np=np0;  
  
  if(Kprimary==1)
  {
  np=nj;
  }

}
		
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




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

float LBTclass::ran0(long *idum)

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
	
//.........................................................................
double LBTclass::alphas0(int &Kalphas,double temp0){
      
  double X;
  if(Kalphas==1) 
	{      
      X=0.3;      
	}
  return X;
}
//.........................................................................
double LBTclass::DebyeMass2(int &Kqhat0,double alphas,double temp0){

  double Y;
  if(Kqhat0==1)
	{      
      Y=4.0*pi*alphas*pow(temp0,2);      
	}
  if(Kqhat0==2)
	{      
      Y=(3.0/2.0)*4.0*pi*alphas*pow(temp0,2);      
	}
  if(Kqhat0==3)
	{      
      Y=1.0;      
	}
  return Y;

}
//.........................................................................

//.........................................................................	
void LBTclass::titau(double ti,double vf[4],double vp[4],double p0[4],double &Vx,double &Vy,double &Veta,double &Xtau){   	

  //..............................................................test part
  //		  cout<<"ti"<<" "<<ti<<" "<<"vf"<<" "<<vf[1]<<" "<<vf[2]<<" "<<vf[3]<<endl;
  //		  cout<<"vp[4]"<<" "<<vp[1]<<" "<<vp[2]<<" "<<vp[3]<<" "<<vp[0]<<endl;	
  //		  cout<<"p0[4]"<<" "<<p0[1]<<" "<<p0[2]<<" "<<p0[3]<<" "<<p0[0]<<endl;		  
  //..............................................................test part

  //....notice the form of vf
  double gamma=1.0/sqrt(1-(vf[1]*vf[1]+vf[2]*vf[2]));
  double mt=sqrt(p0[1]*p0[1]+p0[2]*p0[2]);
  double Yp=1.0/2.0*log((p0[0]+p0[3])/(p0[0]-p0[3]));	
  double etas=vp[3];
  double etaf=atanh(vf[3])+etas;
  double pper=sqrt(p0[1]*p0[1]+p0[2]*p0[2]);
  double vper=sqrt(vf[1]*vf[1]+vf[2]*vf[2]);
  double pvper=p0[1]*vf[1]+p0[2]*vf[2];

  Vx=p0[1]/pper/cosh(Yp-etas);	  
  Vy=p0[2]/pper/cosh(Yp-etas);
  Veta=(1.0/ti)*tanh(Yp-etas);
	  
  Xtau=(gamma*mt*cosh(Yp-etaf)-pvper*vper)/(mt*cosh(Yp-etas));

  //..............................................................test part
  //		  cout<<"gamma"<<" "<<gamma<<" "<<"mt"<<" "<<mt<<" "<<"Yp"<<" "<<Yp<<endl;
  //		  cout<<"etas"<<" "<<etas<<" "<<"etaf"<<" "<<etaf<<" "<<"pper"<<" "<<pper<<endl;	
  //		  cout<<"vper"<<" "<<vper<<" "<<"pvper"<<" "<<pvper<<endl;
  //		  cout<<"Vx"<<" "<<Vx<<" "<<"Vy"<<" "<<Vy<<" "<<"Veta"<<" "<<Veta<<endl;
  //		  cout<<"Xtau"<<" "<<Xtau<<endl;		  
  //..............................................................test part	  
}
//.........................................................................	
	
//.........................................................................


/*	
void LBTclass::lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){   
	  
  double dtemp=0.02;
  iT1=floor((T-0.1)/dtemp);
  iT2=iT1+1;
  iE1=floor(log(E)+2);
  iE2=iE1+1;

  T1=0.12+(iT1-1)*0.02;
  T2=T1+dtemp;
  E1=exp(iE1-2.0);
  E2=exp(iE2-2.0);
	  
  if(KATT0==21)
	{	
	  double RTE1=(Rg[iT2][iE1]-Rg[iT1][iE1])*(T-T1)/(T2-T1)+Rg[iT1][iE1];
	  double RTE2=(Rg[iT2][iE2]-Rg[iT1][iE2])*(T-T1)/(T2-T1)+Rg[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
	}

  if(KATT0!=21)
	{    
	  double RTE1=(Rq[iT2][iE1]-Rq[iT1][iE1])*(T-T1)/(T2-T1)+Rq[iT1][iE1];
	  double RTE2=(Rq[iT2][iE2]-Rq[iT1][iE2])*(T-T1)/(T2-T1)+Rq[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	}
	   
}
*/


void LBTclass::lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){   
	  
  double dtemp=0.02;
  iT1=floor((T-0.1)/dtemp);
  iT2=iT1+1;
  iE1=floor(log(E)+2);
  iE2=iE1+1;

  T1=0.12+(iT1-1)*0.02;
  T2=T1+dtemp;
  E1=exp(iE1-2.0);
  E2=exp(iE2-2.0);
	  
  if(KATT0==21) {	
	  double RTE1=(Rg[iT2][iE1]-Rg[iT1][iE1])*(T-T1)/(T2-T1)+Rg[iT1][iE1];
	  double RTE2=(Rg[iT2][iE2]-Rg[iT1][iE2])*(T-T1)/(T2-T1)+Rg[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
  } else if (KATT0==4||KATT0==-4) { // add heavy quark channel
	  double RTE1=(RHQ[iT2][iE1]-RHQ[iT1][iE1])*(T-T1)/(T2-T1)+RHQ[iT1][iE1];
	  double RTE2=(RHQ[iT2][iE2]-RHQ[iT1][iE2])*(T-T1)/(T2-T1)+RHQ[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
  } else {    
	  double RTE1=(Rq[iT2][iE1]-Rq[iT1][iE1])*(T-T1)/(T2-T1)+Rq[iT1][iE1];
	  double RTE2=(Rq[iT2][iE2]-Rq[iT1][iE2])*(T-T1)/(T2-T1)+Rq[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
//          cout<<"RTE2,RTE1,E,E1,E2,RTE: "<<RTE2<<"  "<<RTE1<<"  "<<E<<"  "<<E1<<"  "<<E2<<"  "<<RTE<<endl;
  }

}



























//.........................................................................

void LBTclass::flavor(int &CT,int &KATT0,int &KATT2,int &KATT3,double RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){  

  double RTEg;
  double RTEg1;
  double RTEg2;
  double RTEg3;
  double RTEq;
  double RTEq3;
  double RTEq4;
  double RTEq5;
  double RTEq6;
  double RTEq7;
  double RTEq8;

  int vb[7]={0};
  int b=0;
  int KATT00=KATT0;


  vb[1]=1; 
  vb[2]=2; 
  vb[3]=3; 
  vb[4]=-1; 
  vb[5]=-2; 
  vb[6]=-3;

	  
  //.....
  linear(KATT0,E,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2,RTEg,RTEg1,RTEg2,RTEg3,RTEq,RTEq3,RTEq4,RTEq5,RTEq6,RTEq7,RTEq8);

  //.....for gluon
  if(KATT00==21)
	{
	  double R0=RTE;
	  double R1=RTEg1;
	  double R2=RTEg2;
	  double R3=RTEg3;
	   	
	  double a=ran0(&NUM1);
	   
	  if(a<=R1/R0)
        { 
	      CT=1;
	      KATT3=21;
	      KATT2=21;
	      KATT0=21;
	    }	  
	  
	  if(a>R1/R0 && a<=(R1+R2)/R0)
	    {
	      CT=2;
	      b=floor(ran0(&NUM1)*6+1);
	      if(b==7)
		    {
			  b=6;
	        }
	      KATT3=21;
	      KATT2=vb[b];
	      KATT0=-KATT2;
	    }	  
	  
	  if(a>(R1+R2)/R0 && a<=1.0)
	    {
	      CT=3;
	      b=floor(ran0(&NUM1)*6+1);
	      if(b==7)
		    { 
			  b=6;
	        }
	      KATT3=vb[b];
	      KATT2=KATT3;
	      KATT0=21;
	    }
	}
	  
  //.....for quark and antiquark
  if(KATT00!=21)
	{      
	  double R00=RTE;
	  double R3=RTEq3;
	  double R4=RTEq4;
	  double R5=RTEq5;
	  double R6=RTEq6;
	  double R7=RTEq7;
	  double R8=RTEq8;

	  double a=ran0(&NUM1);
	  if(a<=R3/R00)
	    { 
	      CT=13;
	      KATT3=21;
          KATT2=21; 
	      KATT0=KATT0;
	    }
	  	  
	  if(a>R3/R00 && a<=(R3+R4)/R00)
		{ 
	      CT=4;
		f1:	      b=floor(ran0(&NUM1)*6+1);
		  if(b==7)
		    {
			  b=6;
	        }
		  KATT3=vb[b]; 
		  if(KATT3==KATT0)
			{
			  goto f1;
	        }
		  KATT2=KATT3;
		  //	      KATT0=KATT0
	    }
	  	  
	  if(a>(R3+R4)/R00 && a<=(R3+R4+R5)/R00)
	    { 
	      CT=5;
	      KATT3=KATT0;
	      KATT2=KATT0;
	    }
	  //.....the only difference between quark and antiquark	   	   
	    
	  if(a>(R3+R4+R5)/R00 && a<=(R3+R4+R5+R6)/R00)
		{
          CT=6;
	      KATT3=-KATT0;  
		f2:	      b=floor(ran0(&NUM1)*3+1);
		  if(b==4)
		    {
			  b=3;
	        }
		  KATT2=-KATT0/abs(KATT0)*vb[b]; 
		  if(abs(KATT2)==abs(KATT3))
			{
			  goto f2;
            }
		  KATT0=-KATT2;
	    }
	   	   
	  if(a>(R3+R4+R5+R6)/R00 && a<=(R3+R4+R5+R6+R7)/R00)
        {	   
	      CT=7;
	      KATT3=-KATT0; 
	      KATT2=KATT3;
		  //	      KATT0=KATT0
	    }	   
	   
	  if(a>(R00-R8)/R00 && a<=1.0)
	    {
		  CT=8;
		  KATT3=-KATT0; 
		  KATT2=21; 
		  KATT0=21;
	    }	   	  
	}	
}
 	
  

//.........................................................................

void LBTclass::linear(int KATT,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2,double &RTEg,double &RTEg1,double &RTEg2,double &RTEg3,double &RTEq,double &RTEq3,double &RTEq4,double &RTEq5,double &RTEq6,double &RTEq7,double &RTEq8){
  if(KATT==21)
	{
	  double RTE1=(Rg1[iT2][iE1]-Rg1[iT1][iE1])*(T-T1)/(T2-T1)+Rg1[iT1][iE1];
	  double RTE2=(Rg1[iT2][iE2]-Rg1[iT1][iE2])*(T-T1)/(T2-T1)+Rg1[iT1][iE2];
	  RTEg1=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	
	  RTE1=(Rg2[iT2][iE1]-Rg2[iT1][iE1])*(T-T1)/(T2-T1)+Rg2[iT1][iE1];
	  RTE2=(Rg2[iT2][iE2]-Rg2[iT1][iE2])*(T-T1)/(T2-T1)+Rg2[iT1][iE2];
	  RTEg2=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	
	
	  RTE1=(Rg3[iT2][iE1]-Rg3[iT1][iE1])*(T-T1)/(T2-T1)+Rg3[iT1][iE1];
	  RTE2=(Rg3[iT2][iE2]-Rg3[iT1][iE2])*(T-T1)/(T2-T1)+Rg3[iT1][iE2];
	  RTEg3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;		
     
	}

  if(KATT!=21)
	{
	  double RTE1=(Rq3[iT2][iE1]-Rq3[iT1][iE1])*(T-T1)/(T2-T1)+Rq3[iT1][iE1];
	  double RTE2=(Rq3[iT2][iE2]-Rq3[iT1][iE2])*(T-T1)/(T2-T1)+Rq3[iT1][iE2];
	  RTEq3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq4[iT2][iE1]-Rq4[iT1][iE1])*(T-T1)/(T2-T1)+Rq4[iT1][iE1];
	  RTE2=(Rq4[iT2][iE2]-Rq4[iT1][iE2])*(T-T1)/(T2-T1)+Rq4[iT1][iE2];
	  RTEq4=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq5[iT2][iE1]-Rq5[iT1][iE1])*(T-T1)/(T2-T1)+Rq5[iT1][iE1];
	  RTE2=(Rq5[iT2][iE2]-Rq5[iT1][iE2])*(T-T1)/(T2-T1)+Rq5[iT1][iE2];
	  RTEq5=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq6[iT2][iE1]-Rq6[iT1][iE1])*(T-T1)/(T2-T1)+Rq6[iT1][iE1];
	  RTE2=(Rq6[iT2][iE2]-Rq6[iT1][iE2])*(T-T1)/(T2-T1)+Rq6[iT1][iE2];
	  RTEq6=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq7[iT2][iE1]-Rq7[iT1][iE1])*(T-T1)/(T2-T1)+Rq7[iT1][iE1];
	  RTE2=(Rq7[iT2][iE2]-Rq7[iT1][iE2])*(T-T1)/(T2-T1)+Rq7[iT1][iE2];
	  RTEq7=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq8[iT2][iE1]-Rq8[iT1][iE1])*(T-T1)/(T2-T1)+Rq8[iT1][iE1];
	  RTE2=(Rq8[iT2][iE2]-Rq8[iT1][iE2])*(T-T1)/(T2-T1)+Rq8[iT1][iE2];
	  RTEq8=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	}	
}

//.........................................................................

void LBTclass::twflavor(int &CT,int &KATT0,int &KATT2,double E,double T){  

  double RTEg;
  double RTEg1;
  double RTEg2;
  double RTEg3;
  double RTEq;
  double RTEq3;
  double RTEq4;
  double RTEq5;
  double RTEq6;
  double RTEq7;
  double RTEq8;

  int vb[7]={0};
  int b=0;
  int KATT00=KATT0;
  int KATT20=KATT2;
	  
  vb[1]=1;
  vb[2]=2; 
  vb[3]=3; 
  vb[4]=-1; 
  vb[5]=-2; 
  vb[6]=-3; 

  twlinear(KATT0,E,T,RTEg1,RTEg2,RTEq6,RTEq7,RTEq8);

  //.....for gluon
  if(KATT00==21)
	{
	  //	R0  =RTE
	  double R1  =RTEg1;
	  double R2  =RTEg2;
	  //	R3  =RTEg3

	  if(KATT20==21)
		{
		  double a=ran0(&NUM1);
		  if(a<=R1/(R1+R2))
			{
			  CT=1;
			  //	        KATT3=KATT2
			  //		    KATT2=21 
			  //		    KATT0=21
			}
	     
		  if(a>R1/(R1+R2))
			{
			  CT=2;
			  //	        KATT3=KATT2 
			  b=floor(ran0(&NUM1)*6+1);
			  if(b==7)
				{
				  b=6;
				}
			  KATT2=vb[b]; 
			  KATT0=-KATT2;
			}
		}	  
         
	  if(KATT20!=21)
		{ 
	      CT=3;
		  //	      KATT3=KATT2 
		  //		  KATT2=KATT2 
		  //		  KATT0=21
		}	     
	}   	            
	  
  //.....for quark and antiquark
  if(KATT00!=21)
	{      

	  //	R00 =RTE
	  //	R3  =RTEq3
	  //	R4  =RTEq4
	  //	R5  =RTEq5
	  double R6  =RTEq6;
	  double R7  =RTEq7;
	  double R8  =RTEq8;
	  double R00 =R6+R7+R8;

	  if(KATT20==21)
		{
	      CT=13;
		  //	      KATT3=KATT2 
		  //		  KATT2=21 
		  //		  KATT0=KATT0
		}
	  	  
	  if(KATT20!=21)
	    {

	      if(abs(KATT20)!=abs(KATT00))
			{
			  CT=4;
			  //	      KATT3=KATT2 
			  //		  KATT2=KATT3 
			  //		  KATT0=KATT0
			}
	  	  
	      if(KATT20==KATT00)
			{
			  CT=5;
			  //	      KATT3=KATT2
			  //	      KATT0=KATT0
			}
   	      
	      if(KATT20==-KATT00)
			{
			  double a=ran0(&NUM1);
			  if(a<=(R6)/R00)
				{
				  CT=6;
				  //	         KATT3=KATT2
				tf2:	     b=floor(ran0(&NUM1)*3+1);
				  if(b==4)
					{
					  b=3;
					}  
				  KATT2=-KATT0/abs(KATT0)*vb[b]; 
				  if(abs(KATT2)==abs(KATT0))
					{
					  goto tf2;			
					} 
				  KATT0=-KATT2;	     
				} 
	   	   
			  if(a>(R6)/R00 && a<=(R6+R7)/R00)
				{
				  CT=7;
				  //	         KATT3=KATT2 
				  //		     KATT2=KATT3 
				  //		     KATT0=KATT0
				}	   
	   
			  if(a>(R6+R7)/R00 && a<=1.0)
				{
				  CT=8;
				  //	         KATT3=KATT2 
				  KATT2=21; 
				  KATT0=21;
				}	      
			}	   	  	   
	    } 
    }      	

}  

//.........................................................................

void LBTclass::twlinear(int KATT,double E,double T,double &RTEg1,double &RTEg2,double &RTEq6,double &RTEq7,double &RTEq8){ 

  //.....    
  double dtemp=0.02;
  int iT1=floor((T-0.1)/dtemp);
  int iT2=iT1+1;
  int iE1=floor(log(E)+2);
  int iE2=iE1+1;
  //
  double T1=0.12+(iT1-1)*0.02;
  double T2=T1+dtemp;
  double E1=exp(iE1-2.0);
  double E2=exp(iE2-2.0);
  //
  if(KATT==21)
	{
	  double RTE1=(Rg1[iT2][iE1]-Rg1[iT1][iE1])*(T-T1)/(T2-T1)+Rg1[iT1][iE1];
	  double RTE2=(Rg1[iT2][iE2]-Rg1[iT1][iE2])*(T-T1)/(T2-T1)+Rg1[iT1][iE2];
	  RTEg1=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	
	  RTE1=(Rg2[iT2][iE1]-Rg2[iT1][iE1])*(T-T1)/(T2-T1)+Rg2[iT1][iE1];
	  RTE2=(Rg2[iT2][iE2]-Rg2[iT1][iE2])*(T-T1)/(T2-T1)+Rg2[iT1][iE2];
	  RTEg2=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	
	
	  //	   RTE1=(Rg3[iT2][iE1]-Rg3[iT1][iE1])*(T-T1)/(T2-T1)+Rg3[iT1][iE1];
	  //	   RTE2=(Rg3[iT2][iE2]-Rg3[iT1][iE2])*(T-T1)/(T2-T1)+Rg3[iT1][iE2];
	  //	   RTEg3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;		
     
	}

  if(KATT!=21)
	{
	  //	   RTE1=(Rq3[iT2][iE1]-Rq3[iT1][iE1])*(T-T1)/(T2-T1)+Rq3[iT1][iE1];
	  //	   RTE2=(Rq3[iT2][iE2]-Rq3[iT1][iE2])*(T-T1)/(T2-T1)+Rq3[iT1][iE2];
	  //	   RTEq3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  //	   RTE1=(Rq4[iT2][iE1]-Rq4[iT1][iE1])*(T-T1)/(T2-T1)+Rq4[iT1][iE1];
	  //	   RTE2=(Rq4[iT2][iE2]-Rq4[iT1][iE2])*(T-T1)/(T2-T1)+Rq4[iT1][iE2];
	  //	   RTEq4=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1

	  //	   RTE1=(Rq5[iT2][iE1]-Rq5[iT1][iE1])*(T-T1)/(T2-T1)+Rq5[iT1][iE1];
	  //	   RTE2=(Rq5[iT2][iE2]-Rq5[iT1][iE2])*(T-T1)/(T2-T1)+Rq5[iT1][iE2];
	  //	   RTEq5=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  double RTE1=(Rq6[iT2][iE1]-Rq6[iT1][iE1])*(T-T1)/(T2-T1)+Rq6[iT1][iE1];
	  double RTE2=(Rq6[iT2][iE2]-Rq6[iT1][iE2])*(T-T1)/(T2-T1)+Rq6[iT1][iE2];
	  RTEq6=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq7[iT2][iE1]-Rq7[iT1][iE1])*(T-T1)/(T2-T1)+Rq7[iT1][iE1];
	  RTE2=(Rq7[iT2][iE2]-Rq7[iT1][iE2])*(T-T1)/(T2-T1)+Rq7[iT1][iE2];
	  RTEq7=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq8[iT2][iE1]-Rq8[iT1][iE1])*(T-T1)/(T2-T1)+Rq8[iT1][iE1];
	  RTE2=(Rq8[iT2][iE2]-Rq8[iT1][iE2])*(T-T1)/(T2-T1)+Rq8[iT1][iE2];
	  RTEq8=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	}	
}	

//.........................................................................

void LBTclass::trans(double v[4],double p[4]){	
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);	
  double ga=1.0/sqrt(1.0-vv*vv);
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]-ppar);
  p[1]=p[1]+v[1]*gavv;
  p[2]=p[2]+v[2]*gavv;
  p[3]=p[3]+v[3]*gavv;
}

void LBTclass::transback(double v[4],double p[4]){
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  double ga=1.0/sqrt(1.0-vv*vv); 
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(-ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]+ppar);
  p[1]=p[1]-v[1]*gavv;
  p[2]=p[2]-v[2]*gavv;
  p[3]=p[3]-v[3]*gavv;
}	

	
//.........................................................................
void LBTclass::colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt){
  //
  //    p0 initial jet momentum, output to final momentum
  //    p2 final thermal momentum,p3 initial termal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	

  //    transform to local comoving frame of the fluid
//  cout << endl;
//  cout << "flow  "<< v0[1] << " " << v0[2] << " " << v0[3] << " "<<" Elab " << p0[0] << endl;
  
  trans(v0,p0);
//  cout << p0[0] << " " << sqrt(qhat0ud) << endl;
  
//  cout << sqrt(pow(p0[1],2)+pow(p0[2],2)+pow(p0[3],2)) << " " << p0[1] << " " << p0[2] << " " << p0[3] << endl;

  //************************************************************
  trans(v0,p4);
  //************************************************************


  //    sample the medium parton thermal momentum in the comoving frame


  double xw;
  double razim;
  double rcos;
  double rsin;
	  
  double ss;  
  double tmin;
  double tmid;
  double tmax;
	  
  double rant;
  double tt;
	  
  double uu;	  
  double ff;
  double rank;
	  
  double mmax;
  double msq;
	  
  double f1;
  double f2;
  
  double p0ex[4]={0.0};

  //    Initial 4-momentum of jet
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	  
  
  int ic=0;

  do
	{
      do{	  
	  th:   xw=15.0*ran0(&NUM1);
		razim=2.0*pi*ran0(&NUM1);
		rcos=1.0-2.0*ran0(&NUM1);
		rsin=sqrt(1.0-rcos*rcos);
		//
		p2[0]=xw*temp;
		p2[3]=p2[0]*rcos;
		p2[1]=p2[0]*rsin*cos(razim);
		p2[2]=p2[0]*rsin*sin(razim);
	  
		f1=pow(xw,3 )/(exp(xw)-1)/1.4215;
		f2=pow(xw,3)/(exp(xw)+1)/1.2845;
		//
		//    cms energy
		//
		ss=2.0*(p0[0]*p2[0]-p0[1]*p2[1]-p0[2]*p2[2]-p0[3]*p2[3]);
	  
		//	if(ss.lt.2.d0*qhat0ud) goto 14

		tmin=qhat0ud;
		tmid=ss/2.0;
		tmax=ss-qhat0ud;
	  
		//    use (s^2+u^2)/(t+qhat0ud)^2 as scattering cross section in the
		//
		rant=ran0(&NUM1);
		tt=rant*ss;	  

//		ic+=1;
//		cout << p0[0] << "  " << p2[0] <<  endl;
//		cout << tt << "  " << ss <<  "" << qhat0ud <<endl;
//		cout << ic << endl;

	  }while((tt<qhat0ud) || (tt>(ss-qhat0ud)));

	  uu=ss-tt;	  
	  
	  if(CT==1)
		{
		  ff=f1;	  
		  mmax=4.0/pow(ss,2)*(3.0-tmin*(ss-tmin)/pow(ss,2)+(ss-tmin)*ss/pow(tmin,2)+tmin*ss/pow((ss-tmin),2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(3.0-tt*uu/pow(ss,2)+uu*ss/pow(tt,2)+tt*ss/pow(uu,2))/mmax;
		}

	
	  if(CT==2)
		{
		  ff=f1;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		}

	
	  if(CT==3)
		{
		  ff=f2;
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
			}
		  else
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		}

	
	  if(CT==13)
		{
		  ff=f1;
	  
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
			}
		  else 
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		}

	
	  if(CT==4)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(ss,2)+pow(uu,2))/pow(tt,2))/mmax;
		}

	
	  if(CT==5)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(ss,2)+pow(tmin,2))/pow((ss-tmin),2)-2.0/3.0*pow(ss,2)/tmin/(ss-tmin));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*((pow(ss,2)+pow(uu,2))/pow(tt,2)+(pow(ss,2)+pow(tt,2))/pow(uu,2)-2.0/3.0*pow(ss,2)/tt/uu))/mmax;	  
		}

	
	  if(CT==6)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+0.5);
		}

	
	  if(CT==7)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2)+2.0/3.0*pow((ss-tmin),2)/ss/tmin);
		  msq=(pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(((pow(ss,2)+pow(uu,2))/pow(tt,2))+(pow(tt,2)+pow(uu,2))/pow(ss,2)+2.0/3.0*pow(uu,2)/ss/tt)))/mmax;	  
		}

	
	  if(CT==8)
		{
		  ff=f2;	 
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);	  
		}
	
	  rank=ran0(&NUM1);
	}while(rank>(msq*ff));
	
  //
  p3[1]=p2[1]; 
  p3[2]=p2[2];
  p3[3]=p2[3];
  p3[0]=p2[0];

  //    velocity of the center-of-mass
  //   
  
  vc[1]=(p0[1]+p2[1])/(p0[0]+p2[0]);
  vc[2]=(p0[2]+p2[2])/(p0[0]+p2[0]);
  vc[3]=(p0[3]+p2[3])/(p0[0]+p2[0]);
  //
  //    transform into the cms frame
  //
  trans(vc,p0);
  trans(vc,p2);  
  //
  //    cm momentum
  //
  double pcm=p2[0];
  //
  //    sample transverse momentum transfer with respect to jet momentum
  //    in cm frame
  //
  double ranp=2.0*pi*ran0(&NUM1);
  //
  //    transverse momentum transfer
  //

  //
  if(pow(pcm,2) < pow((tt/2.0/pcm-pcm),2))
  {
  goto th;
  }  
  //
  
  qt=sqrt(pow(pcm,2)-pow((tt/2.0/pcm-pcm),2));
  double qx=qt*cos(ranp);
  double qy=qt*sin(ranp);      
  //
  //    longitudinal momentum transfer
  //
  double qpar=tt/2.0/pcm;
  //
  //    qt is perpendicular to pcm, need to rotate back to the cm frame
  //
  double upt=sqrt(p2[1]*p2[1]+p2[2]*p2[2])/p2[0];
  double upx=p2[1]/p2[0];
  double upy=p2[2]/p2[0];
  double upz=p2[3]/p2[0];
  //
  //    momentum after collision in cm frame
  // 
  p2[1]=p2[1]-qpar*upx;
  p2[2]=p2[2]-qpar*upy; 
  
  if(upt!=0.0) 
	{
      p2[1]=p2[1]+(upz*upx*qy+upy*qx)/upt;
      p2[2]=p2[2]+(upz*upy*qy-upx*qx)/upt;
	}
 s2:  p2[3]=p2[3]-qpar*upz-upt*qy;

  p0[1]=-p2[1];
  p0[2]=-p2[2];
  p0[3]=-p2[3];
  //
  //    transform from cm back to the comoving frame
  //         
  transback(vc,p2);
  transback(vc,p0);    
  
  //************************************************************
  //
  //     calculate qt in the rest frame of medium
  //
  if(p0[4]>p2[4])
	{
	  rotate(p4[1],p4[2],p4[3],p0,1);
	  qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
	  rotate(p4[1],p4[2],p4[3],p0,-1);
	}
  else
	{
	  rotate(p4[1],p4[2],p4[3],p2,1);
	  qt=sqrt(pow(p2[1],2)+pow(p2[2],2));
	  rotate(p4[1],p4[2],p4[3],p2,-1);
	}
  //************************************************************	  	  
	  
  //
  //    transform from comoving frame to the lab frame
  //
  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p3);
	  
  //************************************************************
  transback(v0,p4);
  //************************************************************



////////////////////////////////////////////////////////////////////////////...2019

/*
    if(switchphasecut == 1)
	{
    double collener = p2[0];
	
    double colltheta = acos((p2[1]*p0[1]+p2[2]*p0[2]+p2[3]*p0[3])/sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3])/sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]));
    double colltheta = acos(p0[2]/sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]));

    //double radener = p2[0];
    //double radtheta = acos((p2[1]*p0[1]+p2[2]*p0[2]+p2[3]*p0[3])/sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3])/sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]));
    //double radtheta = acos(p0[2]/sqrt(p2[1]*p2[1]+p2[2]*p2[2]+p2[3]*p2[3]));

    if(collener>collenercut)
    {
    goto th;
    }

    if(colltheta>collthetacut)
    {
    goto th;
    }
	}
*/

////////////////////////////////////////////////////////////////////////////...2019






  
}
	
void LBTclass::twcoll(int CT,double qhat0ud,double v0[4],double p0[4],double p2[4]){    
  //	
  //     p0 initial jet momentum, output to final momentum
  //     p2 final thermal momentum,p3 initial thermal energy
  //
  //     amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //     transform to local comoving frame of the fluid
       
  trans(v0,p0);
  trans(v0,p2);

  //    velocity of the center-of-mass
			

  vc[1]=(p0[1]+p2[1])/(p0[0]+p2[0]);
  vc[2]=(p0[2]+p2[2])/(p0[0]+p2[0]);
  vc[3]=(p0[3]+p2[3])/(p0[0]+p2[0]);
  //
  //    transform into the cms frame
  //

  trans(vc,p0);
  trans(vc,p2);	   

  //
  //     cm momentum
  //
  double pcm=p2[0];	   	   
  //
  //     sample transverse momentum transfer with respect to jet momentum
  //     in cm frame
  //
  //
  //     Gaussian distribution
  //
  //     qt=sqrt(-dt*qhat0ud*log(1-rant+rant*exp(-scm/(4.d0*dt*qhat0ud))))
  //
  //     static potential distribution
  //	   
  double ss=4.0*pow(pcm,2);
  //
  double tmin=qhat0ud;
  double tmid=ss/2.0;
  double tmax=ss-qhat0ud;
	   
  double rant;
  double tt;
  double uu;
  double mmax;
  double msq;
  double rank;	

  /////////////////////////////////////////////
  double ranp;
  double qt;
  double qx;
  double qy;
  double qpar;
  double upt;
  double upx;
  double upy;
  double upz;
  /////////////////////////////////////////////	   
	   
  //
  //	   CT is a variable notated different collision types.
  //
        
  do
	{	  
	tw:     rant=ran0(&NUM1);
	  tt=rant*ss;
		
	  if((tt<qhat0ud) || (tt>(ss-qhat0ud)))
		break;
	  //      if((tt<qhat0ud) || (tt>(ss-qhat0ud)))
	  //      {		
	  //		goto t59;
	  //      }		
	  uu=ss-tt;

	  //	   gg to gg
	  if(CT==1)
		{
		  //		
		  mmax=3.0-tmin*(ss-tmin)/pow(ss,2)+(ss-tmin)*ss/pow(tmin,2)+tmin*ss/pow((ss-tmin),2);
		  msq=(3.0-tt*uu/pow(ss,2)+uu*ss/pow(tt,2)+tt*ss/pow(uu,2))/mmax;
		  //
		}
	

	  //	   gg to qqbar
	  if(CT==2)
		{
		  //
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		  //
		}

	  //	   gq to gq, gqbar to gqbar
	  if(CT==3)
		{
		  //
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin);
			}
		  else
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		  //
		}
	   
	  //	   qg to qg, qbarg to qbarg
	  if(CT==13)
		{	   
		  //
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin);		
			}
		  else
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		  //	    
		}

	  //	   qiqj to qiqj, qiqjbar to qiqjbar, qibarqj to qibarqj, qibarqjbar to qibarqjbar
	  //	   for i not equal j
	  if(CT==4)
		{	  
		  //		
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2);
		  msq=(4.0/9.0*(pow(ss,2)+pow(uu,2))/pow(tt,2))/mmax;
		  //
		}

	  //	   qiqi to qiqi, qibarqibar to qibarqibar
	  if(CT==5)
		{
		  //	 
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(ss,2)+pow(tmin,2))/pow((ss-tmin),2)-2.0/3.0*pow(ss,2)/tmin/(ss-tmin);
		  msq=(4.0/9.0*((pow(ss,2)+pow(uu,2))/pow(tt,2)+(pow(ss,2)+pow(tt,2))/pow(uu,2)-2.0/3.0*pow(ss,2)/tt/uu))/mmax;
		  //	  
		}

	  //     qiqibar to qjqjbar for i not equal j
	  if(CT==6)
		{
		  //	 
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+0.5);
		  //
		}

	  //	   qiqibar to qiqibar
	  if(CT==7)
		{
		  //	 
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2)+2.0/3.0*pow((ss-tmin),2)/ss/tmin;
		  msq=(4.0/9.0*(((pow(ss,2)+pow(uu,2))/pow(tt,2))+(pow(tt,2)+pow(uu,2))/pow(ss,2)+2.0/3.0*pow(uu,2)/ss/tt))/mmax;
		  //
		}

	  //	   qqbar to gg
	  if(CT==8)
		{
		  //	 
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		  //
		}
	   
	  rank=ran0(&NUM1);
      

	}while(rank>msq);	   	  
	  	  
  ///////////////////////////////////////////////////////////////////////
  //    transverse momentum transfer

  if((tt>qhat0ud) && (tt<(ss-qhat0ud)))
	{

      ranp=2.0*pi*ran0(&NUM1);
	  //
	  //
	  //
      qt=sqrt(pow(pcm,2)-pow((tt/2.0/pcm-pcm),2));
      qx=qt*cos(ranp);
      qy=qt*sin(ranp);

	  //
	  //    longitudinal momentum transfer
	  //
      qpar=tt/2.0/pcm;
	  //
	  //    qt is perpendicular to pcm, need to rotate back to the cm frame
	  //
      upt=sqrt(p2[1]*p2[1]+p2[2]*p2[2])/p2[0];
      upx=p2[1]/p2[0];
      upy=p2[2]/p2[0];
      upz=p2[3]/p2[0];
	  //
	  //    momentum after collision in cm frame
	  //

      p2[1]=p2[1]-qpar*upx;
      p2[2]=p2[2]-qpar*upy;
      if(upt!=0.0) 
		{
		  p2[1]=p2[1]+(upz*upx*qy+upy*qx)/upt;
		  p2[2]=p2[2]+(upz*upy*qy-upx*qx)/upt;
		}
	s3:   p2[3]=p2[3]-qpar*upz-upt*qy;

      p0[1]=-p2[1];
      p0[2]=-p2[2];
      p0[3]=-p2[3];
	  
	}
  //
  //    transform from cm back to the comoving frame
  //
 t59:  transback(vc,p2);
       transback(vc,p0);
  //
  //    transform from comoving frame to the lab frame
  //

  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p3);
	  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void LBTclass::colljet23( double temp, double qhat0ud, double v0[4], double p0[4],double p2[4], double p3[4], double p4[4], double qt, int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint){
  //    p0 initial jet momentum, output to final momentum
  //    p3 initial thermal momentum
  //    p2 initial thermal momentum, output to final thermal momentum
  //    p4 radiated gluon momentum
  //    qt transverse momentum transfer in the rest frame of medium
  //    q0,ql energy and longitudinal momentum transfer
  //    i=0: 2->3 finished; i=1: can not find a gluon when nloop<=30
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
	  


  double xw,rank,yy,razim,pcm;
  double rant,ranp,tt,ss,qx,qy,qpar,upt,upx,upy,upz;
  double rcos,rsin,dsdt,q0,ql;
  double w,wkt2,dng0,wkt;

  double xwmin,xwmax,wkt2min,wkt2max,ppx,ppxmax,tauf,taufmax,dngtest,dngtestmax;
  double dt1,dt2,dt3;

  double px,py,pz,px0,py0,pz0,abc,E;
	  
  int nloop1,nloop2,nloop3;
	  
  const double pi=3.1415926;
	  
  nloop1=0;
  nloop2=0;
  nloop3=0;
  ic=0;
	  

  //    initial thermal parton momentum in lab frame
  p2[1]=p3[1];
  p2[2]=p3[2];
  p2[3]=p3[3];
  p2[0]=p3[0];

  double Elab=p0[0];

  //    transform to local comoving frame of the fluid
  trans(v0,p0);
  trans(v0,p2);
  //
  px0=p0[1];
  py0=p0[2];
  pz0=p0[3];

  //    rotate to the frame in which jet moves along z-axis
  rotate(px0,py0,pz0,p0,1);
  rotate(px0,py0,pz0,p2,1);
	  
  //
  //    w&kt of radiated gluon in the rest frame of medium and jet along z-direction
  //

if(p0[0] < 2 * sqrt(qhat0ud))
{
    ic=1;
    goto Ben30;
}

xwmin = sqrt(qhat0ud)/p0[0];
xwmax = 1 - xwmin;

xw=xwmin+ran0(&NUM1)*(xwmax-xwmin);

wkt2min = 0.0;
wkt2max = pow(xwmax*p0[0],2);

wkt2=wkt2min+ran0(&NUM1)*(wkt2max-wkt2min);

dt1=(tint-tisint)*Ejpint/Elabint/sctr;
dt2=(tint-tirint)*Ejpint/Elabint/sctr;
dt3=(tisint-tirint)*Ejpint/Elabint/sctr;

Ben20:
    if(isp==0) 
    {
        ppx=(1+pow((1-xw),2));
        ppxmax=(1+pow((1-xwmin),2));
    }
    else if (isp==1)
    {
        ppx=(1-xw)*(1+pow((1-xw),2));
        ppxmax=(1-xwmin)*(1+pow((1-xwmin),2));
    }
    else if (isp==2)
    {
        ppx=pow(1-xw+pow(xw,2), 3);
        ppxmax=pow(1-xwmin+pow(xwmin,2), 3);
    }

  w = xw * p0[0];

  rank=ran0(&NUM1);
  wkt2max = pow(w,2);
  wkt2=1/((1 - rank) * (1/wkt2min) + rank * (1/wkt2max));

  tauf=2*p0[0]*xw*(1-xw)/wkt2;
  taufmax=2*p0[0]*xwmin*(1-xwmin)/wkt2min;

  dngtestmax=2.0*CA*alphas/pi*px*qhat/pow(wkt2min,2)*pow(sin(dtlastrad/(2.0*taufmax)),2);
  dngtest=2.0*CA*alphas/pi*px*qhat/pow(wkt2,2)*pow(sin(dtlastrad/(2.0*tauf)),2);  

    if(dngtest/dngtestmax > 1.0)
    {
        cout << "in colljet23 dngtest=" << dngtest << endl;
        cout << "in colljet23 dngtestmax=" << dngtestmax << endl;
        cout << tint << "  " << tisint  << "  " << tirint << endl;
        cout << Elabint << "  " << Ejpint << "  " << p0[0] << endl;
        cout << dt1 << "  " << dt2  << "  " << dt3 << endl;
        cout << xw << "  " << wkt2 << endl;
        cout << ppx << "  " << ppxmax << endl;
        cout << tauf << "  " << taufmax << endl;
        exit(1);
    }
 
  rank=ran0(&NUM1);
  if(wkt2 > pow(w,2) || rank > dngtest/dngtestmax)
	{
	  nloop3=nloop3+1;
	  if(nloop3>1.0E6)
		{
		  ic=1;
		  goto Ben30;
        }
	  goto Ben20;
	}
	  
  ranp=2.0*pi*ran0(&NUM1);
  p4[1]=sqrt(wkt2)*cos(ranp);
  p4[2]=sqrt(wkt2)*sin(ranp);
  p4[3]=sqrt(pow(w,2)-wkt2);
  p4[0]=w;
	  
  //    energy & momentum transfer in lab frame

  ranp=2.0*pi*ran0(&NUM1);
  qx=qt*cos(ranp);
  qy=qt*sin(ranp);

	  
  abc=pow((-p4[3]+p2[3]+p0[3]),2)*(pow(p4[1],4)+pow(p4[2],4)+pow(p4[3],4)+2*pow(p4[3],2)\
								   *pow(p2[1],2)+pow(p2[1],4)+2*pow(p4[3],2)*pow(p2[2],2)+2*pow(p2[1],2)*pow(p2[2],2)+\
								   pow(p2[2],4)-4*pow(p4[3],3)*p2[3]-4*p4[3]*pow(p2[1],2)*p2[3]-4*p4[3]*\
								   pow(p2[2],2)*p2[3]+6*pow(p4[3],2)*pow(p2[3],2)+2*pow(p2[1],2)*pow(p2[3],2)+2*\
								   pow(p2[2],2)*pow(p2[3],2)-4*p4[3]*pow(p2[3],3)+pow(p2[3],4)-2*pow(p4[3],2)*\
								   pow(p2[0],2)-2*pow(p2[1],2)*pow(p2[0],2)-2*pow(p2[2],2)*pow(p2[0],2)+4*p4[3]*\
								   p2[3]*pow(p2[0],2)-2*pow(p2[3],2)*pow(p2[0],2)+pow(p2[0],4)+2*pow(p4[3],2)*\
								   pow(p0[1],2)-2*pow(p2[1],2)*pow(p0[1],2)-2*pow(p2[2],2)*pow(p0[1],2)-4*p4[3]*\
								   p2[3]*pow(p0[1],2)+2*pow(p2[3],2)*pow(p0[1],2)-2*pow(p2[0],2)*pow(p0[1],2)+\
								   pow(p0[1],4)+2*pow(p4[3],2)*pow(p0[2],2)-2*pow(p2[1],2)*pow(p0[2],2)-2*pow(p2[2],2)\
								   *pow(p0[2],2)-4*p4[3]*p2[3]*pow(p0[2],2)+2*pow(p2[3],2)*pow(p0[2],2)-2*\
								   pow(p2[0],2)*pow(p0[2],2)+2*pow(p0[1],2)*pow(p0[2],2)+pow(p0[2],4)-4*pow(p4[3],3)\
								   *p0[3]-4*p4[3]*pow(p2[1],2)*p0[3]-4*p4[3]*pow(p2[2],2)*p0[3]+12*\
								   pow(p4[3],2)*p2[3]*p0[3]+4*pow(p2[1],2)*p2[3]*p0[3]+4*pow(p2[2],2)*p2[3]\
								   *p0[3]-12*p4[3]*pow(p2[3],2)*p0[3]+4*pow(p2[3],3)*p0[3]+4*p4[3]\
								   *pow(p2[0],2)*p0[3]-4*p2[3]*pow(p2[0],2)*p0[3]-4*p4[3]*pow(p0[1],2)*p0[3]\
								   +4*p2[3]*pow(p0[1],2)*p0[3]-4*p4[3]*pow(p0[2],2)*p0[3]+4*p2[3]*\
								   pow(p0[2],2)*p0[3]+6*pow(p4[3],2)*pow(p0[3],2)+2*pow(p2[1],2)*pow(p0[3],2)+2*\
								   pow(p2[2],2)*pow(p0[3],2)-12*p4[3]*p2[3]*pow(p0[3],2)+6*pow(p2[3],2)*pow(p0[3],2)\
								   -2*pow(p2[0],2)*pow(p0[3],2)+2*pow(p0[1],2)*pow(p0[3],2)+2*pow(p0[2],2)*pow(p0[3],2)-4\
								   *p4[3]*pow(p0[3],3)+4*p2[3]*pow(p0[3],3)+pow(p0[3],4)-4*pow(p4[3],2)*p2[0]\
								   *p0[0]-4*pow(p2[1],2)*p2[0]*p0[0]-4*pow(p2[2],2)*p2[0]*p0[0]+8*p4[3]*\
								   p2[3]*p2[0]*p0[0]-4*pow(p2[3],2)*p2[0]*p0[0]+4*pow(p2[0],3)*p0[0]-4*\
								   p2[0]*pow(p0[1],2)*p0[0]-4*p2[0]*pow(p0[2],2)*p0[0]+8*p4[3]*p2[0]*\
								   p0[3]*p0[0]-8*p2[3]*p2[0]*p0[3]*p0[0]-4*p2[0]*pow(p0[3],2)*p0[0]\
								   -2*pow(p4[3],2)*pow(p0[0],2)-2*pow(p2[1],2)*pow(p0[0],2)-2*pow(p2[2],2)*pow(p0[0],2)\
								   +4*p4[3]*p2[3]*pow(p0[0],2)-2*pow(p2[3],2)*pow(p0[0],2)+6*pow(p2[0],2)*\
								   pow(p0[0],2)-2*pow(p0[1],2)*pow(p0[0],2)-2*pow(p0[2],2)*pow(p0[0],2)+4*p4[3]*p0[3]\
								   *pow(p0[0],2)-4*p2[3]*p0[3]*pow(p0[0],2)-2*pow(p0[3],2)*pow(p0[0],2)+4*p2[0]*\
								   pow(p0[0],3)+pow(p0[0],4)-4*pow(p4[3],2)*p2[1]*qx-4*pow(p2[1],3)*qx-4*p2[1]*\
								   pow(p2[2],2)*qx+8*p4[3]*p2[1]*p2[3]*qx-4*p2[1]*pow(p2[3],2)*qx+4*p2[1]\
								   *pow(p2[0],2)*qx+4*pow(p4[3],2)*p0[1]*qx-4*pow(p2[1],2)*p0[1]*qx-4*\
								   pow(p2[2],2)*p0[1]*qx-8*p4[3]*p2[3]*p0[1]*qx+4*pow(p2[3],2)*p0[1]*qx\
								   -4*pow(p2[0],2)*p0[1]*qx+4*p2[1]*pow(p0[1],2)*qx+4*pow(p0[1],3)*qx+4*p2[1]\
								   *pow(p0[2],2)*qx+4*p0[1]*pow(p0[2],2)*qx+8*p4[3]*p2[1]*p0[3]*qx-8*p2[1]\
								   *p2[3]*p0[3]*qx-8*p4[3]*p0[1]*p0[3]*qx+8*p2[3]*p0[1]*p0[3]*qx\
								   -4*p2[1]*pow(p0[3],2)*qx+4*p0[1]*pow(p0[3],2)*qx+8*p2[1]*p2[0]*p0[0]*\
								   qx-8*p2[0]*p0[1]*p0[0]*qx+4*p2[1]*pow(p0[0],2)*qx-4*p0[1]*pow(p0[0],2)\
								   *qx+4*pow(p4[3],2)*pow(qx,2)+4*pow(p2[1],2)*pow(qx,2)-8*p4[3]*p2[3]*pow(qx,2)+4*\
								   pow(p2[3],2)*pow(qx,2)-4*pow(p2[0],2)*pow(qx,2)+8*p2[1]*p0[1]*pow(qx,2)+4*pow(p0[1],2)\
								   *pow(qx,2)-8*p4[3]*p0[3]*pow(qx,2)+8*p2[3]*p0[3]*pow(qx,2)+4*pow(p0[3],2)*\
								   pow(qx,2)-8*p2[0]*p0[0]*pow(qx,2)-4*pow(p0[0],2)*pow(qx,2)-4*pow(p4[1],3)*\
								   (p0[1]+qx)-4*pow(p4[3],2)*p2[2]*qy-4*pow(p2[1],2)*p2[2]*qy-4*pow(p2[2],3)\
								   *qy+8*p4[3]*p2[2]*p2[3]*qy-4*p2[2]*pow(p2[3],2)*qy+4*p2[2]*\
								   pow(p2[0],2)*qy+4*p2[2]*pow(p0[1],2)*qy+4*pow(p4[3],2)*p0[2]*qy-4*\
								   pow(p2[1],2)*p0[2]*qy-4*pow(p2[2],2)*p0[2]*qy-8*p4[3]*p2[3]*p0[2]*\
								   qy+4*pow(p2[3],2)*p0[2]*qy-4*pow(p2[0],2)*p0[2]*qy+4*pow(p0[1],2)*p0[2]\
								   *qy+4*p2[2]*pow(p0[2],2)*qy+4*pow(p0[2],3)*qy+8*p4[3]*p2[2]*p0[3]*qy\
								   -8*p2[2]*p2[3]*p0[3]*qy-8*p4[3]*p0[2]*p0[3]*qy+8*p2[3]*p0[2]\
								   *p0[3]*qy-4*p2[2]*pow(p0[3],2)*qy+4*p0[2]*pow(p0[3],2)*qy+8*p2[2]\
								   *p2[0]*p0[0]*qy-8*p2[0]*p0[2]*p0[0]*qy+4*p2[2]*pow(p0[0],2)*qy-4\
								   *p0[2]*pow(p0[0],2)*qy+8*p2[1]*p2[2]*qx*qy+8*p2[2]*p0[1]*qx*qy+8\
								   *p2[1]*p0[2]*qx*qy+8*p0[1]*p0[2]*qx*qy+4*pow(p4[3],2)*pow(qy,2)+4*\
								   pow(p2[2],2)*pow(qy,2)-8*p4[3]*p2[3]*pow(qy,2)+4*pow(p2[3],2)*pow(qy,2)-4*\
								   pow(p2[0],2)*pow(qy,2)+8*p2[2]*p0[2]*pow(qy,2)+4*pow(p0[2],2)*pow(qy,2)-8*p4[3]\
								   *p0[3]*pow(qy,2)+8*p2[3]*p0[3]*pow(qy,2)+4*pow(p0[3],2)*pow(qy,2)-8*p2[0]*\
								   p0[0]*pow(qy,2)-4*pow(p0[0],2)*pow(qy,2)-4*pow(p4[2],3)*(p0[2]+qy)+4*\
								   pow(p4[3],2)*p2[0]*p4[0]+4*pow(p2[1],2)*p2[0]*p4[0]+4*pow(p2[2],2)*\
								   p2[0]*p4[0]-8*p4[3]*p2[3]*p2[0]*p4[0]+4*pow(p2[3],2)*p2[0]*\
								   p4[0]-4*pow(p2[0],3)*p4[0]+4*p2[0]*pow(p0[1],2)*p4[0]+4*p2[0]*\
								   pow(p0[2],2)*p4[0]-8*p4[3]*p2[0]*p0[3]*p4[0]+8*p2[3]*p2[0]*\
								   p0[3]*p4[0]+4*p2[0]*pow(p0[3],2)*p4[0]+4*pow(p4[3],2)*p0[0]*p4[0]\
								   +4*pow(p2[1],2)*p0[0]*p4[0]+4*pow(p2[2],2)*p0[0]*p4[0]-8*p4[3]\
								   *p2[3]*p0[0]*p4[0]+4*pow(p2[3],2)*p0[0]*p4[0]-12*pow(p2[0],2)\
								   *p0[0]*p4[0]+4*pow(p0[1],2)*p0[0]*p4[0]+4*pow(p0[2],2)*p0[0]*p4[0]\
								   -8*p4[3]*p0[3]*p0[0]*p4[0]+8*p2[3]*p0[3]*p0[0]*p4[0]+4*\
								   pow(p0[3],2)*p0[0]*p4[0]-12*p2[0]*pow(p0[0],2)*p4[0]-4*pow(p0[0],3)*p4[0]\
								   -8*p2[1]*p2[0]*qx*p4[0]+8*p2[0]*p0[1]*qx*p4[0]-8*p2[1]*p0[0]\
								   *qx*p4[0]+8*p0[1]*p0[0]*qx*p4[0]+8*p2[0]*pow(qx,2)*p4[0]+8*p0[0]\
								   *pow(qx,2)*p4[0]-8*p2[2]*p2[0]*qy*p4[0]+8*p2[0]*p0[2]*qy*p4[0]-8\
								   *p2[2]*p0[0]*qy*p4[0]+8*p0[2]*p0[0]*qy*p4[0]+8*p2[0]*pow(qy,2)*\
								   p4[0]+8*p0[0]*pow(qy,2)*p4[0]-2*pow(p4[3],2)*pow(p4[0],2)-2*pow(p2[1],2)*\
								   pow(p4[0],2)-2*pow(p2[2],2)*pow(p4[0],2)+4*p4[3]*p2[3]*pow(p4[0],2)-2*\
								   pow(p2[3],2)*pow(p4[0],2)+6*pow(p2[0],2)*pow(p4[0],2)-2*pow(p0[1],2)*pow(p4[0],2)\
								   -2*pow(p0[2],2)*pow(p4[0],2)+4*p4[3]*p0[3]*pow(p4[0],2)-4*p2[3]*p0[3]\
								   *pow(p4[0],2)-2*pow(p0[3],2)*pow(p4[0],2)+12*p2[0]*p0[0]*pow(p4[0],2)+6*\
								   pow(p0[0],2)*pow(p4[0],2)+4*p2[1]*qx*pow(p4[0],2)-4*p0[1]*qx*pow(p4[0],2)\
								   -4*pow(qx,2)*pow(p4[0],2)+4*p2[2]*qy*pow(p4[0],2)-4*p0[2]*qy*pow(p4[0],2)\
								   -4*pow(qy,2)*pow(p4[0],2)-4*p2[0]*pow(p4[0],3)-4*p0[0]*pow(p4[0],3)+\
								   pow(p4[0],4)-4*p4[2]*(p0[2]+qy)*(pow(p4[3],2)-pow(p2[1],2)-pow(p2[2],2)\
																	+pow(p2[3],2)-pow(p2[0],2)+pow(p0[1],2)+pow(p0[2],2)+2*p2[3]*p0[3]+\
																	pow(p0[3],2)-2*p4[3]*(p2[3]+p0[3])-2*p2[0]*p0[0]-pow(p0[0],2)+2\
																	*p2[1]*qx+2*p0[1]*qx+2*p2[2]*qy+2*p0[2]*qy+2*p2[0]*p4[0]\
																	+2*p0[0]*p4[0]-pow(p4[0],2))+2*pow(p4[2],2)*(pow(p4[3],2)-pow(p2[1],2)-\
																												 pow(p2[2],2)+pow(p2[3],2)-pow(p2[0],2)+pow(p0[1],2)+3*pow(p0[2],2)+2*p2[3]*p0[3]\
																												 +pow(p0[3],2)-2*p4[3]*(p2[3]+p0[3])-2*p2[0]*p0[0]-pow(p0[0],2)+2*p2[1]\
																												 *qx+2*p0[1]*qx+2*p2[2]*qy+6*p0[2]*qy+2*pow(qy,2)+2*p2[0]*p4[0]+2\
																												 *p0[0]*p4[0]-pow(p4[0],2))-4*p4[1]*(p0[1]+qx)*(pow(p4[2],2)+pow(p4[3],2)\
																																								-pow(p2[1],2)-pow(p2[2],2)+pow(p2[3],2)-pow(p2[0],2)+pow(p0[1],2)+pow(p0[2],2)+2*p2[3]\
																																								*p0[3]+pow(p0[3],2)-2*p4[3]*(p2[3]+p0[3])-2*p2[0]*p0[0]-pow(p0[0],2)+2\
																																								*p2[1]*qx+2*p0[1]*qx+2*p2[2]*qy+2*p0[2]*qy-2*p4[2]*(p0[2]+qy)\
																																								+2*p2[0]*p4[0]+2*p0[0]*p4[0]-pow(p4[0],2))+2*pow(p4[1],2)*(pow(p4[2],2)+\
																																																						   pow(p4[3],2)-pow(p2[1],2)-pow(p2[2],2)+pow(p2[3],2)-pow(p2[0],2)+3*pow(p0[1],2)+\
																																																						   pow(p0[2],2)+2*p2[3]*p0[3]+pow(p0[3],2)-2*p4[3]*(p2[3]+p0[3])-2*\
																																																						   p2[0]*p0[0]-pow(p0[0],2)+2*p2[1]*qx+6*p0[1]*qx+2*pow(qx,2)+2*p2[2]*\
																																																						   qy+2*p0[2]*qy-2*p4[2]*(p0[2]+qy)+2*p2[0]*p4[0]+2*p0[0]\
																																																						   *p4[0]-pow(p4[0],2)));
  //
  if(abc<0.0)
	{
	  nloop1=nloop1+1;
	  nloop3=0;
	  if(nloop1>1.0E4)
		{
		  ic=1;
		  goto Ben30;
		}
	  goto Ben20;
	}   	  
  //

  q0=(-pow(p4[1],2)*p2[0]-pow(p4[2],2)*p2[0]+pow(p4[3],2)*p2[0]+pow(p2[1],2)*p2[0]\
	  +pow(p2[2],2)*p2[0]-2*p4[3]*p2[3]*p2[0]+pow(p2[3],2)*p2[0]-pow(p2[0],3)\
	  +2*p4[1]*p2[0]*p0[1]-p2[0]*pow(p0[1],2)+2*p4[2]*p2[0]*p0[2]-p2[0]\
	  *pow(p0[2],2)-2*p4[3]*p2[0]*p0[3]+2*p2[3]*p2[0]*p0[3]+p2[0]*\
	  pow(p0[3],2)-pow(p4[1],2)*p0[0]-pow(p4[2],2)*p0[0]-pow(p4[3],2)*p0[0]+pow(p2[1],2)\
	  *p0[0]+pow(p2[2],2)*p0[0]+2*p4[3]*p2[3]*p0[0]-pow(p2[3],2)*p0[0]-\
	  pow(p2[0],2)*p0[0]+2*p4[1]*p0[1]*p0[0]-pow(p0[1],2)*p0[0]+2*p4[2]*p0[2]\
	  *p0[0]-pow(p0[2],2)*p0[0]+2*p4[3]*p0[3]*p0[0]-2*p2[3]*p0[3]*p0[0]-\
	  pow(p0[3],2)*p0[0]+p2[0]*pow(p0[0],2)+pow(p0[0],3)+2*p4[1]*p2[0]*qx-2*p2[1]\
	  *p2[0]*qx-2*p2[0]*p0[1]*qx+2*p4[1]*p0[0]*qx-2*p2[1]*p0[0]*qx-2\
	  *p0[1]*p0[0]*qx+2*p4[2]*p2[0]*qy-2*p2[2]*p2[0]*qy-2*p2[0]*p0[2]\
	  *qy+2*p4[2]*p0[0]*qy-2*p2[2]*p0[0]*qy-2*p0[2]*p0[0]*qy+pow(p4[1],2)\
	  *p4[0]+pow(p4[2],2)*p4[0]+pow(p4[3],2)*p4[0]-pow(p2[1],2)*p4[0]-pow(p2[2],2)*\
	  p4[0]-2*p4[3]*p2[3]*p4[0]+pow(p2[3],2)*p4[0]+pow(p2[0],2)*p4[0]-2*p4[1]\
	  *p0[1]*p4[0]+pow(p0[1],2)*p4[0]-2*p4[2]*p0[2]*p4[0]+pow(p0[2],2)*p4[0]-2\
	  *p4[3]*p0[3]*p4[0]+2*p2[3]*p0[3]*p4[0]+pow(p0[3],2)*p4[0]-2*p2[0]\
	  *p0[0]*p4[0]-3*pow(p0[0],2)*p4[0]-2*p4[1]*qx*p4[0]+2*p2[1]*qx*p4[0]\
	  +2*p0[1]*qx*p4[0]-2*p4[2]*qy*p4[0]+2*p2[2]*qy*p4[0]+2*p0[2]*qy\
	  *p4[0]+p2[0]*pow(p4[0],2)+3*p0[0]*pow(p4[0],2)-pow(p4[0],3)+sqrt(abc))/(2*\
																			  (pow(p4[3],2)+pow(p2[3],2)-pow(p2[0],2)+2*p2[3]*p0[3]+pow(p0[3],2)-2*p4[3]*\
																			   (p2[3]+p0[3])-2*p2[0]*p0[0]-pow(p0[0],2)+2*p2[0]*p4[0]+2\
																			   *p0[0]*p4[0]-pow(p4[0],2)));
  //

  ql=(pow(p4[1],2)*pow(p4[3],2)+pow(p4[2],2)*pow(p4[3],2)+pow(p4[3],4)-pow(p4[3],2)*\
	  pow(p2[1],2)-pow(p4[3],2)*pow(p2[2],2)-2*pow(p4[1],2)*p4[3]*p2[3]-2*pow(p4[2],2)\
	  *p4[3]*p2[3]-2*pow(p4[3],3)*p2[3]+2*p4[3]*pow(p2[1],2)*p2[3]+2*p4[3]\
	  *pow(p2[2],2)*p2[3]+pow(p4[1],2)*pow(p2[3],2)+pow(p4[2],2)*pow(p2[3],2)-pow(p2[1],2)\
	  *pow(p2[3],2)-pow(p2[2],2)*pow(p2[3],2)+2*p4[3]*pow(p2[3],3)-pow(p2[3],4)-pow(p4[3],2)\
	  *pow(p2[0],2)+pow(p2[3],2)*pow(p2[0],2)-2*p4[1]*pow(p4[3],2)*p0[1]+4*p4[1]*\
	  p4[3]*p2[3]*p0[1]-2*p4[1]*pow(p2[3],2)*p0[1]+pow(p4[3],2)*pow(p0[1],2)-2\
	  *p4[3]*p2[3]*pow(p0[1],2)+pow(p2[3],2)*pow(p0[1],2)-2*p4[2]*pow(p4[3],2)*p0[2]\
	  +4*p4[2]*p4[3]*p2[3]*p0[2]-2*p4[2]*pow(p2[3],2)*p0[2]+pow(p4[3],2)\
	  *pow(p0[2],2)-2*p4[3]*p2[3]*pow(p0[2],2)+pow(p2[3],2)*pow(p0[2],2)-2*pow(p4[1],2)\
	  *p4[3]*p0[3]-2*pow(p4[2],2)*p4[3]*p0[3]-4*pow(p4[3],3)*p0[3]+2*p4[3]\
	  *pow(p2[1],2)*p0[3]+2*p4[3]*pow(p2[2],2)*p0[3]+2*pow(p4[1],2)*p2[3]*p0[3]\
	  +2*pow(p4[2],2)*p2[3]*p0[3]+6*pow(p4[3],2)*p2[3]*p0[3]-2*pow(p2[1],2)*p2[3]\
	  *p0[3]-2*pow(p2[2],2)*p2[3]*p0[3]-2*pow(p2[3],3)*p0[3]+2*p4[3]*pow(p2[0],2)\
	  *p0[3]+4*p4[1]*p4[3]*p0[1]*p0[3]-4*p4[1]*p2[3]*p0[1]*p0[3]-2*\
	  p4[3]*pow(p0[1],2)*p0[3]+2*p2[3]*pow(p0[1],2)*p0[3]+4*p4[2]*p4[3]*p0[2]\
	  *p0[3]-4*p4[2]*p2[3]*p0[2]*p0[3]-2*p4[3]*pow(p0[2],2)*p0[3]+2*p2[3]\
	  *pow(p0[2],2)*p0[3]+pow(p4[1],2)*pow(p0[3],2)+pow(p4[2],2)*pow(p0[3],2)+6*pow(p4[3],2)*\
	  pow(p0[3],2)-pow(p2[1],2)*pow(p0[3],2)-pow(p2[2],2)*pow(p0[3],2)-6*p4[3]*p2[3]*\
	  pow(p0[3],2)-pow(p2[0],2)*pow(p0[3],2)-2*p4[1]*p0[1]*pow(p0[3],2)+pow(p0[1],2)*\
	  pow(p0[3],2)-2*p4[2]*p0[2]*pow(p0[3],2)+pow(p0[2],2)*pow(p0[3],2)-4*p4[3]*\
	  pow(p0[3],3)+2*p2[3]*pow(p0[3],3)+pow(p0[3],4)-2*pow(p4[3],2)*p2[0]*p0[0]+2*\
	  pow(p2[3],2)*p2[0]*p0[0]+4*p4[3]*p2[0]\
	  *p0[3]*p0[0]-2*p2[0]*pow(p0[3],2)\
	  *p0[0]-pow(p4[3],2)*pow(p0[0],2)+pow(p2[3],2)*pow(p0[0],2)+2*p4[3]*p0[3]*\
	  pow(p0[0],2)-pow(p0[3],2)*pow(p0[0],2)-2*p4[1]*pow(p4[3],2)*qx\
	  +2*pow(p4[3],2)*p2[1]\
	  *qx+4*p4[1]*p4[3]*p2[3]*qx-4*p4[3]*p2[1]*p2[3]*qx-2*p4[1]*\
	  pow(p2[3],2)*qx+2*p2[1]*pow(p2[3],2)*qx+2*pow(p4[3],2)*p0[1]*qx-4*p4[3]\
	  *p2[3]*p0[1]*qx+2*pow(p2[3],2)*p0[1]*qx+4*p4[1]*p4[3]*p0[3]*qx-4\
	  *p4[3]*p2[1]*p0[3]*qx-4*p4[1]*p2[3]*\
	  p0[3]*qx+4*p2[1]*p2[3]*p0[3]\
	  *qx-4*p4[3]*p0[1]*p0[3]*qx+4*p2[3]*p0[1]*p0[3]*qx-2*p4[1]*\
	  pow(p0[3],2)*qx+2*p2[1]*pow(p0[3],2)*qx+2*p0[1]*pow(p0[3],2)*qx-2*p4[2]*\
	  pow(p4[3],2)*qy+2*pow(p4[3],2)*p2[2]*qy+4*p4[2]*p4[3]*p2[3]*qy-4*p4[3]\
	  *p2[2]*p2[3]*qy-2*p4[2]*pow(p2[3],2)*qy+2*p2[2]*pow(p2[3],2)*qy+2*\
	  pow(p4[3],2)*p0[2]*qy-4*p4[3]*p2[3]*p0[2]*qy+2*pow(p2[3],2)*p0[2]*qy+4\
	  *p4[2]*p4[3]*p0[3]*qy-4*p4[3]*p2[2]*p0[3]*qy\
	  -4*p4[2]*p2[3]*p0[3]\
	  *qy+4*p2[2]*p2[3]*p0[3]*qy-4*p4[3]*p0[2]*p0[3]*qy+4*p2[3]*p0[2]\
	  *p0[3]*qy-2*p4[2]*pow(p0[3],2)*qy+2*p2[2]*pow(p0[3],2)*qy+2*p0[2]\
	  *pow(p0[3],2)*qy+2*pow(p4[3],2)*p2[0]*p4[0]-2*pow(p2[3],2)*p2[0]*p4[0]-\
	  4*p4[3]*p2[0]*p0[3]*p4[0]+2*p2[0]*pow(p0[3],2)*p4[0]+2*pow(p4[3],2)*\
	  p0[0]*p4[0]-2*pow(p2[3],2)*p0[0]*p4[0]-4*p4[3]*p0[3]*p0[0]*p4[0]\
	  +2*pow(p0[3],2)*p0[0]*p4[0]-pow(p4[3],2)*pow(p4[0],2)+pow(p2[3],2)*pow(p4[0],2)+2\
	  *p4[3]*p0[3]*pow(p4[0],2)-pow(p0[3],2)*pow(p4[0],2)-p2[0]*sqrt(abc)-p0[0]*\
	  sqrt(abc)+p4[0]*sqrt(abc))/(2*(p4[3]-p2[3]-p0[3])*(pow(p4[3],2)\
														 +pow(p2[3],2)-pow(p2[0],2)+2*p2[3]*p0[3]+pow(p0[3],2)-2*p4[3]*(p2[3]\
																														+p0[3])-2*p2[0]*p0[0]-pow(p0[0],2)+2*p2[0]*p4[0]+2\
														 *p0[0]*p4[0]-pow(p4[0],2)));
  //

  if(p0[0]+q0-p4[0]>0 && p2[0]-q0>0)
	{
	  
      p0[1]=p0[1]+qx-p4[1];
      p0[2]=p0[2]+qy-p4[2];
      p0[3]=p0[3]+ql-p4[3];
      p0[0]=p0[0]+q0-p4[0];

      p2[1]=p2[1]-qx;
      p2[2]=p2[2]-qy;
      p2[3]=p2[3]-ql;
      p2[0]=p2[0]-q0;
	}
  else
	{
	  nloop2=nloop2+1;
	  nloop1=0;
	  nloop3=0;
	  if(nloop2>1.0E4)
		{
          ic=1;
          goto Ben30;
        }
	  goto Ben20;
	}

  //    rotate back to lab frame
  rotate(px0,py0,pz0,p0,-1);
  rotate(px0,py0,pz0,p2,-1);
  rotate(px0,py0,pz0,p4,-1);

  //    transform from comoving frame to the lab frame
  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p4);

 Ben30:    double e=2.718;  
 
}


void LBTclass::radiation(double qhat0ud,double v0[4],double P2[4],double P3[4],double P4[4],double Pj0[4], int &ic, double tint, double tisint, double tirint, double dtlastrad, double Elabint, double Ejpint){
  //    work in the rest frame of medium
  //    return the 4-momentum of final states in 1->3 radiation
  //    input: P2(4)-momentum of radiated gluon from 2->3
  //           P3(4)-momentum of daughter parton from 2->3
  //           Pj0(4)-inital momentum of jet before 2->3
  //           v0-local velocity of medium
  //    output:P2(4)-momentum of 1st radiated gluon
  //           P3(4)-momentum of daughter parton
  //           P4(4)-momentum of 2nd radiated gluon
  //           i=1: no radiation; 0:successful radiation

	  
  double rcos,rsin,rcoskt,rsinkt,rank;
  double px0,py0,pz0,E0,px,py,pz;
  double xw,wkt2;
  double xwmin,xwmax,wkt2min,wkt2max,ppx,ppxmax,tauf,taufmax,dngtest,dngtestmax;
  double dt1,dt2,dt3;

  double P50P,Q2,z,w,z1,pt,dltkt,akt,akx,aky;
  double P3m,P2p,P2m,P4p,P4m,P5p,P5m,cospt,sinpt;

  int nloop1,nloop2,nloop3,i;	  

  double P50[4]={0.0};
  double P2i[4]={0.0};
  double P3i[4]={0.0};	  
	  
  ic=0;
  nloop1=0;
  nloop2=0;

  for(int i=0;i<=3;i++)
	{
	  P2i[i]=P2[i];
	  P3i[i]=P3[i];
	}
	  
  //    transform to local comoving frame of the fluid

  trans(v0,P2);
  trans(v0,P3);
  trans(v0,Pj0);

  //
  px0=P2[1]+P3[1];
  py0=P2[2]+P3[2];
  pz0=P2[3]+P3[3];


  //    rotate to the frame in which jet moves along z-axis
  rotate(px0,py0,pz0,P3,1);
  rotate(px0,py0,pz0,P2,1);
  rotate(px0,py0,pz0,Pj0,1);
  
  for(int i=0;i<=3;i++)
	{
	  P50[i]=P2[i]+P3[i];
	}
	  
  //     transfer to light cone momentum

  P50P=P50[0]+P50[3];
  Q2=P50[0]*P50[0]-P50[3]*P50[3];  //! Q^2--virtuality of parton 1
  z =(P2[0]+P2[3])/P50P;

  //    cos of angle beteen kt and ptp
  //    rcos is sampled according to (1-rcos)^2 distribution

  rcos=1.0-2.0*ran0(&NUM1);
  rsin=sqrt(1-rcos*rcos);
  //
  //     w&kt of radiated gluon in the rest frame of medium and jet along z-direction

if(P50[0]*(1-z) < 2 * sqrt(qhat0ud))
{
    ic=1;
    goto Rad30;
}


xwmin = sqrt(qhat0ud)/Pj0[0];
xwmax = 1 - xwmin;

xw=xwmin+ran0(&NUM1)*(xwmax-xwmin);

wkt2min = 0.0;
wkt2max = pow(xwmax*Pj0[0],2);

wkt2=wkt2min+ran0(&NUM1)*(wkt2max-wkt2min);

dt1=(tint-tisint)*Ejpint/Elabint/sctr;
dt2=(tint-tirint)*Ejpint/Elabint/sctr;
dt3=(tisint-tirint)*Ejpint/Elabint/sctr;

Rad20:
    if(isp==0) 
    {
        ppx=(1+pow((1-xw),2));
        ppxmax=(1+pow((1-xwmin),2));
    }
    else if (isp==1)
    {
        ppx=(1-xw)*(1+pow((1-xw),2));
        ppxmax=(1-xwmin)*(1+pow((1-xwmin),2));
    }
    else if (isp==2)
    {
        ppx=pow(1-xw+pow(xw,2), 3);
        ppxmax=pow(1-xwmin+pow(xwmin,2), 3);
    }

  w = xw * Pj0[0];
  if(xw >= 1) goto Rad20;

  tauf=2*Pj0[0]*xw*(1-xw)/wkt2;
  taufmax=2*Pj0[0]*xwmin*(1-xwmin)/wkt2min;

  dngtestmax=2.0*CA*alphas/pi*px*qhat/pow(wkt2min,2)*pow(sin(dtlastrad/(2.0*taufmax)),2);
  dngtest=2.0*CA*alphas/pi*px*qhat/pow(wkt2,2)*pow(sin(dtlastrad/(2.0*tauf)),2); 




    if(dngtest/dngtestmax > 1.0)
    {
        cout << "in radiation dngtest=" << dngtest << endl;
        cout << "in radiation dngtestmax=" << dngtestmax << endl;
        cout << tint << "  " << tisint  << "  " << tirint << endl;
        cout << Elabint << "  " << Ejpint << "  " << Pj0[0] << endl;
        cout << dt1 << "  " << dt2  << "  " << dt3 << endl;
        cout << "isp = " << isp << endl;
        cout << xw << "  " << wkt2 << endl;
        cout << ppx << "  " << ppxmax << endl;
        cout << tauf << "  " << taufmax << endl;

        exit(1);
    }
 
	  
  rank=ran0(&NUM1);
  if(wkt2 > pow(w,2) || rank > dngtest/dngtestmax)
	{
	  nloop2=nloop2+1;
	  if(nloop2>1.0E6)
		{
		  ic=1;
		  goto Rad30;
        }
	  goto Rad20;
	}
	  
  z1=xw*(Pj0[0]+Pj0[3])/(P50P*(1-z));

if(z1>=1)
{
    goto Rad20;
}

  pt=sqrt(wkt2);
	  
  dltkt=(pow((4*rcos*pt*z-8*rcos*pt*z*z1),2)-4*(-z-4*z1+4*z*z1+4*pow(z1,2)-4*z*pow(z1,2))*(-4*pow(pt,2)*z+4*Q2*z*z1-4*Q2*pow(z,2)*z1-4*Q2*z*z1*z1+4*Q2*pow(z,2)*pow(z1,2)));
  if(dltkt>=0.0)
	{
	  akt=max( (-4*rcos*pt*z + 8*rcos*pt*z*z1 +sqrt(dltkt))/(2*(-z - 4*z1 + 4*z*z1 + 4*pow(z1,2) - 4*z*pow(z1,2))),(-4*rcos*pt*z + 8*rcos*pt*z*z1 -sqrt(dltkt))/(2*(-z - 4*z1 + 4*z*z1 + 4*pow(z1,2) - 4*z*pow(z1,2))) );
	 
  if(akt<0)
	{
	  goto Rad20;
	}
    }
  else
	{
	  nloop1=nloop1+1;
	  if(nloop1>1.0E4)
		{
		  ic=1;
		  goto Rad30;
        }
	  goto Rad20;
	}
    
	  
  //     decide the direction of kt on transverse plane

  rcoskt=1.0-2.0*ran0(&NUM1);
  rsinkt=sqrt(1-rcoskt*rcoskt);
  akx=akt*rcoskt;
  aky=akt*rsinkt;
  px=(rcos*akx - aky*rsin)*pt/akt;
  py=(rcos*aky + akx*rsin)*pt/akt;
  //
  P3m=(pt*pt+akt*akt/4+(2*z1-1)*pt*akt*rcos)/(z1*(1-z1)*(1-z)*P50P);
  P2p=z*P50P;
  P2m=Q2/P50P-P3m;
  P4p=z1*(1-z)*P50P;
  P4m=(pt*pt+akt*akt/4-pt*akt*rcos)/P4p;
  P5p=(1-z1)*(1-z)*P50P;
  P5m=(pt*pt+akt*akt/4+pt*akt*rcos)/P5p;
  cospt=rcos*rcoskt-rsin*rsinkt;
  sinpt=rsin*rcoskt+rcos*rsinkt;
  
  P2[1]=akx;
  P2[2]=aky;
  P2[3]=(P2p-P2m)/2;
  P2[0]=(P2p+P2m)/2;
  P4[1]=pt*cospt-akx/2;
  P4[2]=pt*sinpt-aky/2;
  P4[3]=(P4p-P4m)/2;
  P4[0]=(P4p+P4m)/2;
  P3[1]=-pt*cospt-akx/2;
  P3[2]=-pt*sinpt-aky/2;
  P3[3]=(P5p-P5m)/2;
  P3[0]=(P5p+P5m)/2;

  rotate(px0,py0,pz0,P4,-1);
  transback(v0,P4);
	  
 Rad30:rotate(px0,py0,pz0,P3,-1);
  rotate(px0,py0,pz0,P2,-1);
  rotate(px0,py0,pz0,Pj0,-1);
	  
  //     transfer back to lab frame

  transback(v0,P2);
  transback(v0,P3);
  transback(v0,Pj0);

      if(P2[0] < sqrt(qhat0ud) || P2[0] > Pj0[0] - sqrt(qhat0ud)  ||  P3[0] < sqrt(qhat0ud) || P3[0] > Pj0[0] - sqrt(qhat0ud)  ||  P4[0] < sqrt(qhat0ud) || P4[0] > Pj0[0] - sqrt(qhat0ud))
      {
          ic=1;
          for(int i=0;i<=3;i++)
        	{
              P4[i]=0;
        	  P2[i]=P2i[i];
        	  P3[i]=P3i[i];
        	}
      }

}

	  
	  
void LBTclass::rotate(double px,double py,double pz,double pr[4],int icc){
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

  if(pt==0)
	{
	  cosa=1;
	  sina=0;
	} 
  else
	{
	  cosa=px/pt;
	  sina=py/pt;
	}

  cosb=pz/E;
  sinb=pt/E;

  if(icc==1)
	{ 
	  wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
	  wy1=-wx*sina+wy*cosa;
	  wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
	}  

  else
	{
	  wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
	  wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
	  wz1=-wx*sinb+wz*cosb;
	}

  wx=wx1;
  wy=wy1;
  wz=wz1;

  pr[1]=wx;
  pr[2]=wy;
  pr[3]=wz;      

  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}
	  
	  
void LBTclass::dngintt(double dtint,double tint,double tisint,double tirint,double dtlastrad,double qhat0int,double Elabint,double Ejpint){
  //....table of dng0(100,100) dn/dxw/dkt2, integrate t\in [tirad(k),ti]
  //    xw \in [mu/Ejp,(Ejp-mu)/Ej], kt2 \in [0,wmax]
  	  
//  double vc0int[4]={0.0};
  //	  double xwm[3]={0.0};
  double xwint,wkt2int; 
  double dngmax;
	  
  xwm[1]=sqrt(qhat0int)/Ejpint;
  xwm[2]=(Ejpint-sqrt(qhat0int))/Ejpint;
  wkt2m=pow((Ejpint-sqrt(qhat0int)),2);
	  
  for(int nxw=1;nxw<=100;nxw++)
	{
	  for(int nkt2=1;nkt2<=100;nkt2++)
		{
		  dng0[nxw][nkt2]=0.0;
		}
	}
	  
	  for(int nxw=1;nxw<=100;nxw++)		 
		{   
		  xwint=xwm[1]+(nxw-0.5)*(xwm[2]-xwm[1])/100.0;
		  for(int nkt2=1;nkt2<=100;nkt2++)
		    {
			  wkt2int=(nkt2-0.5)*wkt2m/100.0;
			  dng0[nxw][nkt2]=dng0[nxw][nkt2]+dngrad(xwint,wkt2int,Elabint,Ejpint,tint,tisint,tirint,qhat0int,dtlastrad);     //the order of Ejp and Elab
            }
		}		 

  //    normalization of dng0

  dngmax=0.0;

  for(int nxw=1;nxw<=100;nxw++)
	{
	  for(int nkt2=1;nkt2<=100;nkt2++)
		{		 
		  dngmax=max(dngmax,dng0[nxw][nkt2]);
		}
	}
	  

  for(int nxw=1;nxw<=100;nxw++)
	{
	  for(int nkt2=1;nkt2<=100;nkt2++)
		{
		  dng0[nxw][nkt2]=dng0[nxw][nkt2]/dngmax;
		}
	}
      
}
	  
	  
	  
int LBTclass::KPoisson(double alambda){
  //....Generate numbers according to Poisson distribution
  //    P(lambda)=lambda**k exp(-lambda)/k!
  //    input: average number of radiated gluons lambda
  //    output: number of radiated gluons in this event

  double target,p;		

  double KKPoisson=0;
  target=exp(-alambda);
  p=ran0(&NUM1);
		
  while(p>target)
	{
	  p=p*ran0(&NUM1);
	  KKPoisson=KKPoisson+1;
	}		
  return KKPoisson;
}	  
	  
/*
double angluon(double dtint,double tint,double tisint,double tirint,double qhat0int,double Elabint,double Ejpint){
  //....number of radiated gluon in [tirad,ti]       

//  double vc0int[4]={0.0};
  //	  double Elabint,Ejpint;
  double xwint,wkt2int;

  double aangluon=0.0;		  
	  
	  for(int nxw=1;nxw<=100;nxw++)
		{
		  xwint=sqrt(qhat0)/Ejpint+(nxw-0.5)*(Ejpint-2*sqrt(qhat0))/100.0/Ejpint;
		  for(int nkt2=1;nkt2<=100;nkt2++)
			{
			  wkt2int=(nkt2-0.5)*pow((xwint*Ejpint),2)/100.0;			   

			  aangluon=aangluon+dngrad(xwint,wkt2int,Elabint,Ejpint,tint,tisint,tirint,qhat0int)*pow((xwint*Ejpint),2)/100.0*(Ejpint-2*sqrt(qhat0))/100.0/Ejpint;			   
            }
		}
	 

  return aangluon;
}	  
*/








  
double LBTclass::angluon(double dtint,double tint,double tisint,double tirint,double qhat0int,double Elabint,double Ejpint,double dtlrf,double dtlastrad){
  //....number of radiated gluon in [tirad,ti]       

//  double vc0int[4]={0.0};
  //	  double Elabint,Ejpint;
  double xwint,wkt2int;

  double aangluon=0.0;		  
	  
  int nnxw=100;
  int nnkt2=100;
	  
	  for(int nxw=1;nxw<=nnxw;nxw++)
		{
		  xwint=sqrt(qhat0)/Ejpint+(nxw-0.5)*(Ejpint-2*sqrt(qhat0))/nnxw/Ejpint;
		  for(int nkt2=1;nkt2<=nnkt2;nkt2++)
			{
			  wkt2int=(nkt2-0.5)*pow((xwint*Ejpint),2)/nnkt2;			   

			  aangluon=aangluon+dngrad(xwint,wkt2int,Elabint,Ejpint,tint,tisint,tirint,qhat0int,dtlastrad)*pow((xwint*Ejpint),2)/nnkt2*(Ejpint-2*sqrt(qhat0))/nnxw/Ejpint;			   
            }
		}
		
		aangluon=aangluon*dtlrf;

  return aangluon;
}	
	  
	  
double LBTclass::dngrad(double xw0, double wkt20, double Elab0, double Ejp0, double tint0, double tis0, double tir0, double qhat0int0,double dtlastrad){  
  //    dng/dx/dkt2 in Lab frame	  
  double str,px,tauf,delt;
  double dt1,dt2,dt3;
  double CF;
  int CA,NC;  

  double ddngrad=0.0;
  str=0.19732;
  CF=4.0/3;
  CA=3;
  NC=3;
	  
  if(isp==0)
  {
      px=(1+pow((1-xw0),2))/xw0;
  }
  else if(isp==1)
  {
      px=1.0/2.0*(1-xw0)*(1+pow((1-xw0),2))/xw0;
  }
  else if(isp==2)
  {
      px=pow(1-xw0+pow(xw0,2),3)/xw0/(1-xw0);
  }

  tauf=2*Ejp0*xw0*(1-xw0)/wkt20;
	  
  dt1=(tint0-tis0)*Ejp0/Elab0/str;
  dt2=(tint0-tir0)*Ejp0/Elab0/str;
  dt3=(tis0-tir0)*Ejp0/Elab0/str;

  if(qhat0int0<1e-6)       
	{
	  ddngrad=0.0;
	}
  else
	{
	  //......qhat global variable the value of qhat can be get now 	  
	  //ddngrad=2*CA*alphas/pi*px*qhat/pow(wkt20,2)*0.5*dt1*(1- tauf/dt1*(sin(dt2/tauf)-sin(dt3/tauf)));
	  
	  ddngrad=2.0*CA*alphas/pi*px*qhat/pow(wkt20,2)*pow(sin(dtlastrad/(2.0*tauf)),2);	  
	}
  return ddngrad;
	  
}





















	  
	  
	  
	  
	  
	  
/*	  
double dngrad(double xw0, double wkt20, double Elab0, double Ejp0, double tint0, double tis0, double tir0, double qhat0int0){  
  //    dng/dx/dkt2 in Lab frame	  
  double str,px,tauf,delt;
  double dt1,dt2,dt3;
  double CF;
  int CA,NC;  

  double ddngrad=0.0;
  str=0.19732;
  CF=4.0/3;
  CA=3;
  NC=3;
	  
  if(isp==0)
  {
      px=(1+pow((1-xw0),2))/xw0;
  }
  else if(isp==1)
  {
      px=1.0/2.0*(1-xw0)*(1+pow((1-xw0),2))/xw0;
  }
  else if(isp==2)
  {
      px=pow(1-xw0+pow(xw0,2),3)/xw0/(1-xw0);
  }

  tauf=2*Ejp0*xw0*(1-xw0)/wkt20;
	  
  dt1=(tint0-tis0)*Ejp0/Elab0/str;
  dt2=(tint0-tir0)*Ejp0/Elab0/str;
  dt3=(tis0-tir0)*Ejp0/Elab0/str;

  if(qhat0int0<1e-6)       
	{
	  ddngrad=0.0;
	}
  else
	{
	  //......qhat global variable the value of qhat can be get now 	  
	  ddngrad=2*CA*alphas/pi*px*qhat/pow(wkt20,2)*0.5*dt1*(1- tauf/dt1*(sin(dt2/tauf)-sin(dt3/tauf)));
	}
  return ddngrad;
	  
}
*/

double LBTclass::dng(double x,double y){
  //   interpolation of dn/dxdkt2 as function of xw,kt**2
  //   xw(1) - min. of xw; xw(2) - max. of xw
  //   wkt2 - maxmum of kt**2
  //   x-sampled xw; y-sampled wkt2
   
  double dxw,xwmin,dkt2,wkt2min,sl;
  int i,j;
  double wi[4]={0.0};
  double wj[4]={0.0};
	   
  dxw=(xwm[2]-xwm[1])/100.0;
  xwmin=xwm[1];
  dkt2=wkt2m/100.0;
  wkt2min=0.0;

  sl=(x-xwmin)/dxw+1;
  i=int(sl);
  if(i<1)
	{
	  i=1;
	}
  if(i>98)
	{
	  i=98;
	}
  wi[2]=sl-i;
  wi[3]=wi[2]*(wi[2]-1.0)*0.5;
  wi[1]=1.0-wi[2]+wi[3];
  wi[2]=wi[2]-2.0*wi[3];

  sl=(y-wkt2min)/dkt2+1;
  j=int(sl);
  if(j<1)
	{
	  j=1;
	}
  if(j>98)
	{
	  j=98;
	}
  wj[2]=sl-j;
  wj[3]=wj[2]*(wj[2]-1.0)*0.5;
  wj[1]=1.0-wj[2]+wj[3];
  wj[2]=wj[2]-2.0*wj[3];

  double numg=0.0;
  
  for(int i1=1;i1<=3;i1++)
	{				
	  for(int j1=1;j1<=3;j1++)
		{
		  numg=numg+dng0[i+i1-1][j+j1-1]*wi[i1]*wj[j1];
		}
	}

  return numg;
  	   
}
//

void LBTclass::bulklinear(double tau, double x,double y,double eta,double &ed, double &temp,double &fraction,double &VX,double &VY,double &VZ)
{


  if(switchmedium==0)
  {
  double dxh=0.3;
  double dyh=0.3;
  double detah=0.3;  
  double dtauh=0.3;

  double x0h=-10.8;
  double y0h=-10.8;
  double eta0h=-4.8;  
  double tau0h=0.2;
  
 ///////////////////////////////////////////////////////////////////////////// 
 //cout<< tau <<" "<< x <<" "<< y <<" "<< eta <<" "<< ed <<" "<< temp <<" "<< fraction <<" "<< VX <<" "<< VY <<" "<< VZ << endl;
 ///////////////////////////////////////////////////////////////////////////// 
  
  
  int ix1=floor((x-x0h)/dxh)+1;
  int ix2=ix1+1;
  int iy1=floor((y-y0h)/dyh)+1;
  int iy2=iy1+1;
  int ieta1=floor((eta-eta0h)/detah)+1;
  int ieta2=ieta1+1;  
  int itau1=floor((tau-tau0h)/dtauh)+1;
  int itau2=itau1+1;  
  
  double x1=x0h+(ix1-1)*dxh;
  double x2=x1+dxh;
  double y1=y0h+(iy1-1)*dyh;
  double y2=y1+dyh;  
  double eta1=eta0h+(ieta1-1)*detah;
  double eta2=eta1+detah;  
  double tau1=tau0h+(itau1-1)*dtauh;
  double tau2=tau1+dtauh;  

  double xs1=(x2-x)/(x2-x1);
  double ys1=(y2-y)/(y2-y1);
  double etas1=(eta2-eta)/(eta2-eta1);
  double taus1=(tau2-tau)/(tau2-tau1);

  double xs2=(x-x1)/(x2-x1);
  double ys2=(y-y1)/(y2-y1);
  double etas2=(eta-eta1)/(eta2-eta1);
  double taus2=(tau-tau1)/(tau2-tau1);

 ///////////////////////////////////////////////////////////////////////////// 
 //cout<<"hahahahahahahaha "<<endl;
 /////////////////////////////////////////////////////////////////////////////  
  
  
  
  temp = xs1*ys1*etas1*taus1*temphydro[itau1][ix1][iy1][ieta1];
  temp+= xs1*ys1*etas1*taus2*temphydro[itau2][ix1][iy1][ieta1];
  temp+= xs1*ys1*etas2*taus1*temphydro[itau1][ix1][iy1][ieta2];
  temp+= xs1*ys1*etas2*taus2*temphydro[itau2][ix1][iy1][ieta2];
  temp+= xs1*ys2*etas1*taus1*temphydro[itau1][ix1][iy2][ieta1];
  temp+= xs1*ys2*etas1*taus2*temphydro[itau2][ix1][iy2][ieta1];
  temp+= xs1*ys2*etas2*taus1*temphydro[itau1][ix1][iy2][ieta2];
  temp+= xs1*ys2*etas2*taus2*temphydro[itau2][ix1][iy2][ieta2];
  temp+= xs2*ys1*etas1*taus1*temphydro[itau1][ix2][iy1][ieta1];
  temp+= xs2*ys1*etas1*taus2*temphydro[itau2][ix2][iy1][ieta1];
  temp+= xs2*ys1*etas2*taus1*temphydro[itau1][ix2][iy1][ieta2];
  temp+= xs2*ys1*etas2*taus2*temphydro[itau2][ix2][iy1][ieta2];
  temp+= xs2*ys2*etas1*taus1*temphydro[itau1][ix2][iy2][ieta1];
  temp+= xs2*ys2*etas1*taus2*temphydro[itau2][ix2][iy2][ieta1];
  temp+= xs2*ys2*etas2*taus1*temphydro[itau1][ix2][iy2][ieta2];
  temp+= xs2*ys2*etas2*taus2*temphydro[itau2][ix2][iy2][ieta2];

  fraction = xs1*ys1*etas1*taus1*fractionhydro[itau1][ix1][iy1][ieta1];
  fraction+= xs1*ys1*etas1*taus2*fractionhydro[itau2][ix1][iy1][ieta1];
  fraction+= xs1*ys1*etas2*taus1*fractionhydro[itau1][ix1][iy1][ieta2];
  fraction+= xs1*ys1*etas2*taus2*fractionhydro[itau2][ix1][iy1][ieta2];
  fraction+= xs1*ys2*etas1*taus1*fractionhydro[itau1][ix1][iy2][ieta1];
  fraction+= xs1*ys2*etas1*taus2*fractionhydro[itau2][ix1][iy2][ieta1];
  fraction+= xs1*ys2*etas2*taus1*fractionhydro[itau1][ix1][iy2][ieta2];
  fraction+= xs1*ys2*etas2*taus2*fractionhydro[itau2][ix1][iy2][ieta2];
  fraction+= xs2*ys1*etas1*taus1*fractionhydro[itau1][ix2][iy1][ieta1];
  fraction+= xs2*ys1*etas1*taus2*fractionhydro[itau2][ix2][iy1][ieta1];
  fraction+= xs2*ys1*etas2*taus1*fractionhydro[itau1][ix2][iy1][ieta2];
  fraction+= xs2*ys1*etas2*taus2*fractionhydro[itau2][ix2][iy1][ieta2];
  fraction+= xs2*ys2*etas1*taus1*fractionhydro[itau1][ix2][iy2][ieta1];
  fraction+= xs2*ys2*etas1*taus2*fractionhydro[itau2][ix2][iy2][ieta1];
  fraction+= xs2*ys2*etas2*taus1*fractionhydro[itau1][ix2][iy2][ieta2];
  fraction+= xs2*ys2*etas2*taus2*fractionhydro[itau2][ix2][iy2][ieta2];
  
  VX = xs1*ys1*etas1*taus1*VXhydro[itau1][ix1][iy1][ieta1];
  VX+= xs1*ys1*etas1*taus2*VXhydro[itau2][ix1][iy1][ieta1];
  VX+= xs1*ys1*etas2*taus1*VXhydro[itau1][ix1][iy1][ieta2];
  VX+= xs1*ys1*etas2*taus2*VXhydro[itau2][ix1][iy1][ieta2];
  VX+= xs1*ys2*etas1*taus1*VXhydro[itau1][ix1][iy2][ieta1];
  VX+= xs1*ys2*etas1*taus2*VXhydro[itau2][ix1][iy2][ieta1];
  VX+= xs1*ys2*etas2*taus1*VXhydro[itau1][ix1][iy2][ieta2];
  VX+= xs1*ys2*etas2*taus2*VXhydro[itau2][ix1][iy2][ieta2];
  VX+= xs2*ys1*etas1*taus1*VXhydro[itau1][ix2][iy1][ieta1];
  VX+= xs2*ys1*etas1*taus2*VXhydro[itau2][ix2][iy1][ieta1];
  VX+= xs2*ys1*etas2*taus1*VXhydro[itau1][ix2][iy1][ieta2];
  VX+= xs2*ys1*etas2*taus2*VXhydro[itau2][ix2][iy1][ieta2];
  VX+= xs2*ys2*etas1*taus1*VXhydro[itau1][ix2][iy2][ieta1];
  VX+= xs2*ys2*etas1*taus2*VXhydro[itau2][ix2][iy2][ieta1];
  VX+= xs2*ys2*etas2*taus1*VXhydro[itau1][ix2][iy2][ieta2];
  VX+= xs2*ys2*etas2*taus2*VXhydro[itau2][ix2][iy2][ieta2];

  VY = xs1*ys1*etas1*taus1*VYhydro[itau1][ix1][iy1][ieta1];
  VY+= xs1*ys1*etas1*taus2*VYhydro[itau2][ix1][iy1][ieta1];
  VY+= xs1*ys1*etas2*taus1*VYhydro[itau1][ix1][iy1][ieta2];
  VY+= xs1*ys1*etas2*taus2*VYhydro[itau2][ix1][iy1][ieta2];
  VY+= xs1*ys2*etas1*taus1*VYhydro[itau1][ix1][iy2][ieta1];
  VY+= xs1*ys2*etas1*taus2*VYhydro[itau2][ix1][iy2][ieta1];
  VY+= xs1*ys2*etas2*taus1*VYhydro[itau1][ix1][iy2][ieta2];
  VY+= xs1*ys2*etas2*taus2*VYhydro[itau2][ix1][iy2][ieta2];
  VY+= xs2*ys1*etas1*taus1*VYhydro[itau1][ix2][iy1][ieta1];
  VY+= xs2*ys1*etas1*taus2*VYhydro[itau2][ix2][iy1][ieta1];
  VY+= xs2*ys1*etas2*taus1*VYhydro[itau1][ix2][iy1][ieta2];
  VY+= xs2*ys1*etas2*taus2*VYhydro[itau2][ix2][iy1][ieta2];
  VY+= xs2*ys2*etas1*taus1*VYhydro[itau1][ix2][iy2][ieta1];
  VY+= xs2*ys2*etas1*taus2*VYhydro[itau2][ix2][iy2][ieta1];
  VY+= xs2*ys2*etas2*taus1*VYhydro[itau1][ix2][iy2][ieta2];
  VY+= xs2*ys2*etas2*taus2*VYhydro[itau2][ix2][iy2][ieta2];
  
  VZ = xs1*ys1*etas1*taus1*VZhydro[itau1][ix1][iy1][ieta1];
  VZ+= xs1*ys1*etas1*taus2*VZhydro[itau2][ix1][iy1][ieta1];
  VZ+= xs1*ys1*etas2*taus1*VZhydro[itau1][ix1][iy1][ieta2];
  VZ+= xs1*ys1*etas2*taus2*VZhydro[itau2][ix1][iy1][ieta2];
  VZ+= xs1*ys2*etas1*taus1*VZhydro[itau1][ix1][iy2][ieta1];
  VZ+= xs1*ys2*etas1*taus2*VZhydro[itau2][ix1][iy2][ieta1];
  VZ+= xs1*ys2*etas2*taus1*VZhydro[itau1][ix1][iy2][ieta2];
  VZ+= xs1*ys2*etas2*taus2*VZhydro[itau2][ix1][iy2][ieta2];
  VZ+= xs2*ys1*etas1*taus1*VZhydro[itau1][ix2][iy1][ieta1];
  VZ+= xs2*ys1*etas1*taus2*VZhydro[itau2][ix2][iy1][ieta1];
  VZ+= xs2*ys1*etas2*taus1*VZhydro[itau1][ix2][iy1][ieta2];
  VZ+= xs2*ys1*etas2*taus2*VZhydro[itau2][ix2][iy1][ieta2];
  VZ+= xs2*ys2*etas1*taus1*VZhydro[itau1][ix2][iy2][ieta1];
  VZ+= xs2*ys2*etas1*taus2*VZhydro[itau2][ix2][iy2][ieta1];
  VZ+= xs2*ys2*etas2*taus1*VZhydro[itau1][ix2][iy2][ieta2];
  VZ+= xs2*ys2*etas2*taus2*VZhydro[itau2][ix2][iy2][ieta2]; 
  
  }



  if(switchmedium==1)
  {
  temp=temp0medium;
  fraction=1.0;
  VX=VXmedium;
  VY=VYmedium;
  VZ=VZmedium;
  }














  
 ///////////////////////////////////////////////////////////////////////////// 
 //cout<<"ooooooooooooooooooooo"<<endl;
 /////////////////////////////////////////////////////////////////////////////   
  
  
//  cout << eta << " " << eta0h << " " << detah << " " << eta-eta0h << " " << (eta-eta0h)/detah << " " << floor((eta-eta0h)/detah) << endl;
//  cout << x << "  " << y << " " <<  eta << " " << tau << endl;
//  cout << ix1 << " " << ix2 << " " << iy1 << " " << iy2 << " " << ieta1 << " " << ieta2 << " " << itau1 << " " << itau2 << endl;
//  cout << x1 << " " << x2 << " " << y1 << " " << y2 << " " << eta1 << " " << eta2 << " " << tau1 << " " << tau2 << endl;
//  cout << xs1 << " " << xs2 << " " << ys1 << " " << ys2 << " " << etas1 << " " << etas2 << " " << taus1 << " " << taus2 << endl;

}





void LBTclass::radiationHQ(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt){  

  //    work in the rest frame of medium
  //    return the 4-momentum of final states in 1->3 radiation
  //    input: P2(4)-momentum of radiated gluon from 2->3
  //           P3(4)-momentum of daughter parton from 2->3
  //           Pj0(4)-inital momentum of jet before 2->3
  //           v0-local velocity of medium
  //    output:P2(4)-momentum of 1st radiated gluon
  //           P3(4)-momentum of daughter parton
  //           P4(4)-momentum of 2nd radiated gluon
  //           i=1: no radiation; 0:successful radiation


  double randomX,randomY;
  int count_sample,nloopOut,flagOut;
  double theta_gluon,kperp_gluon;
  double kpGluon[4];
  double HQmass=sqrt(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]);

  double P50[4]={0.0};
  double P2i[4]={0.0};
  double P3i[4]={0.0};	  
	  
  ic=0;

  for(int i=0;i<=3;i++){
      P2i[i]=P2[i];
      P3i[i]=P3[i];
  }
	  
  // transform to local comoving frame of the fluid
  trans(v0,P2);
  trans(v0,P3);
  trans(v0,Pj0);

  xInt=(P3[0]-HQmass)/HQenergy-xLow;
  if(xInt<=0.0) {
      ic=1;
      return;
  }

  double px0=P2[1]+P3[1];
  double py0=P2[2]+P3[2];
  double pz0=P2[3]+P3[3];

  // rotate to the frame in which jet moves along z-axis
  rotate(px0,py0,pz0,P3,1);
  rotate(px0,py0,pz0,P2,1);
  rotate(px0,py0,pz0,Pj0,1);
  
  for(int i=0;i<=3;i++){
      P50[i]=P2[i]+P3[i];
  }
	
  nloopOut=0;
  flagOut=0;
  // sample the second gluon
  do {
      do {
          randomX=xLow+xInt*ran0(&NUM1);
          randomY=ran0(&NUM1);
      } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
     
      count_sample=0;
      while(max_Ng*ran0(&NUM1)>dNg_over_dxdydt(parID,randomX,randomY,HQenergy,HQmass,temp_med,Tdiff)) {
          count_sample=count_sample+1;
          if(count_sample>1e+6) {
              cout << "give up loop at point 1 ..." << endl;
              kpGluon[1]=0.0;
              kpGluon[2]=0.0;
              kpGluon[3]=0.0;
              kpGluon[0]=0.0;
              ic=1;
    //          break;
              return;
          }
    
          do {
             randomX=xLow+xInt*ran0(&NUM1);
             randomY=ran0(&NUM1);
          } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
      }
    
      theta_gluon=2.0*pi*ran0(&NUM1);
      kperp_gluon=randomX*randomY*HQenergy;
      kpGluon[1]=kperp_gluon*cos(theta_gluon);
      kpGluon[2]=kperp_gluon*sin(theta_gluon);
      kpGluon[3]=randomX*HQenergy*sqrt(1.0-randomY*randomY);
      kpGluon[0]=sqrt(kpGluon[1]*kpGluon[1]+kpGluon[2]*kpGluon[2]+kpGluon[3]*kpGluon[3]);
      
      rotate(Pj0[1],Pj0[2],Pj0[3],kpGluon,-1);
    
      if(kpGluon[0]>(P3[0]-HQmass)) { // which should be impossible due to reset of xInt above
          kpGluon[1]=0.0;
          kpGluon[2]=0.0;
          kpGluon[3]=0.0;
          kpGluon[0]=0.0;
          nloopOut++;
          continue;
      }
    
      // update mometum
    
      P4[1]=kpGluon[1];
      P4[2]=kpGluon[2];
      P4[3]=kpGluon[3];
      P4[0]=kpGluon[0];
    
    // solve energy-momentum conservation
    // my notation in the note
    // p0 is re-constructed off-shell parton, correponding to P5
    // p1 is the final heavy quark from 2->3, corresponding to P3 (input). P3 (output) is the final heavy quark after 1->3.
    // k1 is the first gluon, corresponding to P2
    // k2 is the second gluon, corresponding to P4
    // assume k10 and p10 unchanged and modify their other components while p0 and k2 are fixed
      double sp0z=sqrt(P50[1]*P50[1]+P50[2]*P50[2]+P50[3]*P50[3]);
      double sk10=P2[0];
      double sk1zOld=P2[3];
      double sk1z,sk1p,sk1x,sk1y,sktheta; // unknown
      double sp10=P3[0];
      double sk20=P4[0];
      double sk2x=P4[1];
      double sk2y=P4[2];
      double sk2z=P4[3];
      double sk2p=sqrt(sk2x*sk2x+sk2y*sk2y);
    
      double stheta12,cos12Min2;
      double sAA,aaa,bbb,ccc,abc,abc2;
      double sk1z1,sk1z2,sk1p1,sk1p2;
      int nloop1=0;
      int nloop2=0;
      int flagDone=0;
      int yesA,yesB;
      
      sAA=sk10*sk10+sk2p*sk2p+sp0z*sp0z+sk2z*sk2z-2.0*sp0z*sk2z+HQmass*HQmass-(sp10-sk20)*(sp10-sk20);
      cos12Min2=(sAA*sAA/4.0/sk10/sk10-(sp0z-sk2z)*(sp0z-sk2z))/sk2p/sk2p;
    //  cout << "cos^2 min: " << cos12Min2 << endl;
    
      if(cos12Min2>1.0) {
          nloopOut++;
          continue;
    //      cout << "cos^2 min is outside the kinematic range: " << cos12Min2 << endl;
    //      return;
      }       
    
      do {
          stheta12=2.0*pi*ran0(&NUM1); // theta between k1 and k2
          aaa=4.0*((sp0z-sk2z)*(sp0z-sk2z)+sk2p*sk2p*cos(stheta12)*cos(stheta12));
          bbb=-4.0*sAA*(sp0z-sk2z);
          ccc=sAA*sAA-4.0*sk10*sk10*sk2p*sk2p*cos(stheta12)*cos(stheta12);
    
          abc2=bbb*bbb-4.0*aaa*ccc;
    
          if(abc2<0.0) {
              nloop1++;
              continue;
          } else {
              nloop1=0;
          }
    
          abc=sqrt(abc2);
          sk1z1=(-bbb+abc)/2.0/aaa;
          sk1z2=(-bbb-abc)/2.0/aaa;
          sk1p1=sqrt(sk10*sk10-sk1z1*sk1z1);
          sk1p2=sqrt(sk10*sk10-sk1z2*sk1z2);
    
    // Since we have squared both sides of the original equation during solving, the solutions may not satisfy the original equation. Double check is needed.
    
    
    // require time-like of p1 and k10>k1z;
          if(2.0*sk1p1*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z1+sAA<0.000001 && sk10>abs(sk1z1) && sp10*sp10-(sp0z-sk1z1)*(sp0z-sk1z1)-(sk10*sk10-sk1z1*sk1z1)>HQmass*HQmass) {
              yesA=1; 
          } else {
              yesA=0;
          }
          if(2.0*sk1p2*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z2+sAA<0.000001 && sk10>abs(sk1z2) && sp10*sp10-(sp0z-sk1z2)*(sp0z-sk1z2)-(sk10*sk10-sk1z2*sk1z2)>HQmass*HQmass) {
              yesB=1; 
          } else {
              yesB=0;
          }
              
    // select appropriate solution
          if(yesA==0 && yesB==0) {
    //          cout << "Solutions fail ..." << endl;
              nloop2++;
              continue;
          } else if(yesA==1 && yesB==1) {
    //          cout << "Both solutions work!" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              if(abs(sk1z1-sk1zOld)<abs(sk1z2-sk1zOld)) {
                  sk1z=sk1z1;
              } else {
                  sk1z=sk1z2;
              }
          } else if(yesA==1) {
    //          cout << "pass A ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sk1z=sk1z1;
          } else {
    //          cout << "pass B ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sk1z=sk1z2;
          }
    
          flagDone=1;   
    
          sk1p=sqrt(sk10*sk10-sk1z*sk1z);
    //           cout << "check solution: " << pow(sp10-sk20,2)-pow(sk1p,2)-pow(sk2p,2)-2.0*sk1p*sk2p*cos(stheta12)-pow(sp0z-sk1z-sk2z,2)-pow(HQmass,2) <<"  "<< aaa*sk1z*sk1z+bbb*sk1z+ccc <<"  "<<2.0*sk1p*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z+sAA<<"  "<<2.0*sk1p*sk2p*cos(stheta12)+2.0*(sp0z-sk2z)*sk1z-sAA << endl;
    
      } while (flagDone==0 && nloop1<loopN && nloop2<loopN);
    
      if(flagDone==0) { 
   //       ic=1;  // no appropriate solution
   //       cout << "solution fails ...  " << nloop1 <<"  "<<nloop2<< endl;
          nloopOut++;
          continue;
      } else {
          sktheta=atan2(sk2y,sk2x); // theta for k2
          sktheta=sktheta+stheta12; // theta for k1
          sk1p=sqrt(sk10*sk10-sk1z*sk1z);
          sk1x=sk1p*cos(sktheta);
          sk1y=sk1p*sin(sktheta);
    
          P2[0]=sk10;
          P2[1]=sk1x;
          P2[2]=sk1y;
          P2[3]=sk1z;
          
          P3[0]=sp10-sk20;
          P3[1]=-sk1x-sk2x;
          P3[2]=-sk1y-sk2y;
          P3[3]=sp0z-sk1z-sk2z;
    
        // rotate back to the local coordinate
          rotate(px0,py0,pz0,P2,-1);
          rotate(px0,py0,pz0,P3,-1);
          rotate(px0,py0,pz0,P4,-1);
    
        // boost back to the global frame 
          transback(v0,P2);
          transback(v0,P3);
          transback(v0,P4);
        
        // debug: check on-shell condition of P2, P3 and P4
          if(abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3])>0.000001 || abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass)>0.000001 || abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3])>0.000001) {
              cout << "Wrong solution -- not on shell" << "  " << sk10 << "  " << sk1x << "  " << sk1y << "  " << sk1z << "  " << sk20 << "  " << sk2x << "  " << sk2y << "  " << sk2z << "  " << stheta12 << "  " << sp10 << "  " << sp10-sk20 << "  " << -sk1x-sk2x << "  " << -sk1y-sk2y << "  " << sp0z << "  " << sp0z-sk1z-sk2z << "  " <<HQmass<< "  "<<pow(sp10-sk20,2)-pow(sk1x+sk2x,2)-pow(sk1y+sk2y,2)-pow(sp0z-sk1z-sk2z,2)-pow(HQmass,2)<<endl;
              cout << abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3]) <<"  "<< abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass) <<"  "<< abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3]) <<endl;
          }
    
    //       cout<<"in radiation HQ -- final gluon1, gluon2 & HQ: "<<P2[0]<<"  "<<P4[0]<<"  "<<P3[0]<<endl;
    
          flagOut=1;

          // SC: check energy-momentum conservation
          for(int i=0; i<=3; i++) {
              if(ic==0 && abs(P2i[i]+P3i[i]-P2[i]-P3[i]-P4[i])>0.000001) {
                  cout << "Warning: Violation of E.M. conservation!  " << i << " " << abs(P2i[i]+P3i[i]-P2[i]-P3[i]-P4[i]) << endl;
              }
          }
        
          // SC: check on-shell condition
          double shell2=abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3]);
          double shell3=abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass);
          double shell4=abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3]);
          if(ic==0 && (shell2>0.000001 || shell3>0.000001 || shell4>0.000001)) {
              cout << "Warning: Violation of on-shell: " << shell2 << "  " << shell3 << "  " << shell4 << endl;
          }

     }

  } while(flagOut==0 && nloopOut<loopN);

  if(flagOut==0) ic=1;


////////////////////////////////////////////////////////////////////////////...2019

/*
    double collener = ;
    double colltheta = ;

    double radener = ;
    double radtheta = ;

    if(collener>collenercut)
    {
    goto th;
    }

    if(colltheta>collthetacut)
    {
    goto th;
    }

*/

////////////////////////////////////////////////////////////////////////////...2019


 
}



void LBTclass::collHQ22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt){
  //
  //    HQ 2->2 scatterings
  //    p0 initial HQ momentum, output to final momentum
  //    p2 final thermal momentum, p3 initial thermal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************


  // transform to local comoving frame of the fluid
  trans(v0,p0);

  //************************************************************

  //    sample the medium parton thermal momentum in the comoving frame

  double xw;
  double razim;
  double rcos;
  double rsin;
	  
  double ss;  
	  
  double rant;
  double tt;
	  
  double uu;	  
  double ff;
  double rank;
	  
  double msq;
	  
  double e2,theta2,theta4,phi24;   // the four independent variables      
  double e1,e4,p1,cosTheta24,downFactor,sigFactor; // other useful variables
  double HQmass,fBmax,fFmax,fB,fF,maxValue;
  int index_p1,index_T,index_e2;
  int ct1_loop,ct2_loop,flag1,flag2;

  flag1=0;
  flag2=0;

// continue this function for HQ scattering

  HQmass=p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3];
  if(HQmass>1e-12) {
      HQmass=sqrt(HQmass);
  } else {
      HQmass = 0.0;
  }

  //    Initial 4-momentum of HQ
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	  

  p1=sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]);
  index_p1=(int)((p1-min_p1)/bin_p1);
  index_T=(int)((temp-min_T)/bin_T);
  if(index_p1>=N_p1) {
      index_p1=N_p1-1;
      cout << "warning: p1 is over p_max: " << p1 << endl;
  }
  if(index_T>=N_T) {
      index_T=N_T-1;
      cout << "warning: T is over T_max: " << temp << endl;
  } 
  if(index_T<0) {
      index_T=0;
      cout << "warning: T is below T_min: " << temp << endl;
  } 

  fBmax=distFncBM[index_T][index_p1];
  fFmax=distFncFM[index_T][index_p1];  // maximum of f(xw) at given p1 and T

  maxValue=10.0;  // need actual value later
	
  ct1_loop=0;
  do {   // sample p2 (light parton) using distribution integrated over 3 angles
        ct1_loop++;
        if(ct1_loop>1e6) {
//            cout << "cannot sample light parton for HQ scattering ..." << endl;
            flag1=1;
            break;
        }
        xw=max_e2*ran0(&NUM1);
        index_e2=(int)((xw-min_e2)/bin_e2);
        if(index_e2>=N_e2) index_e2=N_e2-1;
        if(CT==11) { // qc->qc
            ff=distFncF[index_T][index_p1][index_e2]/fFmax;
            maxValue=distMaxF[index_T][index_p1][index_e2];
        } else if(CT==12) { // gc->gc
            ff=distFncB[index_T][index_p1][index_e2]/fBmax;
            maxValue=distMaxB[index_T][index_p1][index_e2];
        } else {
            cout << "Wrong HQ channel ID" << endl;
            exit(EXIT_FAILURE);
        }
  } while(ran0(&NUM1)>ff);
      
  e2=xw*temp;
  e1=p0[0];

  // now e2 is fixed, need to sample the remaining 3 variables
  ct2_loop=0;
  do {
          ct2_loop++;
          if(ct2_loop>1e6) {
              cout << "cannot sample final states for HQ scattering ..." << endl;
              flag2=1;
              break;
          }

          theta2=pi*ran0(&NUM1);
          theta4=pi*ran0(&NUM1);
          phi24=2.0*pi*ran0(&NUM1);

          cosTheta24=sin(theta2)*sin(theta4)*cos(phi24)+cos(theta2)*cos(theta4);
          downFactor=e1-p1*cos(theta4)+e2-e2*cosTheta24;
          e4=(e1*e2-p1*e2*cos(theta2))/downFactor;
          sigFactor=sin(theta2)*sin(theta4)*e2*e4/downFactor; 

          // calculate s,t,u, different definition from light quark -- tt, uu are negative
          ss=2.0*e1*e2+HQmass*HQmass-2.0*p1*e2*cos(theta2);
          tt=-2.0*e2*e4*(1.0-cosTheta24);
          uu=2.0*HQmass*HQmass-ss-tt;

          // re-sample if the kinematic cuts are not satisfied
          if(ss<=2.0*qhat0ud || tt>=-qhat0ud || uu>=-qhat0ud) {
	       rank=ran0(&NUM1);
               sigFactor=0.0;
               msq=0.0;
               continue;
          }

	  if(CT==11) {  // qc->qc
               ff=(1.0/(exp(e2/temp)+1.0))*(1.0-1.0/(exp(e4/temp)+1.0));
               sigFactor=sigFactor*ff;
 	       msq=Mqc2qc(ss,tt,HQmass)/maxValue;
	  }

          if(CT==12) {  // gc->gc
               ff=(1.0/(exp(e2/temp)-1.0))*(1.0+1.0/(exp(e4/temp)-1.0));
               sigFactor=sigFactor*ff;
 	       msq=Mgc2gc(ss,tt,HQmass)/maxValue;
	  }

	  rank=ran0(&NUM1);

  } while(rank>(msq*sigFactor));
	
  if(flag1==0 && flag2==0) {

      // pass p2 value to p3 for initial thermal parton
      p3[1]=e2*sin(theta2);
      p3[2]=0.0;
      p3[3]=e2*cos(theta2);
      p3[0]=e2;
    
      // calculate momenta of outgoing particles
      // here p2 is for p4 (light parton) in my note
    
      p2[1]=e4*sin(theta4)*cos(phi24);
      p2[2]=e4*sin(theta4)*sin(phi24);
      p2[3]=e4*cos(theta4);
      p2[0]=e4;
    
      // Because we treated p0 (p1 in my note for heavy quark) as the z-direction, proper rotations are necessary here
      rotate(p4[1],p4[2],p4[3],p2,-1);
      rotate(p4[1],p4[2],p4[3],p3,-1);
	
      p0[1]=p4[1]+p3[1]-p2[1];
      p0[2]=p4[2]+p3[2]-p2[2];
      p0[3]=p4[3]+p3[3]-p2[3];
      p0[0]=sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]+HQmass*HQmass);
    
      // Debug
      if(fabs(p0[0]+p2[0]-p3[0]-p4[0])>0.00001) {
         cout << "Violation of energy conservation in HQ 2->2 scattering:  " << fabs(p0[0]+p2[0]-p3[0]-p4[0]) << endl;
      }
      
      // calculate qt in the rest frame of medium
      rotate(p4[1],p4[2],p4[3],p0,1);
      qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
      rotate(p4[1],p4[2],p4[3],p0,-1);

      // transform from comoving frame to the lab frame
      transback(v0,p2);
      transback(v0,p0);
      transback(v0,p3);
      transback(v0,p4);

  } else { // no scattering
      transback(v0,p0);
      transback(v0,p4);
      qt=0;
      p2[0]=0;
      p2[1]=0;
      p2[2]=0;
      p2[3]=0;
      p3[0]=0;
      p3[1]=0;
      p3[2]=0;
      p3[3]=0;
  }

	  	 	  
}



  double LBTclass::Mqc2qc(double s, double t, double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;
      double MM;

      MM=64.0/9.0*(pow((m2m-u),2)+pow((s-m2m),2)+2.0*m2m*t)/t/t;

      return(MM);

  }



  double LBTclass::Mgc2gc(double s, double t, double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;
      double MM;

      MM=32.0*(s-m2m)*(m2m-u)/t/t;
      MM=MM+64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(s+m2m))/pow((s-m2m),2);
      MM=MM+64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(u+m2m))/pow((u-m2m),2);
      MM=MM+16.0/9.0*m2m*(4.0*m2m-t)/((s-m2m)*(m2m-u));
      MM=MM+16.0*((s-m2m)*(m2m-u)+m2m*(s-u))/(t*(s-m2m));
      MM=MM+16.0*((s-m2m)*(m2m-u)-m2m*(s-u))/(t*(u-m2m));

      return(MM);

  }


  
 void LBTclass::collHQ23(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int &ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt) {

  //    p0 initial jet momentum, output to final momentum
  //    p3 initial thermal momentum
  //    p2 initial thermal momentum, output to final thermal momentum
  //    p4 radiated gluon momentum
  //    qt transverse momentum transfer in the rest frame of medium
  //    q0,ql energy and longitudinal momentum transfer
  //    i=0: 2->3 finished; i=1: can not find a gluon when nloop<=30
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
	
  double randomX,randomY;
  int count_sample,nloopOut,flagOut;
  double theta_gluon,kperp_gluon;
  double kpGluon[4];
  double HQmass=sqrt(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]);

  ic=0;
  nloopOut=0;
  flagOut=0;

  //    initial thermal parton momentum in lab frame
  p2[1]=p3[1];
  p2[2]=p3[2];
  p2[3]=p3[3];
  p2[0]=p3[0];

//  cout<<"check1 p0[0]: "<<p0[0]<<endl;

  //    transform to local comoving frame of the fluid
  trans(v0,p0);
  trans(v0,p2);
 
//  cout<<"check2 p0[0]: "<<p0[0]<<endl;

  do {

      do {
          randomX=xLow+xInt*ran0(&NUM1);
          randomY=ran0(&NUM1);
      } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
     
      count_sample=0;
      while(max_Ng*ran0(&NUM1)>dNg_over_dxdydt(parID,randomX,randomY,HQenergy,HQmass,temp_med,Tdiff)) {
          count_sample=count_sample+1;
          if(count_sample>1e+6) {
              cout << "give up loop at point 1 ..." << endl;
              kpGluon[1]=0.0;
              kpGluon[2]=0.0;
              kpGluon[3]=0.0;
              kpGluon[0]=0.0;
              ic=1;
    //          break;
              return;
          }
    
          do {
             randomX=xLow+xInt*ran0(&NUM1);
             randomY=ran0(&NUM1);
          } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
      }
    
      theta_gluon=2.0*pi*ran0(&NUM1);
      kperp_gluon=randomX*randomY*HQenergy;
      kpGluon[1]=kperp_gluon*cos(theta_gluon);
      kpGluon[2]=kperp_gluon*sin(theta_gluon);
      kpGluon[3]=randomX*HQenergy*sqrt(1.0-randomY*randomY);
      kpGluon[0]=sqrt(kpGluon[1]*kpGluon[1]+kpGluon[2]*kpGluon[2]+kpGluon[3]*kpGluon[3]);
      
      rotate(p4[1],p4[2],p4[3],kpGluon,-1);
    
      if(kpGluon[0]>(HQenergy-HQmass)) {
          kpGluon[1]=0.0;
          kpGluon[2]=0.0;
          kpGluon[3]=0.0;
          kpGluon[0]=0.0;
          nloopOut++;
          continue;
      }
    
      // update mometum
    
      p4[1]=kpGluon[1];
      p4[2]=kpGluon[2];
      p4[3]=kpGluon[3];
      p4[0]=kpGluon[0];
    
    // solve energy-momentum conservation
      double sE1=p0[0];
      double sp1x=p0[1];
      double sp1y=p0[2];
      double sp1z=p0[3];
      double sE2=p2[0];
      double sp2x=p2[1];
      double sp2y=p2[2];
      double sp2z=p2[3];
      double sk0=p4[0];
      double skx=p4[1];
      double sky=p4[2];
      double skz=p4[3];
      double sqx,sqy,sqzA,sq0A,sqzB,sq0B,sqz,sq0,sqtheta;
      double sAA,sBB,sCC,aaa,bbb,ccc,abc,abc2;
      int nloop1=0;
      int nloop2=0;
      int flagDone=0;
      int yesA,yesB;
     
      do {
          sqtheta=2.0*pi*ran0(&NUM1);
          sqx=qt*cos(sqtheta);
          sqy=qt*sin(sqtheta);
          sAA=(sE1+sE2-sk0)/(sp1z+sp2z-skz);
          sBB=(pow(sE2,2)-pow((sE1-sk0),2)+pow((sp1x-sqx-skx),2)+pow((sp1y-sqy-sky),2)+pow((sp1z-skz),2)+pow(HQmass,2)-pow((sp2x+sqx),2)-pow((sp2y+sqy),2)-pow(sp2z,2))/2.0/(sp1z+sp2z-skz);
          aaa=sAA*sAA-1.0;
          bbb=2.0*(sAA*sp2z+sAA*sBB-sE2);
          ccc=pow((sp2x+sqx),2)+pow((sp2y+sqy),2)+sp2z*sp2z+2.0*sp2z*sBB+sBB*sBB-sE2*sE2;
          abc2=bbb*bbb-4.0*aaa*ccc;
     
          if(abc2<0.0) {
              nloop1++;
              continue;
          } else {
              nloop1=0;
          }
     
          abc=sqrt(abc2);
          sq0A=(-bbb+abc)/2.0/aaa;
          sq0B=(-bbb-abc)/2.0/aaa;
          sqzA=sAA*sq0A+sBB;
          sqzB=sAA*sq0B+sBB;
    
    // require space-like and E_final > M;
          if(sq0A*sq0A-sqx*sqx-sqy*sqy-sqzA*sqzA<0 && sE1-sq0A-sk0>HQmass && sE2+sq0A>0) {
              yesA=1; 
          } else {
              yesA=0;
          }
          if(sq0B*sq0B-sqx*sqx-sqy*sqy-sqzB*sqzB<0 && sE1-sq0B-sk0>HQmass && sE2+sq0B>0) {
              yesB=1; 
          } else {
              yesB=0;
          }
                   
    // select appropriate solution
          if(yesA==0 && yesB==0) {
    //          cout << "Solutions fail ..." << endl;
              nloop2++;
              continue;
          } else if(yesA==1 && yesB==1) {
    //          cout << "Both solutions work!" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              if(abs(sq0A)<abs(sq0B)) {
                  sq0=sq0A;
                  sqz=sqzA;
              } else {
                  sq0=sq0B;
                  sqz=sqzB;
              }
          } else if(yesA==1) {
    //          cout << "pass A ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sq0=sq0A;
              sqz=sqzA;
          } else {
    //          cout << "pass B ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sq0=sq0B;
              sqz=sqzB;
          }
    
          flagDone=1;   
      } while (flagDone==0 && nloop1<loopN && nloop2<loopN);

      if(flagDone==0) { 
    //      ic=1;  // no appropriate solution
    //      cout << "solution fails ..." << endl;
          nloopOut++;
          continue;
      } else {
          p0[0]=sE1-sq0-sk0;
          p0[1]=sp1x-sqx-skx;
          p0[2]=sp1y-sqy-sky;
          p0[3]=sp1z-sqz-skz;
     
          p2[0]=sE2+sq0;
          p2[1]=sp2x+sqx;
          p2[2]=sp2y+sqy;
          p2[3]=sp2z+sqz;
        
    //           cout<<"sE1,sE2,p0[0],p2[0],p4[0]: "<<sE1<<"  "<<sE2<<"  "<<p0[0]<<"  "<<p2[0]<<"  "<<p4[0]<<endl;
    
          transback(v0,p0);
          transback(v0,p2);
          transback(v0,p4);
        
    // debug: check on-shell condition
          if(abs(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]-HQmass*HQmass)>0.000001 || abs(p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3])>0.000001) {
              cout << "Wrong solution -- not on shell" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sq0 << "  " << sqx << "  " << sqy << "  " << sqz << endl;
              cout << abs(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]-HQmass*HQmass) << "  " << abs(p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3]) << "  " << HQmass << endl;
          }

//  cout<<"in colljet23 -- init thermal, final jet, final thermal, gluon: "<<p3[0]<<"  "<<p0[0]<<"  "<<p2[0]<<"  "<<p4[0]<<endl;
          flagOut=1;
      }
  } while(flagOut==0 && nloopOut<loopN);

  if(flagOut==0) ic=1;   

////////////////////////////////////////////////////////////////////////////...2019

/*
    double collener = ;
    double colltheta = ;

    double radener = ;
    double radtheta = ;

    if(collener>collenercut)
    {
    goto th;
    }

    if(colltheta>collthetacut)
    {
    goto th;
    }
*/

////////////////////////////////////////////////////////////////////////////...2019


}



  double LBTclass::dNg_over_dxdydt(int parID, double x0g, double y0g, double HQenergy, double HQmass, double temp_med, double Tdiff) {

      double resultFnc,tauFnc,qhatFnc;


      qhatFnc = qhat_over_T3*pow(temp_med,3);   // no longer have CF factor, taken out in splittingP too.
      tauFnc = tau_f(x0g,y0g,HQenergy,HQmass);

      resultFnc = 4.0/pi*CA*alphasHQ(x0g*y0g*HQenergy,temp_med)*splittingP(parID,x0g)*qhatFnc*pow(y0g,5)*pow(sin(Tdiff/2.0/tauFnc/sctr),2)*pow((HQenergy*HQenergy/(y0g*y0g*HQenergy*HQenergy+HQmass*HQmass)),4)/x0g/x0g/HQenergy/HQenergy/sctr;

      return(resultFnc);
  }



  double LBTclass::tau_f(double x0g, double y0g, double HQenergy, double HQmass) {

      double resultFnc;

      resultFnc = 2.0*HQenergy*x0g*(1.0-x0g)/(pow((x0g*y0g*HQenergy),2)+x0g*x0g*HQmass*HQmass);

      return(resultFnc);
  }

  
  
  
  double LBTclass::splittingP(int parID, double z0g) {

      double resultFnc;

//      resultFnc = (2.0-2.0*z0g+z0g*z0g)*CF/z0g;

      if(parID==21) resultFnc = 2.0*pow(1.0-z0g+pow(z0g,2),3)/z0g/(1.0-z0g);
      else resultFnc = (1.0-z0g)*(2.0-2.0*z0g+z0g*z0g)/z0g;

      return(resultFnc);
  }



  double LBTclass::lambdas(double kTFnc) {

      double resultFnc;
      double cMass=1.27;
      double bMass=4.19;

      if (kTFnc<cMass) {
         resultFnc=0.2;
      } else if (kTFnc<bMass) {
         resultFnc=0.172508;
      } else {
         resultFnc=0.130719;
      }

      return(resultFnc);
  }




  double LBTclass::nflavor(double kTFnc) {

      double resultFnc;
      double cMass=1.27;
      double bMass=4.19;

      if (kTFnc<cMass) {
         resultFnc=3.0;
      } else if (kTFnc<bMass) {
         resultFnc=4.0;
      } else {
         resultFnc=5.0;
      }

      return(resultFnc);
  }



  double LBTclass::alphasHQ(double kTFnc, double tempFnc) {

//      double kTEff,error_para,resultFnc;
//
//      error_para=1.0;
//
//      if(kTFnc<pi*tempFnc*error_para) {
//         kTEff=pi*tempFnc*error_para;
//      } else {
//         kTEff=kTFnc;
//      }
//
//      resultFnc=4.0*pi/(11.0-2.0*nflavor(kTEff)/3.0)/2.0/log(kTEff/lambdas(kTEff));
//
//      return(resultFnc);
      return(alphas);
  }




double LBTclass::nHQgluon(int parID,double dtLRF,double &time_gluon,double &temp_med,double &HQenergy,double &max_Ng){
  // gluon radiation probability for heavy quark       

  int time_num,temp_num,HQenergy_num;
  double delta_Ng;

  if(time_gluon>t_max) {
     cout << "accumulated time exceeds t_max" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     time_gluon=t_max;	 
  }

  if(temp_med>temp_max) {
     cout << "temperature exceeds temp_max" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     temp_med=temp_max;
	 
	 //exit(1);
	 
	 
  }

  if(HQenergy>HQener_max) {
     cout << "HQenergy exceeds HQener_max" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     HQenergy=HQener_max;
  }

  if(temp_med<temp_min) {
     cout << "temperature drops below temp_min" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     temp_med=temp_min;
  }

  time_num=(int)(time_gluon/delta_tg+0.5)+1;
  temp_num=(int)((temp_med-temp_min)/delta_temp+0.5);
  HQenergy_num=(int)(HQenergy/delta_HQener+0.5);

  if(parID==21) {
      delta_Ng=dNg_over_dt_g[time_num][temp_num][HQenergy_num];
      max_Ng=max_dNgfnc_g[time_num][temp_num][HQenergy_num];
  } else if(abs(parID)==4) {
      delta_Ng=dNg_over_dt_c[time_num][temp_num][HQenergy_num];
      max_Ng=max_dNgfnc_c[time_num][temp_num][HQenergy_num];
  } else {
      delta_Ng=dNg_over_dt_q[time_num][temp_num][HQenergy_num];
      max_Ng=max_dNgfnc_q[time_num][temp_num][HQenergy_num];
  } 
 
///////////////////////////////////////////////////////////////// 
//  cout<<delta_Ng<<" "<<D2piT<<" "<<dtLRF<<endl;
///////////////////////////////////////////////////////////////// 
 
  delta_Ng=delta_Ng*6.0/D2piT*dtLRF;
  max_Ng=max_Ng*6.0/D2piT;

//  if(delta_Ng>1) {
//     cout << "Warning: Ng greater than 1   " << time_gluon << "  " << delta_Ng << endl;
//  }

  return delta_Ng;

//  write(6,*) temp_num,HQenergy_num,delta_Ng

}	  




  void LBTclass::read_xyMC(int &numXY) {

      numXY=0;

      ifstream fxyMC("readindatafile/mc_glauber.dat");
      if(!fxyMC.is_open()) {
          cout<<"Erro openning date fxyMC!\n";
          exit (EXIT_FAILURE);
      }

      while(true) {
          if(fxyMC.eof()) {
              numXY--;
              break;
          }
          if(numXY>=maxMC) break;
          fxyMC>>initMCX[numXY]>>initMCY[numXY];
          numXY++;
      }

      cout << "Number of (x,y) points from MC-Glauber: " << numXY << endl;

      fxyMC.close();

  }

  void LBTclass::jetInitialize(int numXY) {

      int index_xy;
      double ipx,ipy,ipz,ip0,iwt;
      double pT_len,ipT,phi,rapidity;
      
      for(int i=1; i<=nj; i=i+2) {

         int index_xy = (int)(ran0(&NUM1)*numXY);
         if(index_xy>=numXY) index_xy = numXY-1;

         V[1][i]=initMCX[index_xy]; 
         V[2][i]=initMCY[index_xy];
         V[3][i]=0.0;
                
         V[0][i]=-log(1.0-ran0(&NUM1));

         pT_len=ipTmax-ipTmin;
         ipT=ipTmin+ran0(&NUM1)*pT_len;

         phi=ran0(&NUM1)*2.0*pi;
         ipx=ipT*cos(phi);
         ipy=ipT*sin(phi);
   
         rapidity=2.0*eta_cut*ran0(&NUM1)-eta_cut;
         ipz=sqrt(ipT*ipT+amss*amss)*sinh(rapidity);
         ip0=sqrt(ipT*ipT+ipz*ipz+amss*amss);

         KATT1[i]=Kjet;
         P[1][i]=ipx;
         P[2][i]=ipy;
         P[3][i]=ipz;
         P[0][i]=ip0;
         P[5][i]=ipT;
         WT[i]=pT_len;

         // let jet partons stream freely for tau_0 in x-y plane;
         V[1][i]=V[1][i]+P[1][i]/P[0][i]*tau0;
         V[2][i]=V[2][i]+P[2][i]/P[0][i]*tau0;
                                                 
         for(int j=0;j<=3;j++) Prad[j][i]=P[j][i];

         Vfrozen[0][i]=tau0;
         Vfrozen[1][i]=V[1][i];
         Vfrozen[2][i]=V[2][i];
         Vfrozen[3][i]=V[3][i];
         Tfrozen[i]=0.0;
         vcfrozen[1][i]=0.0;
         vcfrozen[2][i]=0.0;
         vcfrozen[3][i]=0.0;

         // now its anti-particle (if consider pair production)
	 KATT1[i+1]=-KATT1[i];
         P[1][i+1]=-P[1][i];
         P[2][i+1]=-P[2][i];
         P[3][i+1]=-P[3][i];
         P[0][i+1]=P[0][i];
         P[5][i+1]=P[5][i];
	 WT[i+1]=WT[i];
         V[1][i+1]=V[1][i];
         V[2][i+1]=V[2][i];
         V[3][i+1]=V[3][i];
         V[0][i+1]=V[0][i];
         for(int j=0;j<=3;j++) Prad[j][i+1]=P[j][i];
         Vfrozen[0][i+1]=Vfrozen[0][i];
         Vfrozen[1][i+1]=Vfrozen[1][i];
         Vfrozen[2][i+1]=Vfrozen[2][i];
         Vfrozen[3][i+1]=Vfrozen[3][i];
         Tfrozen[i+1]=Tfrozen[i];
         vcfrozen[1][i+1]=vcfrozen[1][i];
         vcfrozen[2][i+1]=vcfrozen[2][i];
         vcfrozen[3][i+1]=vcfrozen[3][i];

      } 

  }



void LBTclass::LBTinitialize(){	  				






//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



for(int i=0;i<50;i++)
{
for(int j=0;j<100;j++)
{
Rg[i][j]=0.0;         //total gluon scattering rate as functions of initial energy and temperature 
Rg1[i][j]=0.0;        //gg-gg              CT1
Rg2[i][j]=0.0;        //gg-qqbar           CT2
Rg3[i][j]=0.0;        //gq-qg              CT3
Rq[i][j]=0.0;         //total gluon scattering rate as functions of initial energy and temperature
Rq3[i][j]=0.0;        //qg-qg              CT13
Rq4[i][j]=0.0;        //qiqj-qiqj          CT4
Rq5[i][j]=0.0;        //qiqi-qiqi          CT5
Rq6[i][j]=0.0;        //qiqibar-qjqjbar    CT6
Rq7[i][j]=0.0;        //qiqibar-qiqibar    CT7
Rq8[i][j]=0.0;        //qqbar-gg           CT8

qhatLQ[i][j]=0.0;
qhatG[i][j]=0.0;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



RHQ[i][j]=0.0;        //total scattering rate for heavy quark
RHQ11[i][j]=0.0;      //Qq->Qq
RHQ12[i][j]=0.0;      //Qg->Qg
qhatHQ[i][j]=0.0;     //qhat of heavy quark      
}
}




//for(int i=0;i<50000;i++)
for(int i=0;i<50000;i++)
{

tirad[i]=0.0;
tiscatter[i]=0.0;
tiform[i]=0.0;	 //pythia 



Tint_lrf[i]=0.0;          //for heavy quark
radng[i]=0.0;	



for(int j=0;j<4;j++) 
{
V[j][i]=0.0;           //parton position 
P[j][i]=0.0;           //parton 4-momentum
V0[j][i]=0.0;          //negative parton position
P0[j][i]=0.0;          //negative parton 4-momentum



Prad[j][i]=0.0;

Vfrozen[j][i]=0.0;     //parton final 4 coordinate
Vfrozen0[j][i]=0.0;    //negative parton final 4 coordinate 



VV[j][i]=0.0;
VV0[j][i]=0.0;
PP[j][i]=0.0;
PP0[j][i]=0.0;


vcfrozen[j][i]=0.0;


vcfrozen0[j][i]=0.0;

		
}




////////////////////////////////////////////////////////////////////////////...2019
NS[i]=0;                   //splitting rank

NR0[i]=0;                  //scattering rank
NS0[i]=0;                  //splitting rank

Mother[i]=0;               //ID of initial parton of creation
////////////////////////////////////////////////////////////////////////////...2019





NR[i]=0;                  //scattering rank 
KATT1[i]=0;               //parton flavor  
KATT10[i]=0;              //negative parton flavor

tjp[i]=0.0;

CAT[i]=0;
CAT0[i]=0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
WT[i]=0.0;
WT0[i]=0.0;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////


Tfrozen[i]=0.0;
Tfrozen0[i]=0.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////



tempHydro[i]=0.0;
vxHydro[i]=0.0;
vyHydro[i]=0.0;
vzHydro[i]=0.0;




}



for(int i=0;i<4;i++)
{
Vtemp[i]=0.0;

vf[i]=0.0;              //flow velocity in tau-eta coordinate
vfcar[i]=0.0;           //flow velocity in t-z coordinate

vp[i]=0.0;              //position of particle		
vc0[i]=0.0;             //flow velocity     

//...dimensions in subroutine colljet and twcoll
vc[i]=0.0;              
pc0[i]=0.0;
pc2[i]=0.0;
pc3[i]=0.0;
pc4[i]=0.0;		
p0[i]=0.0;
p2[i]=0.0;
p3[i]=0.0;		
p4[i]=0.0;

pc00[i]=0.0;
pc30[i]=0.0;	
		
pc01[i]=0.0;
pb[i]=0.0;
		
Pj0[i]=0.0;
		
PGm[i]=0.0;



}


/*
const int N_p1=100;
const int N_T=50;
const int N_e2=75;
*/

for(int i=0;i<50;i++)
{
for(int j=0;j<100;j++)
{
for(int k=0;k<75;k++)
{
distFncB[i][j][k]=0.0;
distFncF[i][j][k]=0.0;
distMaxB[i][j][k]=0.0;
distMaxF[i][j][k]=0.0;
}
distFncBM[i][j]=0.0;
distFncFM[i][j]=0.0;
}
}





for(int i=0;i<75+2;i++)
{
for(int j=0;j<100+1;j++)
{
for(int k=0;k<400+1;k++)
{

dNg_over_dt_c[i][j][k]=0.0;
dNg_over_dt_q[i][j][k]=0.0;
dNg_over_dt_g[i][j][k]=0.0;
max_dNgfnc_c[i][j][k]=0.0;
max_dNgfnc_q[i][j][k]=0.0;
max_dNgfnc_g[i][j][k]=0.0;

}
}
}


for(int i=0;i<101;i++)
{
for(int j=0;j<101;j++)
{
dng0[i][j]=0.0;	 //table of dn/dkperp2/dx 
}
}


//const int maxMC=2000000;
for(int i=0;i<2000000;i++)
{
initMCX[i]=0.0;
initMCY[i]=0.0;
}


for(int i=0;i<40;i++)
{
for(int j=0;j<80;j++)
{
for(int k=0;k<80;k++)
{
for(int l=0;l<80;l++)
{
temphydro[i][j][k][l]=0.0;
fractionhydro[i][j][k][l]=0.0;
VXhydro[i][j][k][l]=0.0;
VYhydro[i][j][k][l]=0.0;
VZhydro[i][j][k][l]=0.0;
}
}
}
}


for(int i=0;i<3;i++)
{

xwm[i]=0.0;

}


ntest22=0;
ntestrad=0;








//  const double  pi=3.1415926;  //conflict with fastjet  ???
//const double    CA=3.0; 
//const double    sctr=0.197;	     // 1 GeV^-1 = 0.197 fm
//const double    pi=3.1415926; 

//...input with current machine time
//...random number seed (any negative integer)
	  
//   long  NUM1 = -33;
                  
NUM1 = -33;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate


qhat_over_T3=0.0;        //qhat/T^3 for heavy quark as fnc of (T,p)
D2piT=0.0;

// for heavy quark radiation table


delta_tg=t_max/t_gn;
delta_temp=(temp_max-temp_min)/temp_gn;
delta_HQener=HQener_max/HQener_gn;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate



//....Hydro profile


//...input information
ncall=0;
nprint=0;		

nrank=0;                    //NOT NECESSARY

dt=0.0;                       //dtime when tauswitch is turned off (0) / dtau when tauswitch is turned on (1)
ti0=0.0;                      //initial time
timend=0.0;                   //end time or tau RENAME
al=0.0;                       //2 dimensional plots range NOT NECESSARY

ener=0.0;                     //initial energy of the jet parton/heavy quark
amss=0.0;                     //initial mass of the jet parton/heavy quark 
temp=0.0;                    
temp0=0.0;                  
temp00=0.0; 
Ecut=0.0;                     //energy cut of the recoiled partons
		
Kjet=0;                     //initial flavor of the jet parton
Kqhat0=0;                   //Debye screening mass switch RENAME
Kalphas=0;                  //alphas switch		
KINT=0;                     //radiation switch
Kradiation=0;                     //radiation switch
Kprimary=0;
Kfishbone=0;
tauswitch=0;                //coordinate switch

LBTswitch=0;                //switch of the Linear Boltzmann Transport

switchCoLBT_Hydro=0;        //switch of coupled LBT Hydro simulation
		
//...time-tau 				
dtau=0.0;		
tau0=0.0;		
tauend=0.0;		
time0=0.0;

//...variables with switch		
alphas=0.0;
qhat0=0.0;                    //Debye mass RENAME
qhat00=0.0;		

np=0;                          //number of partons at this time step ti				
		
//...radiation block		
icl22=0;                   
icl23=0;                    //the numerical switch in colljet23
iclrad=0;                   //the numerical switch in radiation 	  
isp=0;                      //the splitting function switch
nj=0;                       //number of leading shower partons 

//...Hydroprofile
tauhydroend=0.0;
tauhydrofile=0.0;
tauhydro0=0.0;
dtauhydro=0.0;
ntauhydro=0;
nxhydro=0;
nyhydro=0;
netahydro=0;

//...global variable qhat



	
counth100=0;
  
//double qhat;	                 //transport parameter

qhat=0.0;	                 //transport parameter


//...time system

eGluon=0.0;
nGluon=0.0;











	  

wkt2m=0.0;
	  

Ecut=0.0;                   //energy cut for free streaming   

//........................................................................................................NCUT

		
ncut=0;
ncut0=0;

//........................................................................................................NCUTEND

//////////////
//...test
		
n_coll22=0;		
n_coll23=0;
ng_coll23=0;
ng_nrad=0;
n_radiation=0;
ng_radiation=0;
n_gluon=0;
n_sp1=0;
n_sp2=0;
//      ofstream datEg("./Eg.dat");
//////////////	  



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate

// Variables for HQ 2->2
min_p1=0.0;
min_T=0.0;
min_e2=0.0;
max_p1=0.0;
max_T=0.0;
max_e2=0.0;
bin_p1=0.0;
bin_T=0.0;
bin_e2=0.0;



//int loopN=10000;
loopN=1000;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////...HQupdate







////////////////////////////////////////////////////////////////////////////...2019
switchtwcoll = 0;
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////...2019
NRcut = 0;                   //scattering rank
NScut = 0;                   //splitting rank

NR0cut = 0;                  //scattering rank for negative
NS0cut = 0;                  //splitting rank for negative

dMD = 0.0;                     //distance between mother and daughter parton
dMDcut = 0.0;                  //distance of decoherence


collenercut = 0.0;
collthetacut = 0.0;

radenercut = 0.0;
radthetacut = 0.0;
////////////////////////////////////////////////////////////////////////////...2019


Singlestepswitch = 0;
Force22 = 0;
Force23 = 0;
Force2n = 0;

ExchangeIDswitch = 0;
////////////////////////////////////////////////////////////////////////////...2019
































































//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>




//...readin LBT input

     char buf[1024];
	 
	 char readin_route[1024];
	 
	 
     ifstream f0("LBTinput.txt");
     if(!f0.is_open())
     {
      cout<<"Erro openning date file1!\n";
     }
     else
     {
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;	 	 
     f0 >> LBTswitch;
     cout << LBTswitch <<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;


	 
     f0.getline(buf, 1024);
     cout << buf << endl;	 	 
     f0 >> Kradiation;
     cout << Kradiation <<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;
	 
     f0.getline(buf, 1024);
     cout << buf << endl;	 
     f0 >> Kprimary;
     cout << Kprimary <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;	 
     f0 >> Kfishbone;
     cout << Kfishbone <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;
	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> tauswitch;
     cout << tauswitch <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> switchCoLBT_Hydro;
     cout << switchCoLBT_Hydro <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> Ecut;
     cout << Ecut <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> switchtwcoll;
     cout << switchtwcoll <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> Kqhat0;
     cout << Kqhat0 <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> alphas;
     cout << alphas <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> ncall >> nprint >> NUM1;
     cout << ncall << "   " << nprint << "      " << NUM1 <<endl;	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> time0 >> timend >> dt;
     cout << time0 << "         " << timend << "          " << dt <<endl;	
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> switchsingle;
     cout << switchsingle <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> NRcut >> NScut;
     cout << NRcut << "      " << NScut <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl;

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> switchphasecut >> collenercut >> collthetacut >> radenercut >> radthetacut;
     cout << switchphasecut << "                  " << collenercut << "             " << collthetacut << "             " << radenercut << "           " << radthetacut <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl;	

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> Singlestepswitch >> Force22 >> Force23 >> Force2n >> radlength;
     cout << Singlestepswitch << "                  " << Force22 << "             " << Force23 << "             " << Force2n << "           " << radlength <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl; 

     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> ExchangeIDswitch;
     cout << ExchangeIDswitch <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl; 
	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> ener >> mass >> Kjet >> px0 >> py0 >> pz0 >> en0 >> Xinitial >> Yinitial >> Zinitial;
     cout << ener << "      " << mass << "        " << Kjet << "       " << px0 << "       " << py0 << "     " << pz0 << "       " << en0 << "     " << Xinitial << "           " << Yinitial << "            " << Zinitial <<endl;	
     f0.getline(buf, 1024);
     cout << buf << endl;
	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> switchmedium;
     cout << switchmedium <<endl; 	 
     f0.getline(buf, 1024);
     cout << buf << endl;
	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> temp0medium >> VXmedium >> VYmedium >> VZmedium;	 
     //f0 >> temp0 >> vc0[1] >> vc0[2] >> vc0[3];
     cout << temp0medium << "            " << VXmedium << "           " << VYmedium << "            " << VZmedium<<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;
	 	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> tau0 >> tauend >> dtau;
     cout << tau0 << "      " << tauend << "         " << dtau <<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;
	 	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> tauhydro0 >> nxhydro >> nyhydro >> netahydro >> dtauhydro;
     cout << tauhydro0 << "            " << nxhydro << "           " << nyhydro << "           " << netahydro << "             " << dtauhydro <<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;
	 
     f0.getline(buf, 1024);
     f0.getline(buf, 1024);
	 
     f0.getline(buf, 1024);
     cout << buf << endl;
     f0 >> readin_route;
     cout << readin_route <<endl;
     f0.getline(buf, 1024);
     cout << buf << endl;	 
	 
     }
	  
     f0.close();
	 
	 
//...alphas!	
     //alphas=alphas0(Kalphas,temp0);
	 
//...Debye Mass square
     qhat0=DebyeMass2(Kqhat0,alphas,temp0);	

//


/////.....................test


//    exit(1);


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	 

//...initial particle list input end





















     cout<<"Loading datafile......"<<endl;
	
//........................................................................................................................................................................	

//......scattering rate reading alphas!	
     int it,ie;
     int n=360;	 
	 
//...Hydroprofile
     tauhydro0=0.2;
     nxhydro=73;
     nyhydro=73;
     netahydro=33;
     dtauhydro=0.3;	 

	 
	 
/*	 
///////////////////////////////////////////////////////////////////////////////////////////	 
    //pythia.readString(Form("PhaseSpace:pTHatMin = %f",ptbinmin));	 
	 
	 //char readin_route[1024]="readindatafile";
	 //double ptbinmin=100.0;

	char routeline[1024];
	sprintf(tmpbuf,"PhaseSpace:pTHatMin = %s",readin_route);
    cout<<tmpbuf<<endl;	
	//cout<<Form("PhaseSpace:pTHatMin = %f",ptbinmin)<<endl;
	 
	 exit(1);
///////////////////////////////////////////////////////////////////////////////////////////
*/




	 char ratedata_route[1024];
     sprintf(ratedata_route,"%s/ratedata",readin_route);   
     ifstream f1(ratedata_route);
	 cout<<ratedata_route<<endl;
	 
     //ifstream f1("readindatafile/ratedata");
     if(!f1.is_open())
     {
      cout<<"Erro openning date file1!\n";
     }
     else
     {
      for(int i=1;i<=n;i++)
	  {
	  f1>>it>>ie;	 	  
	  //f1>>Rg[it][ie]>>Rg1[it][ie]>>Rg2[it][ie]>>Rg3[it][ie]>>Rq[it][ie]>>Rq3[it][ie]>>Rq4[it][ie]>>Rq5[it][ie]>>Rq6[it][ie]>>Rq7[it][ie]>>Rq8[it][ie];
	  f1>>qhatG[it][ie]>>Rg[it][ie]>>Rg1[it][ie]>>Rg2[it][ie]>>Rg3[it][ie]>>qhatLQ[it][ie]>>Rq[it][ie]>>Rq3[it][ie]>>Rq4[it][ie]>>Rq5[it][ie]>>Rq6[it][ie]>>Rq7[it][ie]>>Rq8[it][ie];	  
	  }
     }
     f1.close();	 

	 char Hydroprofile_Tau_route[1024];
     sprintf(Hydroprofile_Tau_route,"%s/Hydroprofile/tau.dat",readin_route);   
     ifstream f21(Hydroprofile_Tau_route);
	 cout<<Hydroprofile_Tau_route<<endl;	 
	 
     //ifstream f21("readindatafile/Hydroprofile/tau.dat");
     if(!f21.is_open())
     {
      cout<<"Erro openning date file2!\n";
     }
     else
     {
	  f21>>tauhydroend;	       
     }
     f21.close();

     ntauhydro=floor((tauhydroend-tauhydro0)/dtauhydro)+1;
	 tauhydrofile=tauhydro0+floor((tauhydroend-tauhydro0)/dtauhydro)*dtauhydro;
     
	 double tauh,xh,yh,etah,ed;



	 char Hydroprofile_bulk3D_route[1024];
     sprintf(Hydroprofile_bulk3D_route,"%s/Hydroprofile/bulk3D.dat",readin_route);   
     ifstream f2(Hydroprofile_bulk3D_route);
	 cout<<Hydroprofile_bulk3D_route<<endl;

	 
     //ifstream f2("readindatafile/Hydroprofile/bulk3D.dat");
     if(!f2.is_open())
     {
      cout<<"Erro openning date file2!\n";
     }
     else
     {
      for(int itau=1;itau<=ntauhydro;itau++)
      {
      for(int ix=1;ix<=nxhydro;ix++)
      {
      for(int iy=1;iy<=nyhydro;iy++)
      {
      for(int ieta=1;ieta<=netahydro;ieta++)
      {
	  
	  f2>>tauh>>ed>>temphydro[itau][ix][iy][ieta]>>fractionhydro[itau][ix][iy][ieta]>>VXhydro[itau][ix][iy][ieta]>>VYhydro[itau][ix][iy][ieta]>>VZhydro[itau][ix][iy][ieta];
	  
      }
      }
      }
      }

     }
     f2.close();







//////////////////////////////////////////////////////////////////////////////////////////////////...READIN NEW
	 
// duplicate for heavy quark


	 char ratedata_HQ_route[1024];
     sprintf(ratedata_HQ_route,"%s/ratedata-HQ",readin_route);   
     ifstream f11(ratedata_HQ_route);
	 cout<<ratedata_HQ_route<<endl;


     //ifstream f11("readindatafile/ratedata-HQ");
     if(!f11.is_open())
     {
      cout<<"Erro openning HQ data file!\n";
     }
     else
     {
      for(int i=1;i<=n;i++)
	  {
	  f11>>it>>ie;
	  f11>>RHQ[it][ie]>>RHQ11[it][ie]>>RHQ12[it][ie]>>qhatHQ[it][ie];	      
	  }
     }
     f11.close();

//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	 


//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	 
	 
	 
// read radiation table for heavy quark


	 char dNg_over_dt_cD6_route[1024];
     sprintf(dNg_over_dt_cD6_route,"%s/dNg_over_dt_cD6.dat",readin_route);   
     ifstream f12(dNg_over_dt_cD6_route);
	 cout<<dNg_over_dt_cD6_route<<endl;

	 char dNg_over_dt_qD6_route[1024];
     sprintf(dNg_over_dt_qD6_route,"%s/dNg_over_dt_qD6.dat",readin_route);   
     ifstream f13(dNg_over_dt_qD6_route);
	 cout<<dNg_over_dt_qD6_route<<endl;

	 char dNg_over_dt_gD6_route[1024];
     sprintf(dNg_over_dt_gD6_route,"%s/dNg_over_dt_gD6.dat",readin_route);   
     ifstream f14(dNg_over_dt_gD6_route);
	 cout<<dNg_over_dt_gD6_route<<endl;

     //ifstream f12("readindatafile/dNg_over_dt_cD6.dat");
     //ifstream f13("readindatafile/dNg_over_dt_qD6.dat");
     //ifstream f14("readindatafile/dNg_over_dt_gD6.dat");
     if(!f12.is_open()||!f13.is_open()||!f14.is_open())
     {
      cout<<"Erro openning HQ radiation table file!\n";
     }
     else
     {
     for(int k=1; k<=t_gn; k++) {
         char dummyChar[100];
         long double dummyD;
         f12 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         f13 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         f14 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         for(int i=1; i<=temp_gn; i++) {
            for(int j=1; j<=HQener_gn; j++) {
                f12 >> dummyD >> dummyD >> dummyD >> dNg_over_dt_c[k+1][i][j] >> max_dNgfnc_c[k+1][i][j];
                f13 >> dummyD >> dummyD >> dummyD >> dNg_over_dt_q[k+1][i][j] >> max_dNgfnc_q[k+1][i][j];
                f14 >> dummyD >> dummyD >> dummyD >> dNg_over_dt_g[k+1][i][j] >> max_dNgfnc_g[k+1][i][j];
            }
         }
     }
     }
     cout << dNg_over_dt_c[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_c[t_gn+1][temp_gn][HQener_gn] << endl;
     cout << dNg_over_dt_q[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_q[t_gn+1][temp_gn][HQener_gn] << endl;
     cout << dNg_over_dt_g[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_g[t_gn+1][temp_gn][HQener_gn] << endl;
     f12.close();
     f13.close();
     f14.close();
	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////
	 
//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	 
	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////	 

     for(int i=1; i<=temp_gn; i++) {
         for(int j=1; j<=HQener_gn; j++) {
             dNg_over_dt_c[1][i][j]=0.0;
             dNg_over_dt_q[1][i][j]=0.0;
             dNg_over_dt_g[1][i][j]=0.0;
             max_dNgfnc_c[1][i][j]=0.0;
             max_dNgfnc_q[1][i][j]=0.0;
             max_dNgfnc_g[1][i][j]=0.0;
         }
     }

//////////////////////////////////////////////////////////////////////////////////////////////////	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////...NEW


//////////////////////////////////////////////////////////////////////////////////////////////////	 
	 
// preparation for HQ 2->2
     min_p1=0.0;
     max_p1=200.0;
     bin_p1=(max_p1-min_p1)/N_p1;
     min_T=0.1;
     max_T=0.6;
     bin_T=(max_T-min_T)/N_T;
     min_e2=0.0;
     max_e2=15.0;
     bin_e2=(max_e2-min_e2)/N_e2;


	 char distB_route[1024];
     sprintf(distB_route,"%s/distB.dat",readin_route);   
     ifstream fileB(distB_route);
	 cout<<distB_route<<endl;
	 
     //ifstream fileB("readindatafile/distB.dat");
     if(!fileB.is_open()) {
        cout << "Erro openning data file distB.dat!" << endl;
     } else {
        for(int i=0;i<N_T;i++) {
           for(int j=0;j<N_p1;j++) {
              double dummy_T,dummy_p1;
              fileB>>dummy_T>>dummy_p1;
              if(fabs(min_T+(0.5+i)*bin_T-dummy_T)>1.0e-5 || fabs(min_p1+(0.5+j)*bin_p1-dummy_p1)>1.0e-5) {
                  cout << "Erro in reading data file distB.dat!" << endl;
                  exit (EXIT_FAILURE);
              }
              fileB>>distFncBM[i][j];
              for(int k=0;k<N_e2;k++) fileB>>distFncB[i][j][k];
              for(int k=0;k<N_e2;k++) fileB>>distMaxB[i][j][k];
           }
        }
     }
     fileB.close();
	 
//////////////////////////////////////////////////////////////////////////////////////////////////	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	 
	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////	 
 


	 char distF_route[1024];
     sprintf(distF_route,"%s/distF.dat",readin_route);   
     ifstream fileF(distF_route);
	 cout<<distF_route<<endl;


 
     //ifstream fileF("readindatafile/distF.dat");
     if(!fileF.is_open()) {
        cout << "Erro openning data file distF.dat!" << endl;
     } else {
        for(int i=0;i<N_T;i++) {
           for(int j=0;j<N_p1;j++) {
              double dummy_T,dummy_p1;
              fileF>>dummy_T>>dummy_p1;
              if(fabs(min_T+(0.5+i)*bin_T-dummy_T)>1.0e-5 || fabs(min_p1+(0.5+j)*bin_p1-dummy_p1)>1.0e-5) {
                  cout << "Erro in reading data file distF.dat!" << endl;
                  exit (EXIT_FAILURE);
              }
              fileF>>distFncFM[i][j];
              for(int k=0;k<N_e2;k++) fileF>>distFncF[i][j][k];
              for(int k=0;k<N_e2;k++) fileF>>distMaxF[i][j][k];
           }
        }
     }
     fileF.close();

//////////////////////////////////////////////////////////////////////////////////////////////////	 
	 
//////////////////////////////////////////////////////////////////////////////////////////////////...NEW	



     cout<<"datafile loaded"<<endl;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//...random seed setting
     srand((unsigned)time(NULL));
//     NUM1=-1*rand();
   //NUM1=-6;

     struct tm *local_start;
     time_t time_start;
     time_start=time(NULL);
     local_start=localtime(&time_start);

//...the program starts	 	 
     char buf1[80];
     strftime(buf1,80,"Current Time: %Y-%m-%d %H:%M:%S",local_start);
     cout << "the program starts at:" <<endl;
     cout << buf1 << endl;	 	 
	 
//...................................................................................................................................................Input end.	 


//output data block

//..datafile of all partons for recombination
//    ofstream positiveLBT("./positive.dat");
//    ofstream negativeLBT("./negative.dat");

//    ofstream numpositiveLBT("./numptpositive.dat");
//    ofstream numnegativeLBT("./numptnegative.dat");

}


void LBTclass::LBTclear(){	        
        for(int i=1;i<=50000;i++)		//clear the dimensions for each event
        {
            P[1][i]=0.0;
            P[2][i]=0.0;
            P[3][i]=0.0;
            P[0][i]=0.0;
		    //
            P0[1][i]=0.0;
            P0[2][i]=0.0;
            P0[3][i]=0.0;
            P0[0][i]=0.0;
            //
            PP[1][i]=0.0;
            PP[2][i]=0.0;
            PP[3][i]=0.0;
            PP[0][i]=0.0;
            //
            PP0[1][i]=0.0;
            PP0[2][i]=0.0;
            PP0[3][i]=0.0;
            PP0[0][i]=0.0;
            //
            V[1][i]=0.0;
            V[2][i]=0.0;
            V[3][i]=0.0;
            //
            V0[1][i]=0.0;
            V0[2][i]=0.0;
            V0[3][i]=0.0;
            //
            VV[1][i]=0.0;
            VV[2][i]=0.0;
            VV[3][i]=0.0;
            //
            VV0[1][i]=0.0;
            VV0[2][i]=0.0;
            VV0[3][i]=0.0;
            //
            CAT[i]=0;
            CAT0[i]=0;
            //
            radng[i]=0.0;
            tirad[i]=0.0;
            tiscatter[i]=0.0;
			
	
			//
            Vfrozen[0][i]=0.0; 
            Vfrozen[1][i]=0.0; 
            Vfrozen[2][i]=0.0; 
            Vfrozen[3][i]=0.0; 
            Vfrozen0[0][i]=0.0;   
            Vfrozen0[1][i]=0.0;   
            Vfrozen0[2][i]=0.0;   
            Vfrozen0[3][i]=0.0;   
            Tfrozen[i]=0.0;
            Tfrozen0[i]=0.0;
            vcfrozen[0][i]=0.0;
            vcfrozen[1][i]=0.0;
            vcfrozen[2][i]=0.0;
            vcfrozen[3][i]=0.0;
            vcfrozen0[0][i]=0.0;
            vcfrozen0[1][i]=0.0;
            vcfrozen0[2][i]=0.0;
            vcfrozen0[3][i]=0.0;

            Tint_lrf[i]=0.0; //for heavy quark			

            ////////////////////////////////////////////////////////////////////////////...2019
            NS[i]=0;                   //splitting rank

            NR0[i]=0;                  //scattering rank
            NS0[i]=0;                  //splitting rank

            Mother[i]=0;               //ID of initial parton of creation
////////////////////////////////////////////////////////////////////////////...2019





            NR[i]=0;                  //scattering rank 
			tiform[i]=0.0;


			
        }		
	

//...initial particle list input

//...single jet parton input

     if(switchsingle==1)
	 {

	 P[1][1]=px0;
	 P[2][1]=py0;
	 P[3][1]=pz0;
	 P[0][1]=en0;	

     KATT1[1]=Kjet;	 

	 nj=1;
	 
	 V[1][1]=Xinitial;
	 V[2][1]=Yinitial;	 
	 V[3][1]=Zinitial;	 	 
	 
	 }
//...initial particle list input end

	
//......single initial jet parton 
        nj=1;		

//......initial parton shower
        np=nj;

}


