!****************************************************
!     Running recombination module                  !
!****************************************************
      program main
      IMPLICIT DOUBLE PRECISION(a-H, O-Z)
      common/const/pi,hbarc
      common/parm/sigma,sigmapro,sigmaK,maxf      

      DIMENSION P(5000,4),XV(5000,4)
      DIMENSION Kp(5000,4),idsh(5000),Kpx(5000)
      DIMENSION PH(5000,4),Kst(5000),Qscale(5000)
      DIMENSION KfH(5000), Pth(5000,4),XH(5000,4)
      DIMENSION Ppx(5000,4), Xpx(5000,4)
      dimension Xth(5000,4),Kth(5000),qfscale(5000)
      dimension xnbin(100,2), xnbinK(100,2)
      dimension xnbinpro(100,2)
      DIMENSION NCTKr(5000), CTPr(5000,4),CTXr(5000,4)

C.... input files .........................
!      open (unit=10, file='shower.dat', status='unknown')
!      open (unit=9, file='numsh.dat', status='unknown')
!      open (unit=15, file='thermal.dat', status='unknown')

!      open (unit=30, file='parton/LBT_result', status='unknown')
!      open (unit=45, file='parton/oscar.dat', status='unknown')




c..........................................LBT-Pythia	 
	 
      open (unit=71, file='../../outputdatafile/positive.dat',
     & status='unknown')
      open (unit=72, file='../../outputdatafile/numptpositive.dat',
     & status='unknown')	 	 

      open (unit=81, file='../../outputdatafile/thermal.dat', 
     & status='unknown')	 
	 
c..........................................LBT-Pythia















C.... output files ........................
!      open (unit=31, file='ana_pion.dat', status='unknown')
!      open (unit=32, file='ana_kaon.dat', status='unknown')
!      open (unit=33, file='ana_proton.dat', status='unknown')
      open (unit=51, file='thermal_jet.dat', status='unknown')!output file !noted by wenbin 2018.11.30
      open (unit=52, file='remnant_jet_parton.dat', status='unknown')!output file !noted by wenbin 2018.11.30
	  
      open (unit=58, file='num_remnant.dat', status='unknown')!output file !noted by wenbin 2018.11.30	  
      open (unit=59, file='remnant.dat', status='unknown')!output file !noted by wenbin 2018.11.30	  
	  
      open (unit=53, file='coalesced_thermal_partons.dat'
     .             ,status='unknown')!output file !noted by wenbin 2018.11.30


c..........................................LBT-Pythia
      open (unit=91, file='./output/hadron.dat', status='unknown')
      open (unit=92, file='./output/num.dat', status='unknown') 
c..........................................LBT-Pythia



!      open (unit=217, file='numofreco.dat', status='unknown')
      

      open (unit=88, file='input.txt', status='unknown')
!.....Read input parameters ....... 
      READ(88,*) NEV
      READ(88,*) sigma
      READ(88,*) sigmaK
      READ(88,*) sigmapro
      READ(88,*) maxf

!      NEV = 50000
!.........................      
      step = 30./40.
      pi = 2.*asin(1.)
      
      ntherm = 0
      xnevt = float(NEV)
      nremq = 0
      numq = 0
      kpi = 0
      kK = 0
      kpro = 0
      massud=0.25
      masss=0.43
C.... initialization of binning ..........
      do ini=1, 40
         do jni=1, 2
            xnbin(ini,jni) = 0.
            xnbinK(ini,jni) = 0.
            xnbinpro(ini,jni) = 0.
         enddo
      enddo

!      do itherm=1, ntherm
!         READ(15,*) Kth(itherm),Pth(itherm,1),Pth(itherm,2),
!     .        Pth(itherm,3),pth(itherm,4),Xth(itherm,1),
!     .        Xth(itherm,2),Xth(itherm,3),Xth(itherm,4)       !KTh thermal parton ID  
!      enddo
      
      DO IEV=1, NEV             !...event loop ..............
	  
	  
	  
	  
	  
	  
	  
!...test energy
        pxtotal=0.d0
        pytotal=0.d0
        pztotal=0.d0
        Etotal=0.d0
        pTtotal=0.d0

        pxout=0.d0
        pyout=0.d0
        pzout=0.d0
        Eout=0.d0	  
        pTout=0.d0
		
        pxRemnant=0.d0
        pyRemnant=0.d0
        pzRemnant=0.d0
        ERemnant=0.d0	  
        pTRemnant=0.d0	  
	  
	  
	  
	  
	  
	  
!..... read thermal parton from hydro !added by wenbin 2018.11.26...
!         READ(45,*) mid,ntherm,mid,mid
!      do itherm=1, ntherm
!         READ(45,*)mid,Kth(itherm),Pth(itherm,1),
!     .        Pth(itherm,2),Pth(itherm,3),Pth(itherm,4),amid,
!     .  Xth(itherm,1),Xth(itherm,2),Xth(itherm,3),Xth(itherm,4)       !KTh thermal parton ID  
!      enddo
!..... read LBT parton from LBT result !added by wenbin 2018.11.26

!...test energy
        !pxtotal=pxtotal+Pth(itherm,1)
        !pytotal=pytotal+Pth(itherm,2)
        !pztotal=pztotal+Pth(itherm,3)
        !Etotal=Etotal+Pth(itherm,4)


c..........................................LBT-Pythia
         write(*,*) IEV
         READ(72,*) kev, numparton
         write(*,*) IEV, kev, numparton	
c..........................................LBT-Pythia



         ish = 0

!.....number of partons .........................
         !READ(30,*)kev,  numparton !kev the event id, numparton: number of partons wenbin 2018.11.23
         nchth = 0
!.....  Shower partons ..........................
         DO ish0=1, numparton


c..........................................LBT-Pythia
            READ(71,*) kev0,kID,Pxxx,Pyyy,
     .           Pzzz,amidE,XVxxx,XVyyy,XVzzz,
     .           XVttt,iCAT
c..........................................LBT-Pythia

                 !write(*,*) "----------------------1"

            etap=1.0/2.0*DLOG((amidE+Pzzz)/(amidE-Pzzz))

                 !write(*,*) "----------------------2"

            if(abs(etap).le.4.8)then
	
            ish=ish+1
	
            kev=kev0
            idsh(ish)=kID
            P(ish,1)=Pxxx
            P(ish,2)=Pyyy
            P(ish,3)=Pzzz
	        amid=amidE
            XV(ish,1)=XVxxx
            XV(ish,2)=XVyyy
            XV(ish,3)=XVzzz
            XV(ish,4)=XVttt		

c..........................................LBT-Pythia
!            READ(71,*) kev,idsh(ish),P(ish,1),P(ish,2),
!     .           P(ish,3),amid,XV(ish,1),XV(ish,2),XV(ish,3),
!     .           XV(ish,4)
c..........................................LBT-Pythia



!            write(*,*) kev,idsh(ish),P(ish,1),P(ish,2),
!     .           P(ish,3),amid,XV(ish,1),XV(ish,2),XV(ish,3),
!     .           XV(ish,4)





!            READ(30,*) idsh(ish),P(ish,1),P(ish,2),
!     .           P(ish,3),amid,mid,XV(ish,1),XV(ish,2),XV(ish,3),
!     .           XV(ish,4),Qscale(ish)!,amid,amid,amid,amid,amid,amid

          !write(*,*)idsh(ish),P(ish,1),P(ish,2)
         if((abs(idsh(ish)).eq.1).or.(abs(idsh(ish)).eq.2))then
              P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .         + P(ish,3)*P(ish,3)+0.0*0.0)
         !endif
         else
              if(abs(idsh(ish)).eq.3)then
                P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .           + P(ish,3)*P(ish,3)+0.0*0.0)
              else
                P(ish,4)=sqrt(P(ish,1)*P(ish,1)+P(ish,2)*P(ish,2)
     .           + P(ish,3)*P(ish,3)+0.0*0.0)
              endif
         endif	 
		 
!...test energy
        pxtotal=pxtotal+P(ish,1)
        pytotal=pytotal+P(ish,2)
        pztotal=pztotal+P(ish,3)
        Etotal=Etotal+P(ish,4)		 
		 
		 
		 
		 
		 
        !write(*,*)ntherm,numparton
!!.....number of partons .........................
!         READ(9,*) kev, numparton !kev the event id, numparton: number of partons wenbin 2018.11.23
!         nchth = 0
!!.....  Shower partons ..........................
!         DO ish=1, numparton
!            READ(10,*) kev,idsh(ish),Kp(ish,3),P(ish,1),P(ish,2),
!     .           P(ish,3),P(ish,4),XV(ish,1),XV(ish,2),XV(ish,3),
!     .           XV(ish,4)
!idsh:id, 

!..... counting number of quarks ..................
            if(idsh(ish).eq.21) then
               numq = numq + 2
            else
               numq = numq + 1
            endif

         endif
		 
         numpt=ish

         enddo
		 

                 !write(*,*) "----------------------3"		 
		 
         numparton=numpt
		 
		 
        !write(*,*)ntherm,numparton

          !write(*,*)"zzzzzzzz"
!         call hadronization(ntherm,Kth,Xth,Pth,numparton,idsh,XV,P,nH,
!     .        KfH,XH,PH,Kst)
         call hadronization(ntherm,Kth,Xth,Pth,numparton,idsh,XV,P,nH,
     .        KfH,XH,PH,Kst,Npt,Kpx,Ppx,Xpx, !added by wenbin for outputing remnant jet partons
     .        NCTnr, NCTKr,CTPr,CTXr,Qscale,qfscale)!CTnr, NCTKr,CTPr,CTXr the information of the coalesced hadrons Wenbin
!********* write the hadrons ********
           !write(*,*)"NCTnr,CTPr=1 ",NCTnr,CTPr(NCTnr,1) 
            write(51,17)IEV,nH
			
c..........................................LBT-Pythia			
            write(92,17)IEV,nH
c..........................................LBT-Pythia			
			
C.... Binning the produced hadrons ........................
         do ihad=1, nH
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
            write(51,16) IEV,KfH(ihad),PH(ihad,1),PH(ihad,2),
     .           PH(ihad,3),PH(ihad,4),XH(ihad,1),XH(ihad,2),XH(ihad,3),
     .           XH(ihad,4),Kst(ihad)!Kst is origin of the hadrons(th-th:0, sh-th:1, sh-sh:2,frag:3)


c..........................................LBT-Pythia
            write(91,26) IEV,Kst(ihad),KfH(ihad),PH(ihad,1),PH(ihad,2),
     .           PH(ihad,3),PH(ihad,4)
c..........................................LBT-Pythia

!...test energy
        pxout=pxout+PH(ihad,1)
        pyout=pyout+PH(ihad,2)
        pzout=pzout+PH(ihad,3)
        Eout=Eout+PH(ihad,4)








         enddo

!********* write the Remnant partons ******** 
            write(52,17)IEV,Npt
            write(58,17)IEV,Npt			
C.... Binning the writing ........................
         do ihad=1, Npt
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
            if (qfscale(ihad).gt.20.0)qfscale(ihad)=20.0
            write(52,18) IEV,Kpx(ihad),Ppx(ihad,1),Ppx(ihad,2),
     .           Ppx(ihad,3),Ppx(ihad,4),Xpx(ihad,1),Xpx(ihad,2),
     .           Xpx(ihad,3),Xpx(ihad,4),qfscale(ihad)

            write(59,18) IEV,Kpx(ihad),Ppx(ihad,1),Ppx(ihad,2),
     .           Ppx(ihad,3),Ppx(ihad,4)
	 
!                Xpx(ihad,3),Xpx(ihad,4),qfscale(ihad)	 
            !write(*,*) Ppx(ihad,3),Ppx(ihad,4)




!...test energy
        pxout=pxout+Ppx(ihad,1)
        pyout=pyout+Ppx(ihad,2)
        pzout=pzout+Ppx(ihad,3)
        Eout=Eout+Ppx(ihad,4)

        pxRemnant=pxRemnant+Ppx(ihad,1)
        pyRemnant=pyRemnant+Ppx(ihad,2)
        pzRemnant=pzRemnant+Ppx(ihad,3)
        ERemnant=ERemnant+Ppx(ihad,4)


         enddo

!********* write the coaleced thermal  partons ********
            write(53,*)"# ",IEV,NCTnr
C.... Binning the produced hadrons ........................
        if(NCTnr.eq.0)then
             write(53,*)1,0,0,0,0,0,0,0,0,0
             write(53,*)1,0,0,0,0,0,0,0,0,0

        endif   

         do ihad=1, NCTnr
!            pT = sqrt(PH(ihad,1)**2.+PH(ihad,2)**2.)
         !....... transfer the t,z to tau, etas...wenbin
           ttau=sqrt(CTXr(ihad,4)*CTXr(ihad,4)-CTXr(ihad,3)*CTXr(ihad,3)
     .              )
           eetas=0.5*log((CTXr(ihad,4)+CTXr(ihad,3))/(CTXr(ihad,4)
     .                   -CTXr(ihad,3)))
            write(53,18) IEV,NCTKr(ihad),CTPr(ihad,1),CTPr(ihad,2),
     .           CTPr(ihad,3),CTPr(ihad,4),CTXr(ihad,1),CTXr(ihad,2),
     .           eetas,ttau
         enddo





!...test energy

        pTtotal=sqrt(pxtotal**2+pytotal**2)
        pTout=sqrt(pxout**2+pyout**2)		
        pTRemnant=sqrt(pxRemnant**2+pyRemnant**2)
		
        write(*,*) pxtotal,pxRemnant,pxout,(pxout-pxtotal)/pxtotal
        write(*,*) pytotal,pyRemnant,pyout,(pyout-pytotal)/pytotal
        write(*,*) pztotal,pzRemnant,pzout,(pzout-pztotal)/pztotal
        write(*,*) Etotal,ERemnant,Eout,(Eout-Etotal)/Etotal
        write(*,*) pTtotal,pTRemnant,pTout,(pTout-pTtotal)/pTtotal 






       ENDDO

 15   format(1(f6.2),4(e21.8))
 16   format(1(I9),1(I7),8(f21.8),1(I9))
 17   format(1(I9),1(I9))
 18   format(1(I9),1(I7),9(f21.8))
 
c..........................................LBT-Pythia
 25   format(1(f6.2),4(e21.8))
 26   format(1(I9),1(I3),1(I5),8(f21.8))
c..........................................LBT-Pythia

       end

