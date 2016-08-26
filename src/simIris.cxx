#include <stdlib.h>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TStopwatch.h"

#include "header.h"
#include "detHits.h"

params prm, prm_inp;
YYHit yd;
CsIHit csi;
S3Hit sd1, sd2;

Double_t mA;	
Double_t ma;	
Double_t mB;	
Double_t mBR;	
Double_t mb;	
Double_t mc;
Double_t md;

PTrack hP, lP, decP;

Double_t Qgen, Qdet;
Double_t beamE, beamBeta, beamGamma, beamEcm;

void clearEvt()
{
	TargdE[0]=0.; TargdE[1]=0.;
	Qgen=0.; Qdet=0.;
	mBR=0.;
	lP.Clear();
	hP.Clear();
	decP.Clear();
	yd.Clear();
	csi.Clear();
	sd1.Clear();
	sd2.Clear();
	
	return;
}

int main(int argc, char *argv[])
{
	Bool_t load_params_from_file=kFALSE;
	char *endptr;
	Int_t nsim = 1E6;
	char *paramsname =NULL;
	char *dedxpath =NULL;
	char *outputname =NULL;
	
	if (argc > 1){
		for(int i=0; i<argc; i++){
			if(strncmp(argv[i],"--output=",9)==0){
				outputname = argv[i]+9;
			}
			else if(strncmp(argv[i],"--dedx_dir=",11)==0){
				dedxpath = argv[i]+11;
			}
			else if(strncmp(argv[i],"--params=",9)==0){
				load_params_from_file=kTRUE;
				paramsname = argv[i]+9;
			}
			else if(strncmp(argv[i],"--events=",9)==0){
	  			nsim = strtol(argv[i]+9,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--N=",4)==0){
	  			prm.N = strtol(argv[i]+4,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--A=",4)==0){
	  			prm.A = argv[i]+4;
			}
			else if(strncmp(argv[i],"--a=",4)==0){
	  			prm.a = argv[i]+4;
			}
			else if(strncmp(argv[i],"--b=",4)==0){
	  			prm.b = argv[i]+4;
			}
			else if(strncmp(argv[i],"--c=",4)==0){
	  			prm.c = argv[i]+4;
			}
			else if(strncmp(argv[i],"--d=",4)==0){
	  			prm.d = argv[i]+4;
			}
			else if(strncmp(argv[i],"--B=",4)==0){
	  			prm.B = argv[i]+4;
			}
			else if(strncmp(argv[i],"--E=",4)==0){
	  			prm.E  = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--R=",4)==0){
	  			prm.R = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--W=",4)==0){
	  			prm.W = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--Tt=",5)==0){
	  			prm.Tt = strtod(argv[i]+5,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--Bs=",5)==0){
	  			prm.Bs = strtod(argv[i]+5,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--DYY=",6)==0){
	  			prm.DYY = strtod(argv[i]+6,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--DS3=",6)==0){
	  			prm.DS3 = strtod(argv[i]+6,&endptr);//converting string to number
			}
		}
	}
	TStopwatch timer;
	timer.Start();
	
	if(load_params_from_file==kTRUE){
		printf("Loading parameters from file %s.\n",paramsname);
		prm.Load(paramsname);
	}
	prm.Print();

	TRandom3 *rndm = new TRandom3(0);
	Int_t Evnt=0; 
	Double_t chck, chck2;
	Double_t wght, wght2;

	Double_t ICdE;	
	YYHit *ipyd = &yd; 
	CsIHit *ipcsi = &csi; 
	S3Hit *ipsd1 = &sd1; 
	S3Hit *ipsd2 = &sd2; 
	PTrack *iplP = &lP;
	PTrack *iphP = &hP;
	PTrack *ipdecP = &decP;
	nucleus A, a, B, b, c, d, decB,decc,decd;
	
	A.getInfo(prm.A);
	a.getInfo(prm.a);
	B.getInfo(prm.B);
	b.getInfo(prm.b);
	c.getInfo(prm.c);
	d.getInfo(prm.d);

	mA = A.mass/1000.;	
	ma = a.mass/1000.;	
	mB = B.mass/1000.+prm.R/1000.;	
	mBR = mB;	
	mb = b.mass/1000.;	
	mc = c.mass/1000.;
	md = d.mass/1000.;

// Check for sequential decays ****************************************
	Bool_t seqdec=kFALSE;
	Double_t S_low=50.;
	Int_t seqdecN=0;
	Int_t pick=0;
	Double_t mBdec=0., mcdec=0., mddec=0.;
	Double_t masses2[3];

	if(B.Sn!=0.&&B.Sn<S_low){
		S_low=B.Sn;
		pick=1;
	}
	if(B.Sp!=0.&&B.Sp<S_low){
		S_low=B.Sp;
		pick=2;
	}
	if(B.S2n!=0.&&B.S2n<S_low){
		S_low=B.S2n;
		pick=3;
	}
	if(B.S2p!=0.&&B.S2p<S_low){
		S_low=B.S2p;
		pick=4;
	}
	printf("\nResonance Energy: %.2lf\tlowest threshold: %.2lf\n",prm.R,S_low);

	if(S_low<prm.R){
		switch(pick){
			case 1 : 
				seqdec = kTRUE;
				printf("\nSequential 1n- decay!\n\n"); 
				seqdecN =2;
				decB.getInfo(B.N-1,B.Z);
				mBdec=decB.mass/1000.;
				decc.getInfo("n");
				mcdec=decc.mass/1000.;
				break;
			case 2 : 
				seqdec = kTRUE;
				printf("\nSequential 1p- decay!\n\n"); 
				seqdecN =2;
				decB.getInfo(B.N,B.Z-1);
				mBdec=decB.mass/1000.;
				decc.getInfo("p");
				mcdec=decc.mass/1000.;
				break;
			case 3 : 
				seqdec = kTRUE;
				printf("\nSequential 2n- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(B.N-2,B.Z);
				mBdec=decB.mass/1000.;
				decc.getInfo("n");
				mcdec=decc.mass/1000.;
				decd.getInfo("n");
				mddec=decd.mass/1000.;
				break;
			case 4 : 
				seqdec = kTRUE;
				printf("\nSequential 2p- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(B.N,B.Z-2);
				mBdec=decB.mass/1000.;
				decc.getInfo("p");
				mcdec=decc.mass/1000.;
				decd.getInfo("p");
				mddec=decd.mass/1000.;
				break;
			default : 
				seqdec = kFALSE;
				break;
		}
		masses2[0] = mBdec;
		masses2[1] = mcdec;
		masses2[2] = mddec;
	}
// *******************************************************************

	// Set up output file and tree
	TFile *f = new TFile(outputname,"RECREATE");
	TTree *iris = new TTree("iris","iris simulation");
	
	iris->Branch("Evnt",&Evnt,"Evnt/I"); 
	iris->Branch("beamE",&beamE,"beamE/D"); 
	iris->Branch("beamBeta",&beamBeta,"beamBeta/D"); 
	iris->Branch("beamGamma",&beamGamma,"beamGamma/D"); 
	iris->Branch("beamEcm",&beamEcm,"beamEcm/D"); 
	iris->Branch("lP",&iplP,32000,99); 
	iris->Branch("hP",&iphP,32000,99); 
	iris->Branch("decP",&ipdecP,32000,99); 
	// iris->Branch("E",Etree,"Etree[2]/D"); 
	// iris->Branch("E2",En,"En[2]/D"); 
	// iris->Branch("Ecm",Ecm,"Ecm[2]/D"); 
	// iris->Branch("Theta",Td,"Td[2]/D"); 
	// iris->Branch("Theta2",Td2,"Td2[2]/D"); 
	// iris->Branch("Thetacm",Tcm,"Tcm[2]/D"); 
	// iris->Branch("Phi",Pd,"Pd[2]/D"); 
	// iris->Branch("Phi2",Pd2,"Pd2[2]/D"); 
	iris->Branch("TargdE",TargdE,"TargdE[2]/D"); 
	iris->Branch("wght",&wght,"wght/D"); 
	iris->Branch("Qgen",&Qgen,"Qgen/D"); 
	iris->Branch("Qdet",&Qdet,"Qdet/D"); 
	iris->Branch("mBR",&mBR,"mBR/D"); 
	iris->Branch("ICdE",&ICdE,"ICdE/D"); 
	iris->Branch("yd",&ipyd,32000,99); 
	iris->Branch("csi",&ipcsi,32000,99); 
	iris->Branch("sd1",&ipsd1,32000,99); 
	iris->Branch("sd2",&ipsd2,32000,99); 


	std::string dedxstr = dedxpath;
	if(!seqdec) loadAllELoss(dedxstr,A,B,b);
	else loadAllELoss(dedxstr,A,decB,b);

	Double_t targetTh=prm.Tt*0.0867*0.1; //mu*g/cm^3*0.1
	Double_t BeamSpot=prm.Bs/2.355; // FWHM->sigma 
	const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr
	const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	//const Double_t AgFoil=5.44*10.473*0.1; //mu*g/cm^3*0.1
	const Double_t AgFoil=4.355*10.473*0.1; //mu*g/cm^3*0.1

	Bool_t LEHit, HEHit;

	// Calculate energy loss up to center of the target
	Double_t EA = prm.E;	
   	EA -= simEloss(A,0.5,EA,ICWindow1,eASi3N4, dedxASi3N4);
   	ICdE = simEloss(A,0.586,EA,ICLength,eAC4H10, dedxAC4H10);
   	EA -= ICdE;
   	EA -= simEloss(A,0.5,EA,ICWindow2,eASi3N4, dedxASi3N4);
   	EA -= simEloss(A,47./108.,EA,AgFoil,eAAg, dedxAAg);
   	EA -= simEloss(A,1.,EA,targetTh/2.,eAH, dedxAH);
 
	EA = EA/1000.; // convert to GeV for TGenPhaseSpace
	Double_t PA = sqrt(EA*EA+2*EA*mA);

	TLorentzVector target(0.0, 0.0, 0.0, ma);
	TLorentzVector beam(0.0, 0.0, PA, mA+EA);
	TLorentzVector Sys = beam + target;
	TVector3 boostvect = Sys.BoostVector();

	beamE = EA*1000.;
	beamBeta = Sys.Beta();
	beamGamma = Sys.Gamma();
	beamEcm = EA*ma*1000./(mA+ma);
	printf("\n\nEnergy at center of target: %.2lf MeV\n", beamE);
	printf("\nBeta at center of target: %.3lf \n", beamBeta);
	printf("\nGamma at center of target: %.3lf \n", beamGamma);
	printf("\nCM Energy at center of target: %.2lf MeV\n\n", beamEcm);

	Double_t masses[4] = { mb, mB, mc, md};

	TGenPhaseSpace PS0, PS1;
	Bool_t allowed = PS0.SetDecay(Sys, prm.N, masses);

	if(!allowed){
		printf("Impossible decay!\n");
		printf("Exiting...\n");
		exit(0);
	}

	Double_t wght_max=PS0.GetWtMax();
	printf("%lf\n",wght_max);
	Double_t width = prm.W/2.355/1000.;
	printf("%lf\t%lf\n",mB,width);

	while(Evnt<nsim) 
	{
		wght = 0.;
		clearEvt();
		mBR = rndm->Gaus(mB,width);
		masses[1] =mBR;
		PS0.SetDecay(Sys, prm.N, masses);

		do{	
			wght = PS0.Generate();
			chck = rndm->Uniform(0,1);

		}while(wght<chck);

		TLorentzVector *LVb  = PS0.GetDecay(0);
		TLorentzVector *LVB  = PS0.GetDecay(1);
		TLorentzVector *LVBdec;

		TLorentzVector Frag  = Sys-*LVb;

		Qgen= (mA+ma-mb-Frag.M())*1000.;	
	
		lP.T=LVb->Theta();	
		hP.T=LVB->Theta();
		lP.E=(LVb->E()-mb)*1000.; 	
		hP.E=(LVB->E()-mB)*1000.;
		lP.P=LVb->Phi();	
		hP.P=LVB->Phi();	
		
		// Convert angles to degrees for root file
		lP.Tdeg=RadToDeg()*lP.T;
		hP.Tdeg=RadToDeg()*hP.T;
		lP.Pdeg=RadToDeg()*lP.P;
		hP.Pdeg=RadToDeg()*hP.P;
		// Etree[0] = lP.E;
		// Etree[1] = hP.E;

		if(seqdec)
		{
			PS1.SetDecay(*LVB, seqdecN, masses2);
			do{
				wght2 = PS1.Generate();			
				chck2 = rndm->Uniform(0,1);
				LVBdec  = PS1.GetDecay(0);
			}while(wght2<chck2);
			decP.T=LVBdec->Theta();
			decP.E=(LVBdec->E()-mBdec)*1000.;
			decP.P=LVBdec->Phi();	
			decP.Tdeg=RadToDeg()*hP.T;
			decP.Pdeg=RadToDeg()*hP.P;
		}
	
		// Position on target	
		Double_t  targetPosX = BeamSpot*rndm->Gaus();
		Double_t  targetPosY = BeamSpot*rndm->Gaus();
		TVector3 targetPos(targetPosX,targetPosY,0);
		
		LEHit = detHitsL(lP, b, targetTh, targetPos);
		
		if(!seqdec){
			HEHit = detHitsH(hP, B, targetTh, targetPos);	
		}
		else{ 
			HEHit = detHitsH(decP, decB, targetTh, targetPos);
		}
		
		if(LEHit && yd.dE[0]>0.){
			Double_t Pb = LVb->P();
			Double_t Eb = LVb->E()-mb;	
 			Qdet = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EA)*(mb+Eb)+2.*PA*Pb*cos(yd.fThetaCalc[0]*DegToRad())+2.*(EA+mA+ma-Eb-mb)*ma);
			Qdet =Qdet*1000.;
		}

		lP.Ecm = (LVb->E()-mb)*ma*1000./(mA+ma);
		hP.Ecm = (LVB->E()-mB)*ma*1000./(mA+ma);
		LVb->Boost(-boostvect);
		LVB->Boost(-boostvect);
		lP.Tcm = RadToDeg()*(Pi()-LVb->Theta());
		hP.Tcm = RadToDeg()*LVB->Theta();

		printf("%.5d Events processed..\r",Evnt);
		Evnt++;
		iris->Fill();
	}

	iris->AutoSave();
	f->Close();
	Double_t time=timer.RealTime();
	printf("\nDone. %lf s\n",time);
	
	return 0;
}
