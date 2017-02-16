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
S3Hit sd1(0,60.), sd2(1,1000.);
PTrack hP, lP;
PTrack decHP, declP1, declP2;
IDet det;

Double_t mA;	
Double_t ma;	
Double_t mB;	
Double_t mBR;	
Double_t mb;	
Double_t mc;
Double_t md;

Double_t Qgen, Qdet, Qalt;
Double_t EB_det, PB_det;
Double_t beamE, beamBeta, beamGamma, beamEcm;
Double_t tdE;
TVector3 reacPos;
Double_t SSBdE;	

void clearEvt()
{
	//TargdE[0]=0.; TargdE[1]=0.;
	Qgen=sqrt(-1.); Qdet=sqrt(-1.); Qalt=sqrt(-1.);
	EB_det=sqrt(-1.); PB_det=sqrt(-1.);
	mBR=0.;
	lP.Clear();
	hP.Clear();
	decHP.Clear();
	SSBdE=0.;
	yd.Clear();
	csi.Clear();
	sd1.Clear();
	sd2.Clear();
	det.Clear();
	
	return;
}

int main(int argc, char *argv[])
{
	Bool_t load_params_from_file=kFALSE;
	Bool_t have_dwba_xsec=kFALSE;
	char *endptr;
	Int_t nsim = 1E6;
	char *paramsname =NULL;
	char *dedxpath =NULL;
	char *outputname =NULL;
	char *dwbaname =NULL;
	Bool_t isAgReac = kFALSE;

	std::string binpath(argv[0]);
	printf("%s\n",binpath.data());

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
			else if(strncmp(argv[i],"--dwba=",7)==0){
				have_dwba_xsec=kTRUE;
				dwbaname = argv[i]+7;
			}
			else if(strncmp(argv[i],"--events=",9)==0){
	  			nsim = strtol(argv[i]+9,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--N=",4)==0){
	  			prm_inp.N = strtol(argv[i]+4,&endptr,10);//converting string to number
			}
			else if(strncmp(argv[i],"--A=",4)==0){
	  			prm_inp.A = argv[i]+4;
			}
			else if(strncmp(argv[i],"--a=",4)==0){
	  			prm_inp.a = argv[i]+4;
			}
			else if(strncmp(argv[i],"--b=",4)==0){
	  			prm_inp.b = argv[i]+4;
			}
			else if(strncmp(argv[i],"--c=",4)==0){
	  			prm_inp.c = argv[i]+4;
			}
			else if(strncmp(argv[i],"--d=",4)==0){
	  			prm_inp.d = argv[i]+4;
			}
			else if(strncmp(argv[i],"--B=",4)==0){
	  			prm_inp.B = argv[i]+4;
			}
			else if(strncmp(argv[i],"--E=",4)==0){
	  			prm_inp.E  = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--R=",4)==0){
	  			prm_inp.R = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--W=",4)==0){
	  			prm_inp.W = strtod(argv[i]+4,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--Tt=",5)==0){
	  			prm_inp.Tt = strtod(argv[i]+5,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--Bs=",5)==0){
	  			prm_inp.Bs = strtod(argv[i]+5,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--DYY=",6)==0){
	  			prm_inp.DYY = strtod(argv[i]+6,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--DS3=",6)==0){
	  			prm_inp.DS3 = strtod(argv[i]+6,&endptr);//converting string to number
			}
			else if(strncmp(argv[i],"--AgReac",9)==0){
				isAgReac=kTRUE;
			}
		}
	}
	TStopwatch timer;
	timer.Start();
	
	if(load_params_from_file==kTRUE){
		printf("Loading parameters from file %s.\n",paramsname);
		prm.Load(paramsname);

		if(!prm_inp.A.empty()) prm.A=prm_inp.A; 
		if(!prm_inp.a.empty()) prm.a=prm_inp.a;
		if(!prm_inp.B.empty()) prm.B=prm_inp.B;
		if(!prm_inp.b.empty()) prm.b=prm_inp.b;
		if(!prm_inp.c.empty()) prm.c=prm_inp.c;
		if(!prm_inp.d.empty()) prm.d=prm_inp.d;
		if(prm_inp.E>0.) prm.E=prm_inp.E;
		if(prm_inp.R>0.) prm.R=prm_inp.R;
		if(prm_inp.W>0.) prm.W=prm_inp.W;
		if(prm_inp.Tt>0.) prm.Tt=prm_inp.Tt; 
		if(prm_inp.Bs>0.) prm.Bs=prm_inp.Bs; 
		if(prm_inp.DYY>0.) prm.DYY=prm_inp.DYY; 
		if(prm_inp.DS3>0.) prm.DS3=prm_inp.DS3; 
		if(prm_inp.N>0) prm.N=prm_inp.N;
	}
	else{
		prm.A=prm_inp.A; 
		prm.a=prm_inp.a;
		prm.B=prm_inp.B;
		prm.b=prm_inp.b;
		prm.c=prm_inp.c;
		prm.d=prm_inp.d;
		prm.E=prm_inp.E;
		prm.R=prm_inp.R;
		prm.W=prm_inp.W;
		prm.Tt=prm_inp.Tt; 
		prm.Bs=prm_inp.Bs; 
		prm.DYY=prm_inp.DYY; 
		prm.DS3=prm_inp.DS3; 
		prm.N=prm_inp.N;
	}
	prm_inp.Print();
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
	PTrack *ipdecHP = &decHP;
	PTrack *ipdeclP1 = &declP1;
	PTrack *ipdeclP2 = &declP2;
	IDet *ipdet = &det;
	nucleus A, a, B, b, c, d, decB,decc,decd;
	Double_t reacX, reacY, reacZ;

	A.getInfo(binpath, prm.A);
	a.getInfo(binpath, prm.a);
	B.getInfo(binpath, prm.B);
	b.getInfo(binpath, prm.b);
	c.getInfo(binpath, prm.c);
	d.getInfo(binpath, prm.d);

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
				decB.getInfo(binpath, B.N-1,B.Z);
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "n");
				mcdec=decc.mass/1000.;
				break;
			case 2 : 
				seqdec = kTRUE;
				printf("\nSequential 1p- decay!\n\n"); 
				seqdecN =2;
				decB.getInfo(binpath, B.N,B.Z-1);
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "p");
				mcdec=decc.mass/1000.;
				break;
			case 3 : 
				seqdec = kTRUE;
				printf("\nSequential 2n- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(binpath, B.N-2,B.Z);
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "n");
				mcdec=decc.mass/1000.;
				decd.getInfo(binpath, "n");
				mddec=decd.mass/1000.;
				break;
			case 4 : 
				seqdec = kTRUE;
				printf("\nSequential 2p- decay!\n\n"); 
				seqdecN =3;
				decB.getInfo(binpath, B.N,B.Z-2);
				mBdec=decB.mass/1000.;
				decc.getInfo(binpath, "p");
				mcdec=decc.mass/1000.;
				decd.getInfo(binpath, "p");
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
	lP.nuc = b;
	hP.nuc = B;

	// Set up output file and tree
	TFile *f = new TFile(outputname,"RECREATE");
	TTree *Iris = new TTree("Iris","Iris simulation");
	
	Iris->Branch("Evnt",&Evnt,"Evnt/I"); 
	Iris->Branch("beamE",&beamE,"beamE/D"); 
	Iris->Branch("beamBeta",&beamBeta,"beamBeta/D"); 
	Iris->Branch("beamGamma",&beamGamma,"beamGamma/D"); 
	Iris->Branch("beamEcm",&beamEcm,"beamEcm/D"); 
	Iris->Branch("reacPos","TVector3",&reacPos); 
	Iris->Branch("lP.",&iplP,32000,99); 
	Iris->Branch("hP.",&iphP,32000,99); 
	Iris->Branch("decHP.",&ipdecHP,32000,99); 
	Iris->Branch("declP1.",&ipdeclP1,32000,99); 
	Iris->Branch("declP2.",&ipdeclP2,32000,99); 
	Iris->Branch("wght",&wght,"wght/D"); 
	Iris->Branch("Qgen",&Qgen,"Qgen/D"); 
	Iris->Branch("Qdet",&Qdet,"Qdet/D"); 
	//Iris->Branch("Qalt",&Qalt,"Qalt/D"); 
	//Iris->Branch("EB_det",&EB_det,"EB_det/D"); 
	//Iris->Branch("PB_det",&PB_det,"PB_det/D"); 
	//Iris->Branch("mBR",&mBR,"mBR/D"); 
	Iris->Branch("ICdE",&ICdE,"ICdE/D"); 
	Iris->Branch("SSBdE",&SSBdE,"SSBdE/D"); 
	Iris->Branch("yd.",&ipyd,32000,99); 
	Iris->Branch("csi.",&ipcsi,32000,99); 
	Iris->Branch("sd1.",&ipsd1,32000,99); 
	Iris->Branch("sd2.",&ipsd2,32000,99); 
	Iris->Branch("det",&ipdet,32000,99); 

	std::string dedxstr = dedxpath;
	A.EL.loadIncomingELoss(dedxstr,A.name.data(),A.mass);
	b.EL.loadOutgoingELoss(dedxstr,b.name.data(),b.mass);
	if(!seqdec) B.EL.loadOutgoingELoss(dedxstr,B.name.data(),B.mass);
	else decB.EL.loadOutgoingELoss(dedxstr,decB.name.data(),decB.mass);

	Double_t targetTh=prm.Tt*0.0867*0.1; //mu*g/cm^3*0.1
	Double_t BeamSpot=prm.Bs/2.355; // FWHM->sigma 
	const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr // 10 Torr!
	const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	//const Double_t AgFoil=5.44*10.473*0.1; //mu*g/cm^3*0.1
	const Double_t AgFoil=4.57*10.473*0.1; //mu*g/cm^3*0.1

	Bool_t LEHit, HEHit;

	Int_t LEHitcntr=0;

	Double_t LEeff;
	Double_t E_before_Tgt;
	Double_t E_center_Tgt;
	Double_t E_after_Tgt;
	Double_t E_before_Ag;
	Double_t E_center_Ag;

	TLorentzVector target, beam, Sys;
	TVector3 boostvect;

	TGenPhaseSpace PS0, PS1;
	Double_t wght_max,width;
	Bool_t allowed;

	// Calculate energy loss up to center of the target
	Double_t EA = prm.E;	
   	EA -= eloss(A,0.5,EA,ICWindow1,A.EL.eSi3N4, A.EL.dedxSi3N4);
   	ICdE = eloss(A,0.586,EA,ICLength,A.EL.eC4H10, A.EL.dedxC4H10);
   	EA -= ICdE;
   	EA -= eloss(A,0.5,EA,ICWindow2,A.EL.eSi3N4, A.EL.dedxSi3N4);
	E_before_Ag = EA;
	if(isAgReac){
		E_center_Ag = EA - eloss(A,47./108.,EA,AgFoil/2.,A.EL.eAg, A.EL.dedxAg);
		EA -= eloss(A,47./108.,EA,AgFoil,A.EL.eAg, A.EL.dedxAg);
   		E_before_Tgt = EA;
	}
	else{
		EA -= eloss(A,47./108.,EA,AgFoil,A.EL.eAg, A.EL.dedxAg);
   		E_before_Tgt = EA;
   		E_center_Tgt = EA - eloss(A,1.,EA,targetTh/2.,A.EL.eH, A.EL.dedxH);
   		E_after_Tgt = EA-eloss(A,1.,EA,targetTh,A.EL.eH, A.EL.dedxH);
		
		reacZ = targetTh/2.;
   		EA -= eloss(A,1.,EA,reacZ,A.EL.eH, A.EL.dedxH);
	}
	
	EA = EA/1000.; // convert to GeV for TGenPhaseSpace
	Double_t PA = sqrt(EA*EA+2*EA*mA);
	target.SetXYZT(0.0, 0.0, 0.0, ma);
	beam.SetXYZT(0.0, 0.0, PA, mA+EA);
	Sys = beam + target;
	beamE = EA*1000.;
	beamBeta = Sys.Beta();
	beamGamma = Sys.Gamma();
	beamEcm = EA*ma*1000./(mA+ma);

	if(isAgReac){
		printf("\n\nEnergy before of silver foil: %.2lf MeV\n", E_before_Ag);
		printf("\n\nEnergy at center of silver foil: %.2lf MeV\n", E_center_Ag);
		printf("\n\nEnergy after silver foil: %.2lf MeV\n", E_before_Tgt);
	}
	else{
		printf("\n\nEnergy before target: %.2lf MeV\n", E_before_Tgt);
		printf("\n\nEnergy at center of target: %.2lf MeV\n", E_center_Tgt);
		printf("\n\nEnergy at behind target: %.2lf MeV\n", E_after_Tgt);
	
		printf("\nBeta at center of target: %.3lf \n", beamBeta);
		printf("\nGamma at center of target: %.3lf \n", beamGamma);
		printf("\nCM Energy at center of target: %.2lf MeV\n\n", beamEcm);
	}

	printf("YY1 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",prm.DYY,yd.ThetaMin(prm.DYY),yd.ThetaMax(prm.DYY)); 
	printf("CsI detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",prm.DYY+11.6,csi.ThetaMin(prm.DYY+11.6),csi.ThetaMax(prm.DYY+11.6)); 
	printf("First S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",prm.DS3,sd1.ThetaMin(prm.DS3),sd1.ThetaMax(prm.DS3)); 
	printf("Second S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",prm.DS3+14.8,sd2.ThetaMin(prm.DS3+14.8),sd2.ThetaMax(prm.DS3+14.8)); 
	Double_t masses[4] = { mb, mB, mc, md};
	
	Double_t tht=0.; 
	Double_t xsec=0.; 
	Double_t xsec_chck=0.; 
	Double_t xsec_max=0.; 
	Double_t dwba_th[181]={0.}; 
	Double_t dwba_xsec[181]={0.}; 
	
	if(have_dwba_xsec==kTRUE){ 
		printf("Using DWBA cross section from %s!\n",dwbaname);
		xsec_max = load_dwba(dwbaname,dwba_th,dwba_xsec); 
	}
	
	allowed = PS0.SetDecay(Sys, prm.N, masses);
	
	if(!allowed){
		printf("Impossible decay!\n");
		printf("Exiting...\n");
		exit(0);
	}
	else{
		printf("Starting...\n");
	}

	wght_max=PS0.GetWtMax();
	printf("%lf\t%lf\n",wght_max,xsec_max);
	width = prm.W/1000.;
	printf("%lf\t%lf\n",mB,width);

	Int_t whilecount;
	// Start of event loop
	while(Evnt<nsim) 
	{
		if(isAgReac){
			reacZ = rndm->Uniform(0,AgFoil);
   			EA = E_before_Ag - eloss(A,47./108.,E_before_Ag,reacZ,A.EL.eAg, A.EL.dedxAg);
		}
		else{	
			//reacZ = rndm->Uniform(0,targetTh);
			reacZ = targetTh/2.;
   			EA = E_before_Tgt - eloss(A,1.,E_before_Tgt,reacZ,A.EL.eH, A.EL.dedxH);
 		}
		EA = EA/1000.; // convert to GeV for TGenPhaseSpace
		PA = sqrt(EA*EA+2*EA*mA);

		target.SetXYZT(0.0, 0.0, 0.0, ma);
		beam.SetXYZT(0.0, 0.0, PA, mA+EA);
		Sys = beam + target;
		boostvect = Sys.BoostVector();

		beamE = EA*1000.;
		beamBeta = Sys.Beta();
		beamGamma = Sys.Gamma();
		beamEcm = EA*ma*1000./(mA+ma);

		wght_max=PS0.GetWtMax();
		//width = prm.W/1000.;

		wght = 0.;
		clearEvt();
		mBR = rndm->BreitWigner(mB,width);
		masses[1] =mBR;
		PS0.SetDecay(Sys, prm.N, masses); //recalculate with resonance energy
			
		TLorentzVector *LTmp;
		whilecount=0;
		do{	
			wght = PS0.Generate();
			chck = rndm->Uniform(0,wght_max);
			if(have_dwba_xsec==kTRUE){
				LTmp = PS0.GetDecay(0);
				LTmp->Boost(-boostvect);
				tht = RadToDeg()*LTmp->Theta();
				xsec=eval_theta(tht,dwba_th,dwba_xsec);
				xsec_chck = rndm->Uniform(0,xsec_max);
				LTmp->Boost(boostvect);
			}
			else{
				xsec=1.;
				xsec_chck=0.;
			}
			whilecount++;
			//printf("%d\t%f\t%f\t%f\n",whilecount,tht,xsec,xsec_chck);
		}while(wght<chck||xsec<xsec_chck);

		TLorentzVector *LVb  = PS0.GetDecay(0);
		TLorentzVector *LVB  = PS0.GetDecay(1);
		TLorentzVector *LVBdec;
		TLorentzVector *LVcdec;
		TLorentzVector *LVddec;

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

		if(seqdec)
		{
			PS1.SetDecay(*LVB, seqdecN, masses2);
			do{
				wght2 = PS1.Generate();			
				chck2 = rndm->Uniform(0,1);
				LVBdec  = PS1.GetDecay(0);
				LVcdec  = PS1.GetDecay(1);
				if(seqdecN>2) LVddec  = PS1.GetDecay(2);
			}while(wght2<chck2);
			decHP.T=LVBdec->Theta();
			decHP.E=(LVBdec->E()-mBdec)*1000.;
			decHP.P=LVBdec->Phi();	
			decHP.Tdeg=RadToDeg()*decHP.T;
			decHP.Pdeg=RadToDeg()*decHP.P;
			declP1.T=LVcdec->Theta();
			declP1.E=(LVcdec->E()-mBdec)*1000.;
			declP1.P=LVcdec->Phi();	
			declP1.Tdeg=RadToDeg()*declP1.T;
			declP1.Pdeg=RadToDeg()*declP1.P;
			if(seqdecN>2){
				declP2.T=LVcdec->Theta();
				declP2.E=(LVcdec->E()-mBdec)*1000.;
				declP2.P=LVcdec->Phi();	
				declP2.Tdeg=RadToDeg()*declP2.T;
				declP2.Pdeg=RadToDeg()*declP2.P;
			}
		}
		else if(prm.N>3) // 4body
		{
			LVcdec = PS0.GetDecay(2);
			LVddec = PS0.GetDecay(3);
		
			declP1.T=LVcdec->Theta();	
			declP2.T=LVddec->Theta();
			declP1.E=(LVcdec->E()-mc)*1000.; 	
			declP2.E=(LVddec->E()-md)*1000.;
			declP1.P=LVcdec->Phi();	
			declP2.P=LVddec->Phi();	
			
			// Convert angles to degrees for root file
			declP1.Tdeg=RadToDeg()*declP1.T;
			declP2.Tdeg=RadToDeg()*declP2.T;
			declP1.Pdeg=RadToDeg()*declP1.P;
			declP2.Pdeg=RadToDeg()*declP2.P;
		}
		else if(prm.N>2) // 3body
		{
			LVcdec = PS0.GetDecay(2);
		
			declP1.T=LVcdec->Theta();	
			declP1.E=(LVcdec->E()-mc)*1000.; 	
			declP1.P=LVcdec->Phi();	
			
			// Convert angles to degrees for root file
			declP1.Tdeg=RadToDeg()*declP1.T;
			declP1.Pdeg=RadToDeg()*declP1.P;
		}

		// Position on target	
		reacX = BeamSpot*rndm->Gaus();
		reacY = BeamSpot*rndm->Gaus();
		reacPos.SetXYZ(reacX,reacY,reacZ);

		Double_t EbAg;		
		if(isAgReac){ 
			EbAg = eloss(b,47./108.,lP.E,(AgFoil-reacZ)/Cos(lP.T),b.EL.eAg,b.EL.dedxAg);	
			lP.TrgtdE = eloss(b,1.,lP.E-EbAg,targetTh/Cos(lP.T),b.EL.eH,b.EL.dedxH);	
			lP.Ebt = lP.E-EbAg-lP.TrgtdE;
		}
		else{ 
			lP.TrgtdE = eloss(b,1.,lP.E,(targetTh-reacZ)/Cos(lP.T),b.EL.eH,b.EL.dedxH);	
			lP.Ebt = lP.E-lP.TrgtdE;
		}
		
		LEHit = detHits(lP, b, reacPos);
		
		if(!seqdec){
			if(isAgReac){ 
				EbAg = eloss(B,47./108.,hP.E,(AgFoil-reacZ)/Cos(hP.T),B.EL.eAg,B.EL.dedxAg);	
				hP.TrgtdE = eloss(B,1.,hP.E-EbAg,targetTh/Cos(hP.T),B.EL.eH,B.EL.dedxH);	
				hP.Ebt = hP.E-EbAg-hP.TrgtdE;
			}
			else {
				hP.TrgtdE = eloss(B,1.,hP.E,(targetTh-reacZ)/Cos(hP.T),B.EL.eH,B.EL.dedxH);	
				hP.Ebt = hP.E-hP.TrgtdE;
			}
			HEHit = detHits(hP, B, reacPos);
		}
		else{ 
			decHP.TrgtdE = eloss(B,1.,decHP.E,(targetTh-reacZ)/Cos(decHP.T),decB.EL.eH,decB.EL.dedxH);	
			decHP.Ebt = decHP.E-decHP.TrgtdE;
			HEHit = detHits(decHP, decB, reacPos);	
			if(decc.Z>0){
				declP1.TrgtdE = eloss(B,1.,declP1.E,(targetTh-reacZ)/Cos(declP1.T),decB.EL.eH,decB.EL.dedxH);	
				declP1.Ebt = declP1.E-declP1.TrgtdE;
				detHits(declP1, decB, reacPos);
			}	
			if(seqdecN>2&&decd.Z>0){
				declP2.TrgtdE = eloss(B,1.,declP2.E,(targetTh-reacZ)/Cos(declP2.T),decB.EL.eH,decB.EL.dedxH);	
				declP2.Ebt = declP2.E-declP2.TrgtdE;
				detHits(declP2, decB, reacPos);
			}
		}
		Double_t E_before_SSB;
		if(isAgReac){
			E_before_SSB = E_before_Tgt - eloss(A,1.,E_before_Tgt,targetTh,A.EL.eH,A.EL.dedxH);
			SSBdE = eloss(A,14./28.,E_before_SSB,1000.*2.3212*0.1,B.EL.eSi,B.EL.dedxSi);
		}
		else{ 
			SSBdE = eloss(A,14./28.,E_before_SSB,1000.*2.3212*0.1,B.EL.eSi,B.EL.dedxSi);
		}
		SSBdE =rndm->Gaus(SSBdE,0.05*SSBdE);
		
		//Calculating "measured" Q-Value
		if(LEHit && yd.dE[0]>0.){
			if(csi.dE[0]>0.) LEHitcntr++;
			Double_t Eb = csi.dE[0];
			Eb= Eb+elossFi(Eb,0.1*1.4*6./Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eMy,b.EL.dedxMy); //Mylar                                                                                  
	      	Eb= Eb+elossFi(Eb,0.1*2.702*0.3/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eAl,b.EL.dedxAl); //0.3 u Al                                                                            
	      	Eb= Eb+elossFi(Eb,0.1*1.88219*0.1/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eP,b.EL.dedxP); // 0.1Phosphorus                                                                      
			Eb+= yd.dE[0]; //use measured Yd // change june28
	      	Eb= Eb+elossFi(Eb,0.1*2.32*0.35/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eSi,b.EL.dedxSi); //0.3 u Al + 1 um B equivalent in 0.35 um Si                                                            
	    	Eb= Eb+elossFi(Eb,targetTh/2./Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eH,b.EL.dedxH); //deuteron energy  in mid target midtarget
		
			Eb= Eb/1000.;

			Double_t Pb = sqrt(Eb*Eb+2.*Eb*mb);	
 			Qalt = mA + ma - mb - sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EA)*(mb+Eb)+2.*PA*Pb*cos(yd.fThetaCalc[0]*DegToRad())+2.*(EA+mA+ma-Eb-mb)*ma);
			Qalt =Qalt*1000.;
			EB_det = EA+mA+ma-Eb-mb;
			PB_det = sqrt(PA*PA+Pb*Pb-2.*PA*Pb*cos(yd.fThetaCalc[0]*DegToRad()));
			Qdet = mA+ma-mb-sqrt(EB_det*EB_det-PB_det*PB_det);
			Qdet =Qdet*1000.;
		}

		lP.Ecm = (LVb->E()-mb)*ma*1000./(mA+ma);
		hP.Ecm = (LVB->E()-mB)*ma*1000./(mA+ma);
		LVb->Boost(-boostvect);
		LVB->Boost(-boostvect);
		lP.Tcm = RadToDeg()*(Pi()-LVb->Theta());
		hP.Tcm = RadToDeg()*LVB->Theta();
		
		if(prm.N>3) // 4body
		{
			declP1.Ecm = (LVcdec->E()-mc)*ma*1000./(mA+ma);
			declP2.Ecm = (LVddec->E()-md)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			LVddec->Boost(-boostvect);
			declP1.Tcm = RadToDeg()*(Pi()-LVcdec->Theta());
			declP2.Tcm = RadToDeg()*LVddec->Theta();
		}
		else if(prm.N>2) // 4body
		{
			declP1.Ecm = (LVcdec->E()-mc)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			declP1.Tcm = RadToDeg()*(Pi()-LVcdec->Theta());
		}
		setIDet(ICdE,SSBdE);
		printf("%.6d Events processed..\r",Evnt);
		Evnt++;
		Iris->Fill();
	}

	Iris->AutoSave();
	f->Close();
	LEeff=Double_t(LEHitcntr)/Double_t(nsim)*100.;
	printf("Acceptance for target-like particles: %.1f\n",LEeff);
	Double_t time=timer.RealTime();
	printf("\nDone. %lf s\n",time);
	printf("\nOutput written to %s \n", outputname);
	
	return 0;
}
