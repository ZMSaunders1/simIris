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

reacParams reacPrm;
geoParams geoPrm;
YYHit yd;
CsIHit csi;
S3Hit sd1, sd2;
PTrack blP, tlP;
PTrack blPdec, tlPdec1, tlPdec2;
IDet det;

TGenPhaseSpace PS0, PS1;
TLorentzVector LorVb;
TLorentzVector LorVB;
TLorentzVector LorVBdec;
TLorentzVector LorVcdec;
TLorentzVector LorVddec;

Double_t mA=0.;	
Double_t ma=0.;	
Double_t mB=0.;	
Double_t mBR=0.;	
Double_t mb=0.;	
Double_t mc=0.;
Double_t md=0.;

Double_t Qgen=sqrt(-1.);
Double_t Qdet=sqrt(-1.);
Double_t EB_det = sqrt(-1.);
Double_t PB_det = sqrt(-1.);
Double_t beamE=0.;
Double_t beamBeta=0.;
Double_t beamGamma=0.;
Double_t beamEcm=0.;
TVector3 reacPos;
Double_t SSBdE=0.;

void clearEvt()
{
	//TargdE[0]=0.; TargdE[1]=0.;
	Qgen=sqrt(-1.); Qdet=sqrt(-1.);
	EB_det=sqrt(-1.); PB_det=sqrt(-1.);
	mBR=0.;
	tlP.Clear();
	blP.Clear();
	blPdec.Clear();
	tlPdec1.Clear();
	tlPdec2.Clear();
	SSBdE=0.;
	yd.Clear();
	csi.Clear();
	sd1.Clear();
	sd2.Clear();
	det.Clear();
	LorVb.Clear();
	LorVB.Clear();
	LorVBdec.Clear();
	LorVcdec.Clear();
	LorVddec.Clear();
	PS0.Clear();
	PS1.Clear();
	
	return;
}

int main(int argc, char *argv[])
{
	Bool_t have_reaction=kFALSE;
	Bool_t have_geometry=kFALSE;
	Bool_t have_dwba_xsec=kFALSE;
	char *endptr;
	Int_t nsim = 1E3;
	char *reacParamsname =NULL;
	char *geoParamsname =NULL;
	char *dedxpath =NULL;
	char *outputname =NULL;
	char *dwbaname =NULL;
	Bool_t isSHTReac = kFALSE;

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
			else if(strncmp(argv[i],"--reaction=",11)==0){
				have_reaction=kTRUE;
				reacParamsname = argv[i]+11;
			}
			else if(strncmp(argv[i],"--geo=",6)==0){
				have_geometry=kTRUE;
				geoParamsname = argv[i]+6;
			}
			else if(strncmp(argv[i],"--dwba=",7)==0){
				have_dwba_xsec=kTRUE;
				dwbaname = argv[i]+7;
			}
			else if(strncmp(argv[i],"--events=",9)==0){
	  			nsim = strtol(argv[i]+9,&endptr,10);//converting string to number
			}
		}
	}
	TStopwatch timer;
	timer.Start();
		
	if(have_reaction==kTRUE){
		printf("Loading parameters from file %s.\n",reacParamsname);
		reacPrm.Load(reacParamsname);
	}
	reacPrm.Print();
	
	isSHTReac=reacPrm.SHT;
	
	if(have_geometry==kTRUE){
		printf("Loading parameters from file %s.\n",geoParamsname);
		geoPrm.Load(geoParamsname);
	}
	geoPrm.Print();

	TRandom3 *rndm = new TRandom3(0);
	Int_t Evnt=0; 
	Double_t chck=0.;
	Double_t chck2=0.;
	Double_t wght=0.;
	Double_t wght2=0.;

	Double_t ICdE;	
	YYHit *ipyd = &yd; 
	CsIHit *ipcsi = &csi; 
	S3Hit *ipsd1 = &sd1; 
	S3Hit *ipsd2 = &sd2; 
	PTrack *iptlP = &tlP;
	PTrack *ipblP = &blP;
	PTrack *ipblPdec = &blPdec;
	PTrack *iptlPdec1 = &tlPdec1;
	PTrack *iptlPdec2 = &tlPdec2;
	IDet *ipdet = &det;
	TLorentzVector *LVb = &LorVb;	
	TLorentzVector *LVB = &LorVB;	
	TLorentzVector *LVBdec = &LorVBdec;	
	TLorentzVector *LVcdec = &LorVcdec;	
	TLorentzVector *LVddec = &LorVddec;	
	
	nucleus A, a, B, b, c, d, decB,decc,decd;
	Double_t reacX, reacY, reacZ;

	A.getInfo(binpath, reacPrm.A);
	a.getInfo(binpath, reacPrm.a);
	B.getInfo(binpath, reacPrm.B);
	b.getInfo(binpath, reacPrm.b);
	c.getInfo(binpath, reacPrm.c);
	d.getInfo(binpath, reacPrm.d);

	mA = A.mass/1000.;	
	ma = a.mass/1000.;	
	mB = B.mass/1000.+reacPrm.R/1000.;	
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
	printf("\nResonance Energy: %.2lf\tlowest threshold: %.2lf\n",reacPrm.R,S_low);

	if(S_low<reacPrm.R){
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
	// Set up output file and tree
	TFile *f = new TFile(outputname,"RECREATE");
	TTree *Iris = new TTree("Iris","Iris simulation");
	
	Iris->Branch("Evnt",&Evnt,"Evnt/I"); 
	Iris->Branch("beamE",&beamE,"beamE/D"); 
	Iris->Branch("beamBeta",&beamBeta,"beamBeta/D"); 
	Iris->Branch("beamGamma",&beamGamma,"beamGamma/D"); 
	Iris->Branch("beamEcm",&beamEcm,"beamEcm/D"); 
	Iris->Branch("reacX",&reacX,"reacX/D"); 
	Iris->Branch("reacY",&reacY,"reacY/D"); 
	Iris->Branch("reacZ",&reacZ,"reacZ/D"); 
	Iris->Branch("tlP.",&iptlP,32000,99); 
	Iris->Branch("blP.",&ipblP,32000,99); 
	Iris->Branch("blPdec.",&ipblPdec,32000,99); 
	Iris->Branch("tlPdec1.",&iptlPdec1,32000,99); 
	Iris->Branch("tlPdec2.",&iptlPdec2,32000,99); 
	Iris->Branch("wght",&wght,"wght/D"); 
	Iris->Branch("Qgen",&Qgen,"Qgen/D"); 
	Iris->Branch("Qdet",&Qdet,"Qdet/D"); 
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
	A.EL.loadIncomingELoss(dedxstr,A.name.data(),geoPrm.MFoil,geoPrm.MTgt,A.mass);
	b.EL.loadOutgoingELoss(dedxstr,b.name.data(),geoPrm.MFoil,geoPrm.MTgt,b.mass);
	if(!seqdec) B.EL.loadOutgoingELoss(dedxstr,B.name.data(),geoPrm.MFoil,geoPrm.MTgt,B.mass);
	else{
	   	decB.EL.loadOutgoingELoss(dedxstr,decB.name.data(),geoPrm.MFoil,geoPrm.MTgt,decB.mass);
	   	if(decc.Z>0) decc.EL.loadOutgoingELoss(dedxstr,decc.name.data(),geoPrm.MFoil,geoPrm.MTgt,decc.mass);
	   	if(seqdecN>2&&decd.Z>0) decd.EL.loadOutgoingELoss(dedxstr,decd.name.data(),geoPrm.MFoil,geoPrm.MTgt,decd.mass);
	}

	Double_t FoilTh=geoPrm.TFoil; //mu*g/cm^3*0.1
	Double_t targetTh=geoPrm.TTgt; //mu*g/cm^3*0.1
	Double_t BeamSpot=geoPrm.Bs/2.355; // FWHM->sigma 
	const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr 
	const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	
	//tlP.nuc = b;
	//blP.nuc = B;

	yd.Init(geoPrm.TYY);
	sd1.Init(0,geoPrm.TS3[0]);
	sd2.Init(1,geoPrm.TS3[1]);

	Bool_t LEHit, HEHit;

	Int_t LEHitcntr=0;
	Int_t HEHitcntr=0;

	Double_t LEeff, HEeff;
	Double_t E_before_Tgt=0.;
	Double_t E_center_Tgt=0.;
	Double_t E_after_Tgt=0.;
	Double_t E_before_Ag=0.;
	Double_t E_center_Ag=0.;
	Double_t E_after_Ag=0.;

	TLorentzVector target, beam, Sys;
	TVector3 boostvect;

	Double_t wght_max,width;
	Bool_t allowed;

	// Calculate energy loss up to center of the target
	Double_t EA = reacPrm.E;	
   	EA -= eloss(A,0.5,EA,ICWindow1,A.EL.eSi3N4, A.EL.dedxSi3N4);
   	ICdE = eloss(A,0.586,EA,ICLength,A.EL.eC4H10, A.EL.dedxC4H10);
   	EA -= ICdE;
   	EA -= eloss(A,0.5,EA,ICWindow2,A.EL.eSi3N4, A.EL.dedxSi3N4);
	E_before_Ag = EA;
	if(!isSHTReac){
		E_center_Ag = EA - eloss(A,47./108.,EA,FoilTh/2.,A.EL.eFoil, A.EL.dedxFoil);
		EA -= eloss(A,47./108.,EA,FoilTh,A.EL.eFoil, A.EL.dedxFoil);
   		E_before_Tgt = EA;
	}
	else{
		EA -= eloss(A,47./108.,EA,FoilTh,A.EL.eFoil, A.EL.dedxFoil);
   		E_after_Ag = EA;
   		E_before_Tgt = EA;
   		E_center_Tgt = EA - eloss(A,1.,EA,targetTh/2.,A.EL.eTgt, A.EL.dedxTgt);
   		E_after_Tgt = EA-eloss(A,1.,EA,targetTh,A.EL.eTgt, A.EL.dedxTgt);
		
		reacZ = targetTh/2.;
   		EA -= eloss(A,1.,EA,reacZ,A.EL.eTgt, A.EL.dedxTgt);
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

	if(!isSHTReac){
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

	printf("YY1 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DYY,yd.ThetaMin(geoPrm.DYY),yd.ThetaMax(geoPrm.DYY)); 
	printf("CsI detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DYY+11.6,csi.ThetaMin(geoPrm.DYY+11.6),csi.ThetaMax(geoPrm.DYY+11.6)); 
	printf("First S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DS3,sd1.ThetaMin(geoPrm.DS3),sd1.ThetaMax(geoPrm.DS3)); 
	printf("Second S3 detector at distance of %.1lf mm from target, covering theta range from %.2lf to %.2lf\n",geoPrm.DS3+14.8,sd2.ThetaMin(geoPrm.DS3+14.8),sd2.ThetaMax(geoPrm.DS3+14.8)); 
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
	
	allowed = PS0.SetDecay(Sys, reacPrm.N, masses);
	
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
	width = reacPrm.W/1000.;
	printf("%lf\t%lf\n",mB,width);

	Int_t whilecount;
	// Start of event loop
	while(Evnt<nsim) 
	{
		if(!isSHTReac){
			reacZ = FoilTh/2.;
			//reacZ = rndm->Uniform(0,FoilTh);
   			EA = E_before_Ag - eloss(A,47./108.,E_before_Ag,reacZ,A.EL.eFoil, A.EL.dedxFoil);
		}
		else{	
			//reacZ = rndm->Uniform(0,targetTh);
			reacZ = targetTh/2.;
   			EA = E_before_Tgt - eloss(A,1.,E_before_Tgt,reacZ,A.EL.eTgt, A.EL.dedxTgt);
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

		//width = reacPrm.W/1000.;

		wght = 0.;
		clearEvt();
		mBR = rndm->BreitWigner(mB,width);
		masses[1] =mBR;
		PS0.SetDecay(Sys, reacPrm.N, masses); //recalculate with resonance energy
		wght_max=PS0.GetWtMax();
	
		TLorentzVector *LTmp;
		whilecount=0;
		do{	
			wght = PS0.Generate();
			if(wght!=wght) continue; // catch NaNs
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

		LVb  = PS0.GetDecay(0);
		LVB  = PS0.GetDecay(1);

		TLorentzVector Frag  = Sys-*LVb;

		Qgen= (mA+ma-mb-Frag.M())*1000.;	
	
		tlP.T=LVb->Theta();	
		blP.T=LVB->Theta();
		tlP.E=(LVb->E()-mb)*1000.; 	
		blP.E=(LVB->E()-mB)*1000.;
		tlP.P=LVb->Phi();	
		blP.P=LVB->Phi();	
		
		// Convert angles to degrees for root file
		tlP.Tdeg=RadToDeg()*tlP.T;
		blP.Tdeg=RadToDeg()*blP.T;
		tlP.Pdeg=RadToDeg()*tlP.P;
		blP.Pdeg=RadToDeg()*blP.P;

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
			blPdec.T=LVBdec->Theta();
			blPdec.E=(LVBdec->E()-mBdec)*1000.;
			blPdec.P=LVBdec->Phi();	
			blPdec.Tdeg=RadToDeg()*blPdec.T;
			blPdec.Pdeg=RadToDeg()*blPdec.P;
			tlPdec1.T=LVcdec->Theta();
			tlPdec1.E=(LVcdec->E()-mcdec)*1000.;
			tlPdec1.P=LVcdec->Phi();	
			tlPdec1.Tdeg=RadToDeg()*tlPdec1.T;
			tlPdec1.Pdeg=RadToDeg()*tlPdec1.P;
			if(seqdecN>2){
				tlPdec2.T=LVcdec->Theta();
				tlPdec2.E=(LVcdec->E()-mddec)*1000.;
				tlPdec2.P=LVcdec->Phi();	
				tlPdec2.Tdeg=RadToDeg()*tlPdec2.T;
				tlPdec2.Pdeg=RadToDeg()*tlPdec2.P;
			}
		}
		else if(reacPrm.N>3) // 4body
		{
			LVcdec = PS0.GetDecay(2);
			LVddec = PS0.GetDecay(3);
		
			tlPdec1.T=LVcdec->Theta();	
			tlPdec2.T=LVddec->Theta();
			tlPdec1.E=(LVcdec->E()-mc)*1000.; 	
			tlPdec2.E=(LVddec->E()-md)*1000.;
			tlPdec1.P=LVcdec->Phi();	
			tlPdec2.P=LVddec->Phi();	
			
			// Convert angles to degrees for root file
			tlPdec1.Tdeg=RadToDeg()*tlPdec1.T;
			tlPdec2.Tdeg=RadToDeg()*tlPdec2.T;
			tlPdec1.Pdeg=RadToDeg()*tlPdec1.P;
			tlPdec2.Pdeg=RadToDeg()*tlPdec2.P;
		}
		else if(reacPrm.N>2) // 3body
		{
			LVcdec = PS0.GetDecay(2);
		
			tlPdec1.T=LVcdec->Theta();	
			tlPdec1.E=(LVcdec->E()-mc)*1000.; 	
			tlPdec1.P=LVcdec->Phi();	
			
			// Convert angles to degrees for root file
			tlPdec1.Tdeg=RadToDeg()*tlPdec1.T;
			tlPdec1.Pdeg=RadToDeg()*tlPdec1.P;
		}

		// Position on target	
		reacX = BeamSpot*rndm->Gaus();
		reacY = BeamSpot*rndm->Gaus();
		reacPos.SetXYZ(reacX,reacY,reacZ);

		if(!isSHTReac){ 
			tlP.AgdE = eloss(b,47./108.,tlP.E,(FoilTh-reacZ)/Cos(tlP.T),b.EL.eFoil,b.EL.dedxFoil);	
			tlP.TrgtdE = eloss(b,1.,tlP.E-tlP.AgdE,targetTh/Cos(tlP.T),b.EL.eTgt,b.EL.dedxTgt);	
			tlP.Ebt = tlP.E-tlP.AgdE-tlP.TrgtdE;
		}
		else{
		   	tlP.AgdE = 0.;	
			tlP.TrgtdE = eloss(b,1.,tlP.E,(targetTh-reacZ)/Cos(tlP.T),b.EL.eTgt,b.EL.dedxTgt);	
			tlP.Ebt = tlP.E-tlP.TrgtdE;
		}
		
		LEHit = detHits(tlP, b, reacPos,geoPrm.Mask,geoPrm.Shield);
		
		if(!seqdec){
			if(!isSHTReac){ 
				blP.AgdE = eloss(B,47./108.,blP.E,(FoilTh-reacZ)/Cos(blP.T),B.EL.eFoil,B.EL.dedxFoil);	
				blP.TrgtdE = eloss(B,1.,blP.E-blP.AgdE,targetTh/Cos(blP.T),B.EL.eTgt,B.EL.dedxTgt);	
			}
			else {
		   		blP.AgdE = 0.;	
				blP.TrgtdE = eloss(B,1.,blP.E,(targetTh-reacZ)/Cos(blP.T),B.EL.eTgt,B.EL.dedxTgt);	
			}
			blP.Ebt = blP.E-blP.AgdE-blP.TrgtdE;
			HEHit = detHits(blP, B, reacPos,geoPrm.Mask,geoPrm.Shield);
		}
		else{ 
		   	blPdec.AgdE = 0.;	
			blPdec.TrgtdE = eloss(decB,1.,blPdec.E,(targetTh-reacZ)/Cos(blPdec.T),decB.EL.eTgt,decB.EL.dedxTgt);	
			blPdec.Ebt = blPdec.E-blPdec.TrgtdE;
			HEHit = detHits(blPdec, decB, reacPos,geoPrm.Mask,geoPrm.Shield);	
			if(decc.Z>0){
		   		tlPdec1.AgdE = 0.;	
				tlPdec1.TrgtdE = eloss(decc,1.,tlPdec1.E,(targetTh-reacZ)/Cos(tlPdec1.T),decc.EL.eTgt,decc.EL.dedxTgt);	
				tlPdec1.Ebt = tlPdec1.E-tlPdec1.TrgtdE;
				detHits(tlPdec1, decc, reacPos,geoPrm.Mask,geoPrm.Shield);
			}	
			if(seqdecN>2&&decd.Z>0){
		   		tlPdec2.AgdE = 0.;	
				tlPdec2.TrgtdE = eloss(decd,1.,tlPdec2.E,(targetTh-reacZ)/Cos(tlPdec2.T),decd.EL.eTgt,decd.EL.dedxTgt);	
				tlPdec2.Ebt = tlPdec2.E-tlPdec2.TrgtdE;
				detHits(tlPdec2, decd, reacPos,geoPrm.Mask,geoPrm.Shield);
			}
		}
		if(!isSHTReac){
			//E_before_SSB = E_before_Tgt - eloss(A,1.,E_after_Ag,targetTh,A.EL.eTgt,A.EL.dedxTgt);
			SSBdE = eloss(A,14./28.,E_after_Ag,500.*2.3212*0.1,B.EL.eSi,B.EL.dedxSi);
		}
		else{ 
			SSBdE = eloss(A,14./28.,E_after_Tgt,500.*2.3212*0.1,B.EL.eSi,B.EL.dedxSi);
		}
		SSBdE =rndm->Gaus(SSBdE,0.05*SSBdE);
		sortEnergies(); // sort detector hits by energy

		//Calculating "measured" Q-Value
		if(LEHit && yd.dE.size()>0 && csi.dE.size()>0){
			if(csi.dE[0]>0. && yd.dE[0]>0.){
			   	LEHitcntr++;
				Double_t Eb = csi.dE[0];
				Eb= Eb+elossFi(Eb,0.1*1.4*6./Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eMy,b.EL.dedxMy); //Mylar
	      		Eb= Eb+elossFi(Eb,0.1*2.702*0.3/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eAl,b.EL.dedxAl); //0.3 u Al
	      		Eb= Eb+elossFi(Eb,0.1*1.88219*0.1/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eP,b.EL.dedxP); // 0.1Phosphorus
				Eb+= yd.dE[0]; //use measured Yd // change june28
	      		Eb= Eb+elossFi(Eb,0.1*2.32*0.35/Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eSi,b.EL.dedxSi); //0.3 u Al + 1 um B equivalent in 0.35 um Si
	    		Eb= Eb+elossFi(Eb,targetTh/2./Cos(yd.fThetaCalc[0]*DegToRad()),b.EL.eTgt,b.EL.dedxTgt); //deuteron energy  in mid target midtarget
		
				Eb= Eb/1000.;

				Double_t Pb = sqrt(Eb*Eb+2.*Eb*mb);	
				EB_det = EA+mA+ma-Eb-mb;
				PB_det = sqrt(PA*PA+Pb*Pb-2.*PA*Pb*cos(yd.fThetaCalc[0]*DegToRad()));
				Qdet = mA+ma-mb-sqrt(EB_det*EB_det-PB_det*PB_det);
				Qdet =Qdet*1000.;
			}
		}

		if(HEHit && sd1.dE.size()>0 && sd2.dE.size()>0.){
			if(sd1.dE[0]>0. && sd2.dE[0]>0.) HEHitcntr++;
		}
		tlP.Ecm = (LVb->E()-mb)*ma*1000./(mA+ma);
		blP.Ecm = (LVB->E()-mB)*ma*1000./(mA+ma);
		LVb->Boost(-boostvect);
		LVB->Boost(-boostvect);
		tlP.Tcm = RadToDeg()*(Pi()-LVb->Theta());
		blP.Tcm = RadToDeg()*LVB->Theta();
		
		if(seqdec)
		{
			blPdec.Ecm = (LVBdec->E()-mBdec)*ma*1000./(mA+ma);
			LVBdec->Boost(-boostvect);
			blPdec.Tcm = RadToDeg()*LVBdec->Theta();
			tlPdec1.Ecm = (LVcdec->E()-mcdec)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			tlPdec1.Tcm = RadToDeg()*(Pi()-LVcdec->Theta());
			if(seqdecN>2){
				tlPdec2.Ecm = (LVddec->E()-mddec)*ma*1000./(mA+ma);
				LVddec->Boost(-boostvect);
				tlPdec2.Tcm = RadToDeg()*LVddec->Theta();
			}
		}
		else if(reacPrm.N>3) // 4body
		{
			tlPdec1.Ecm = (LVcdec->E()-mc)*ma*1000./(mA+ma);
			tlPdec2.Ecm = (LVddec->E()-md)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			LVddec->Boost(-boostvect);
			tlPdec1.Tcm = RadToDeg()*(Pi()-LVcdec->Theta());
			tlPdec2.Tcm = RadToDeg()*LVddec->Theta();
		}
		else if(reacPrm.N>2) // 3body
		{
			tlPdec1.Ecm = (LVcdec->E()-mc)*ma*1000./(mA+ma);
			LVcdec->Boost(-boostvect);
			tlPdec1.Tcm = RadToDeg()*(Pi()-LVcdec->Theta());
		}
		setIDet(ICdE,SSBdE);
		printf("Writing %s: %.6d of %.6d events processed..\r",outputname,Evnt,nsim);
		Evnt++;
		Iris->Fill();
	}

	Iris->AutoSave();
	f->Close();
	LEeff=Double_t(LEHitcntr)/Double_t(nsim)*100.;
	HEeff=Double_t(HEHitcntr)/Double_t(nsim)*100.;
	printf("Acceptance for target-like particles: %.1f\n",LEeff);
	printf("Acceptance for beam-like particles: %.1f\n",HEeff);
	Double_t time=timer.RealTime();
	printf("\nDone. %lf s\n",time);
	printf("\nOutput written to %s \n", outputname);
	
	return 0;
}
