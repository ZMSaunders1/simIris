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

#include "params.h"
#include "nucleus.h"
#include "shieldClear.h"
#include "eloss.h"
#include "YYHit.h"
#include "CsIHit.h"
#include "S3Hit.h"

using namespace TMath;

params prm;
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

Double_t Etree[2], En[2], Ecm[2], ETmp[2]; 
Double_t T[2], Td[2], Td2[2], Tcm[2]; 
Double_t P[2], Pd[2], Pd2[2]; 
Double_t TargdE[2];
Double_t Qgen, Qdet;
 	
Double_t eAH[100], eAC4H10[100], eASi3N4[100], eAAg[100];	
Double_t dedxAH[100], dedxAC4H10[100], dedxASi3N4[100], dedxAAg[100];	
Double_t eBSi[100], eBH[100], eBSiO2[100], eBB[100], eBP[100], eBAl[100], eBMy[100], eBCsI[100];	
Double_t dedxBSi[100], dedxBH[100], dedxBSiO2[100], dedxBB[100], dedxBP[100], dedxBAl[100], dedxBMy[100], dedxBCsI[100];	
Double_t ebH[100], ebSi[100], ebAl[100], ebB[100], ebMy[100], ebP[100], ebCsI[100], ebSiO2[100];	
Double_t dedxbH[100], dedxbSi[100], dedxbAl[100], dedxbB[100], dedxbMy[100], dedxbP[100], dedxbCsI[100], dedxbSiO2[100];	

void clearEvt()
{
	En[0]=0.; En[1]=0.;
	Ecm[0]=0.; Ecm[1]=0.;
	Etree[0]=0.; Etree[1]=0.;
	ETmp[0]=0.; ETmp[1]=0.;
	T[0]=0.; T[1]=0.;
	Td[0]=0.; Td[1]=0.;
	Td2[0]=0.; Td2[1]=0.;
	Tcm[0]=0.; Tcm[1]=0.;
	P[0]=0.; P[1]=0.;
	Pd[0]=0.; Pd[1]=0.;
	Pd2[0]=0.; Pd2[1]=0.;
	TargdE[0]=0.; TargdE[1]=0.;
	Qgen=0.; Qdet=0.;
	mBR=0.;
	yd.Clear();
	csi.Clear();
	sd1.Clear();
	sd2.Clear();
	
	return;
}

void loadAllELoss(std::string path, nucleus A, nucleus B, nucleus b)
{
	printf("\nLoading %s/lise_%s_in_H.txt\n",path.data(),A.name.data());	
	loadELoss(Form("%s/lise_%s_in_H.txt",path.data(),A.name.data()),eAH,dedxAH,A.mass/1000.);	
	printf("Loading %s/lise_%s_in_C4H10.txt\n",path.data(),A.name.data());	
	loadELoss(Form("%s/lise_%s_in_C4H10.txt",path.data(),A.name.data()),eAC4H10,dedxAC4H10,A.mass/1000.);	
	printf("Loading %s/lise_%s_in_Si3N4.txt\n",path.data(),A.name.data());	
	loadELoss(Form("%s/lise_%s_in_Si3N4.txt",path.data(),A.name.data()),eASi3N4,dedxASi3N4,A.mass/1000.);	
	printf("Loading %s/lise_%s_in_Ag.txt\n",path.data(),A.name.data());	
	loadELoss(Form("%s/lise_%s_in_Ag.txt",path.data(),A.name.data()),eAAg,dedxAAg,A.mass/1000.);	
	
	printf("Loading %s/lise_%s_in_Si.txt\n",path.data(),B.name.data());	
	loadELoss(Form("%s/lise_%s_in_Si.txt",path.data(),B.name.data()),eBSi,dedxBSi,B.mass/1000.);	
	printf("Loading %s/lise_%s_in_H.txt\n",path.data(),B.name.data());
	loadELoss(Form("%s/lise_%s_in_H.txt",path.data(),B.name.data()),eBH,dedxBH,B.mass/1000.);
	printf("Loading %s/lise_%s_in_SiO2.txt\n",path.data(),B.name.data());
	loadELoss(Form("%s/lise_%s_in_SiO2.txt",path.data(),B.name.data()),eBSiO2,dedxBSiO2,B.mass/1000.);
	printf("Loading %s/lise_%s_in_Al.txt\n",path.data(),B.name.data());
	loadELoss(Form("%s/lise_%s_in_Al.txt",path.data(),B.name.data()),eBAl,dedxBAl,B.mass/1000.);
	printf("Loading %s/lise_%s_in_B.txt\n",path.data(),B.name.data());
	loadELoss(Form("%s/lise_%s_in_B.txt",path.data(),B.name.data()),eBB,dedxBB,B.mass/1000.);
	printf("Loading %s/lise_%s_in_P.txt\n",path.data(),B.name.data());
	loadELoss(Form("%s/lise_%s_in_P.txt",path.data(),B.name.data()),eBP,dedxBP,B.mass/1000.);
	printf("Loading %s/lise_%s_in_My.txt\n",path.data(),B.name.data());	
	loadELoss(Form("%s/lise_%s_in_My.txt",path.data(),B.name.data()),eBMy,dedxBMy,B.mass/1000.);	
	printf("Loading %s/lise_%s_in_CsI.txt\n\n",path.data(),B.name.data());	
	loadELoss(Form("%s/lise_%s_in_CsI.txt",path.data(),B.name.data()),eBCsI,dedxBCsI,B.mass/1000.);	

	printf("Loading %s/lise_%s_in_Si.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_Si.txt",path.data(),b.name.data()),ebSi,dedxbSi,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_H.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_H.txt",path.data(),b.name.data()),ebH,dedxbH,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_Al.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_Al.txt",path.data(),b.name.data()),ebAl,dedxbAl,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_B.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_B.txt",path.data(),b.name.data()),ebB,dedxbB,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_My.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_My.txt",path.data(),b.name.data()),ebMy,dedxbMy,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_P.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_P.txt",path.data(),b.name.data()),ebP,dedxbP,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_CsI.txt\n",path.data(),b.name.data());	
	loadELoss(Form("%s/lise_%s_in_CsI.txt",path.data(),b.name.data()),ebCsI,dedxbCsI,b.mass/1000.);	
	printf("Loading %s/lise_%s_in_SiO2.txt\n\n",path.data(),b.name.data());
	loadELoss(Form("%s/lise_%s_in_SiO2.txt",path.data(),b.name.data()),ebSiO2,dedxbSiO2,b.mass/1000.);
	return;
}

Bool_t LE_ELoss(nucleus b, Double_t targetTh, TVector3 targetPos) // Calculate energy loss in detectors for light ejectile
{
	const Double_t SdThickness1=60.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t SdThickness2=1000.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t YYThickness[8]={ 112.*2.3212*0.1, 109.*2.3212*0.1, 110.*2.3212*0.1, 106.*2.3212*0.1, 101.*2.3212*0.1, 109.*2.3212*0.1, 111.*2.3212*0.1, 103.*2.3212*0.1 };
	const Double_t CsIThickness=12000.*4.51*0.1; //mu*g/cm^3*0.1

	Bool_t mask = maskClear(T[0],P[0]);
	Bool_t shield = shieldClear(T[0],P[0]);
	Bool_t YYHit = yd.Calculate(T[0],P[0],prm.DYY,targetPos) ;
	Bool_t CsIHit = csi.Calculate(T[0],P[0],prm.DYY+11.6,targetPos) ;
	Bool_t Sd1Hit = sd1.Calculate(T[0],P[0],prm.DS3,targetPos);
	Bool_t Sd2Hit = sd2.Calculate(T[0],P[0],prm.DS3+14.8,targetPos);

	TRandom3 *rndm = new TRandom3(0);
	Int_t mul;

	TargdE[0] = simEloss(b,1.,En[0],targetTh/2./Cos(T[0]),ebH,dedxbH);	
	ETmp[0] = En[0]-TargdE[0];
	
	if(mask && shield && YYHit){
		mul = yd.mul-1;
  		ETmp[0] -= simEloss(b,13./27.,ETmp[0],0.1*2.702*0.1/Cos(T[0]),ebAl, dedxbAl);
  		ETmp[0] -= simEloss(b,5./10.,ETmp[0],0.05*2.3502*0.1/Cos(T[0]),ebB, dedxbB);
  		yd.dE[mul] = simEloss(b,14./28.,ETmp[0],YYThickness[yd.Seg[mul]]/Cos(T[0]),ebSi, dedxbSi);
		ETmp[0] = ETmp[0]-yd.dE[mul];
  		yd.dE[mul] = rndm->Gaus(yd.dE[mul],0.004*yd.dE[mul]);
	}
	if(mask && shield && CsIHit){
		mul = csi.mul-1;
  		ETmp[0] -= simEloss(b,15./31.,ETmp[0],0.1*1.8219*0.1/Cos(T[0]),ebP, dedxbP);
  		ETmp[0] -= simEloss(b,13./27.,ETmp[0],0.3*2.702*0.1/Cos(T[0]),ebAl, dedxbAl);
  		ETmp[0] -= simEloss(b,100./192.,ETmp[0],6.*1.4*0.1/Cos(T[0]),ebMy, dedxbMy);
		csi.dE[mul] = simEloss(b,108./260.,ETmp[0],CsIThickness/Cos(T[0]),ebCsI, dedxbCsI);
		csi.dE[mul] = rndm->Gaus(csi.dE[mul],0.03*csi.dE[mul]);
	}
	if(shield && mask && Sd1Hit){
		mul = sd1.mul-1;
		ETmp[0] = ETmp[0] - simEloss(b,13./27.,ETmp[1],0.1*2.702*1.5/Cos(T[0]),ebAl,dedxbAl); //first metal
		ETmp[0] = ETmp[0] - simEloss(b,30./60.,ETmp[1],0.1*2.65*3.5/Cos(T[0]),ebSiO2,dedxbSiO2); //SiO2
		ETmp[0] = ETmp[0] - simEloss(b,13./27.,ETmp[1],0.1*2.702*0.3/Cos(T[0]),ebAl,dedxbAl); //second metal
		ETmp[0] = ETmp[0] - simEloss(b,5./10.,ETmp[1],0.1*2.3502*0.5/Cos(T[0]),ebB,dedxbB); //boron junction implant 		
		sd1.dE[mul] = simEloss(b,14./28.,ETmp[0],SdThickness1/Cos(T[0]),ebSi,dedxbSi);
   		ETmp[0] = ETmp[0] - sd1.dE[mul];
		sd1.dE[mul] = rndm->Gaus(sd1.dE[mul],0.01*sd1.dE[mul]);
	}
	if(shield && mask && Sd2Hit){
		mul = sd2.mul-1;
		ETmp[0] = ETmp[0] - simEloss(b,15./31.,ETmp[0],0.1*1.8219*0.5/Cos(T[0]),ebP,dedxbP); //phosphorus implant
		ETmp[0] = ETmp[0] - simEloss(b,13./27.,ETmp[0],0.1*2.702*0.3/Cos(T[0]),ebAl,dedxbAl); //metal
		ETmp[0] = ETmp[0] - simEloss(b,13./27.,ETmp[0],0.1*2.702*0.3/Cos(T[0]),ebAl,dedxbAl); //metal
		ETmp[0] = ETmp[0] - simEloss(b,15./31.,ETmp[0],0.1*1.822*0.5/Cos(T[0]),ebP,dedxbP); //phosphorus implant
		sd2.dE[mul] = simEloss(b,14./28.,ETmp[0],SdThickness2/Cos(T[0]),ebSi,dedxbSi);
		sd2.dE[mul] = rndm->Gaus(sd2.dE[mul],0.01*sd2.dE[mul]);
	}

	return (mask && shield && YYHit && CsIHit);
}

Bool_t HE_ELoss(nucleus B, Double_t targetTh, TVector3 targetPos) // Calculate energy loss in detectors for heavy ejectile
{
	const Double_t SdThickness1=60.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t SdThickness2=1000.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t YYThickness[8]={ 112.*2.3212*0.1, 109.*2.3212*0.1, 110.*2.3212*0.1, 106.*2.3212*0.1, 101.*2.3212*0.1, 109.*2.3212*0.1, 111.*2.3212*0.1, 103.*2.3212*0.1 };
	const Double_t CsIThickness=12000.*4.51*0.1; //mu*g/cm^3*0.1
	
	Bool_t mask = maskClear(T[1],P[1]);
	Bool_t shield = shieldClear(T[1],P[1]);
	Bool_t YYHit = yd.Calculate(T[1],P[1],prm.DYY,targetPos) ;
	Bool_t CsIHit = csi.Calculate(T[1],P[1],prm.DYY+11.6,targetPos) ;
	Bool_t Sd1Hit = sd1.Calculate(T[1],P[1],prm.DS3,targetPos);
	Bool_t Sd2Hit = sd2.Calculate(T[1],P[1],prm.DS3+14.8,targetPos);

	TRandom3 *rndm = new TRandom3(0);
	Int_t mul;

	TargdE[1] = simEloss(B,1.,En[1],targetTh/2./Cos(T[1]),eBH,dedxBH);
	ETmp[1] = En[1] - TargdE[1];
	
	if(mask && shield && YYHit){
		mul = yd.mul-1;
  		ETmp[1] -= simEloss(B,13./27.,ETmp[1],0.1*2.702*0.1/Cos(T[1]),eBAl, dedxBAl);
  		ETmp[1] -= simEloss(B,5./10.,ETmp[1],0.05*2.3502*0.1/Cos(T[1]),eBB, dedxBB);
  		yd.dE[mul] = simEloss(B,14./28.,ETmp[1],YYThickness[yd.Seg[mul]]/Cos(T[1]),eBSi, dedxBSi);
		ETmp[1] = ETmp[1]-yd.dE[mul];
  		yd.dE[mul] = rndm->Gaus(yd.dE[mul],0.004*yd.dE[mul]);
	}
	if(mask && shield && CsIHit){
		mul = csi.mul-1;
  		ETmp[1] -= simEloss(B,15./31.,ETmp[1],0.1*1.8219*0.1/Cos(T[1]),eBP, dedxBP);
  		ETmp[1] -= simEloss(B,13./27.,ETmp[1],0.3*2.702*0.1/Cos(T[1]),eBAl, dedxBAl);
  		ETmp[1] -= simEloss(B,100./192.,ETmp[1],6.*1.4*0.1/Cos(T[1]),eBMy, dedxBMy);
		csi.dE[mul] = simEloss(B,108./260.,ETmp[1],CsIThickness/Cos(T[1]),eBCsI, dedxBCsI);
		csi.dE[mul] = rndm->Gaus(csi.dE[mul],0.03*csi.dE[mul]);
	}
	if(shield && mask && Sd1Hit){
		mul = sd1.mul-1;
		ETmp[1] = ETmp[1] - simEloss(B,13./27.,ETmp[1],0.1*2.702*1.5/Cos(T[1]),eBAl,dedxBAl); //first metal
		ETmp[1] = ETmp[1] - simEloss(B,30./60.,ETmp[1],0.1*2.65*3.5/Cos(T[1]),eBSiO2,dedxBSiO2); //SiO2
		ETmp[1] = ETmp[1] - simEloss(B,13./27.,ETmp[1],0.1*2.702*0.3/Cos(T[1]),eBAl,dedxBAl); //second metal
		ETmp[1] = ETmp[1] - simEloss(B,5./10.,ETmp[1],0.1*2.3502*0.5/Cos(T[1]),eBB,dedxBB); //boron junction implant 		
		printf("%s ",B.name.data());
		sd1.dE[mul] = simEloss(B,14./28.,ETmp[1],SdThickness1/Cos(T[1]),eBSi,dedxBSi);
   		ETmp[1] = ETmp[1] - sd1.dE[mul];
		sd1.dE[mul] = rndm->Gaus(sd1.dE[mul],0.01*sd1.dE[mul]);
	}
	if(shield && mask && Sd2Hit){
		mul = sd2.mul-1;
		ETmp[1] = ETmp[1] - simEloss(B,15./31.,ETmp[1],0.1*1.8219*0.5/Cos(T[1]),eBP,dedxBP); //phosphorus implant
		ETmp[1] = ETmp[1] - simEloss(B,13./27.,ETmp[1],0.1*2.702*0.3/Cos(T[1]),eBAl,dedxBAl); //metal
		ETmp[1] = ETmp[1] - simEloss(B,13./27.,ETmp[1],0.1*2.702*0.3/Cos(T[1]),eBAl,dedxBAl); //metal
		ETmp[1] = ETmp[1] - simEloss(B,15./31.,ETmp[1],0.1*1.822*0.5/Cos(T[1]),eBP,dedxBP); //phosphorus implant
		sd2.dE[mul] = simEloss(B,14./28.,ETmp[1],SdThickness2/Cos(T[1]),eBSi,dedxBSi);
		sd2.dE[mul] = rndm->Gaus(sd2.dE[mul],0.01*sd2.dE[mul]);
	}

	return (shield && mask && Sd1Hit && Sd2Hit);
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
	
	YYHit *ipyd = &yd; 
	CsIHit *ipcsi = &csi; 
	S3Hit *ipsd1 = &sd1; 
	S3Hit *ipsd2 = &sd2; 
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

	std::string dedxstr = dedxpath;
	if(!seqdec) loadAllELoss(dedxstr,A,B,b);
	else loadAllELoss(dedxstr,A,decB,b);

	Double_t targetTh=prm.Tt*0.0867*0.1; //mu*g/cm^3*0.1
	Double_t BeamSpot=prm.Bs/2.355; // FWHM->sigma 
	const Double_t ICLength=22.9*0.062; //cm*mg/cm^3 at 19.5 Torr
	const Double_t ICWindow1=0.03*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t ICWindow2=0.05*3.44*0.1; //mu*g/cm^3*0.1
	const Double_t AgFoil=5.44*10.473*0.1; //mu*g/cm^3*0.1

	Bool_t LEHit, HEHit;

	// Calculate energy loss up to center of the target
	Double_t EA = prm.E;	
   	EA -= simEloss(A,0.5,EA,ICWindow1,eASi3N4, dedxASi3N4);
   	EA -= simEloss(A,0.586,EA,ICLength,eAC4H10, dedxAC4H10);
   	EA -= simEloss(A,0.5,EA,ICWindow2,eASi3N4, dedxASi3N4);
   	EA -= simEloss(A,47./108.,EA,AgFoil,eAAg, dedxAAg);
   	EA -= simEloss(A,1.,EA,targetTh/2.,eAH, dedxAH);
 
	EA = EA/1000.; // convert to GeV for TGenPhaseSpace
	Double_t PA = sqrt(EA*EA+2*EA*mA);

	TLorentzVector target(0.0, 0.0, 0.0, ma);
	TLorentzVector beam(0.0, 0.0, PA, mA+EA);
	TLorentzVector Sys = beam + target;
	TVector3 boostvect = Sys.BoostVector();

	printf("\n\nEnergy at center of target: %.2lf MeV\n", EA*1000.);
	printf("\nBeta at center of target: %.3lf \n", Sys.Beta());
	printf("\nGamma at center of target: %.3lf \n", Sys.Gamma());
	printf("\nCM Energy at center of target: %.2lf MeV\n\n", EA*ma*1000./(mA+ma));

	Double_t masses[4] = { mb, mB, mc, md};

	TGenPhaseSpace PS0, PS1;
	Bool_t allowed = PS0.SetDecay(Sys, prm.N, masses);

	if(!allowed){
		printf("Impossible decay!\n");
		printf("Exiting...\n");
		exit(0);
	}

	// Set up output file and tree
	TFile *f = new TFile(outputname,"RECREATE");
	TTree *iris = new TTree("iris","iris simulation");
	
	iris->Branch("Evnt",&Evnt,"Evnt/I"); 
	iris->Branch("E",Etree,"Etree[2]/D"); 
	iris->Branch("E2",En,"En[2]/D"); 
	iris->Branch("Ecm",Ecm,"Ecm[2]/D"); 
	iris->Branch("Theta",Td,"Td[2]/D"); 
	iris->Branch("Theta2",Td2,"Td2[2]/D"); 
	iris->Branch("Thetacm",Tcm,"Tcm[2]/D"); 
	iris->Branch("Phi",Pd,"Pd[2]/D"); 
	iris->Branch("Phi2",Pd2,"Pd2[2]/D"); 
	iris->Branch("TargdE",TargdE,"TargdE[2]/D"); 
	//iris->Branch("csi.dE",&csi.dE,"csi.dE/D"); 
	iris->Branch("wght",&wght,"wght/D"); 
	iris->Branch("Qgen",&Qgen,"Qgen/D"); 
	iris->Branch("Qdet",&Qdet,"Qdet/D"); 
	iris->Branch("mBR",&mBR,"mBR/D"); 
	iris->Branch("yd",&ipyd,32000,99); 
	iris->Branch("csi",&ipcsi,32000,99); 
	iris->Branch("sd1",&ipsd1,32000,99); 
	iris->Branch("sd2",&ipsd2,32000,99); 

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
	
		T[0]=LVb->Theta();	
		T[1]=LVB->Theta();
		En[0]=(LVb->E()-mb)*1000.; 	
		En[1]=(LVB->E()-mB)*1000.;
		P[0]=LVb->Phi();	
		P[1]=LVB->Phi();	
		
		// Convert angles to degrees for root file
		Td[0]=RadToDeg()*T[0];
		Td[1]=RadToDeg()*T[1];
		Pd[0]=RadToDeg()*P[0];
		Pd[1]=RadToDeg()*P[1];
		Etree[0] = En[0];
		Etree[1] = En[1];

		if(seqdec)
		{
			PS1.SetDecay(*LVB, seqdecN, masses2);
			do{
				wght2 = PS1.Generate();			
				chck2 = rndm->Uniform(0,1);
				LVBdec  = PS1.GetDecay(0);
			}while(wght2<chck2);
			T[1]=LVBdec->Theta();
			En[1]=(LVBdec->E()-mBdec)*1000.;
			P[1]=LVBdec->Phi();	
			Td2[0]=RadToDeg()*T[0];
			Td2[1]=RadToDeg()*T[1];
			Pd2[0]=RadToDeg()*P[0];
			Pd2[1]=RadToDeg()*P[1];
		}
		else
		{
			Td2[0]=Td[0];
			Td2[1]=Td[1];
			Pd2[0]=Pd[0];
			Pd2[1]=Pd[1];
		}	
	
		// Position on target	
		Double_t  targetPosX = BeamSpot*rndm->Gaus();
		Double_t  targetPosY = BeamSpot*rndm->Gaus();
		TVector3 targetPos(targetPosX,targetPosY,0);
		
		LEHit = LE_ELoss(b, targetTh, targetPos);
		
		if(!seqdec){
			HEHit = HE_ELoss(B, targetTh, targetPos);	
		}
		else{ 
			HEHit = HE_ELoss(decB, targetTh, targetPos);
		}
		
		if(LEHit && yd.dE[0]>0.){
			Double_t Pb = LVb->P();
			Double_t Eb = LVb->E()-mb;	
 			Qdet = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EA)*(mb+Eb)+2.*PA*Pb*cos(yd.fThetaCalc[0]*DegToRad())+2.*(EA+mA+ma-Eb-mb)*ma);
			Qdet =Qdet*1000.;
		}

		Ecm[0] = (LVb->E()-mb)*ma*1000./(mA+ma);
		Ecm[1] = (LVB->E()-mB)*ma*1000./(mA+ma);
		LVb->Boost(-boostvect);
		LVB->Boost(-boostvect);
		Tcm[0] = RadToDeg()*(Pi()-LVb->Theta());
		Tcm[1] = RadToDeg()*LVB->Theta();

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
