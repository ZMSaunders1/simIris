#include "header.h"
 	
Double_t eAH[100], eAC4H10[100], eASi3N4[100], eAAg[100];	
Double_t dedxAH[100], dedxAC4H10[100], dedxASi3N4[100], dedxAAg[100];	
Double_t eBSi[100], eBH[100], eBSiO2[100], eBB[100], eBP[100], eBAl[100], eBMy[100], eBCsI[100];	
Double_t dedxBSi[100], dedxBH[100], dedxBSiO2[100], dedxBB[100], dedxBP[100], dedxBAl[100], dedxBMy[100], dedxBCsI[100];	
Double_t ebH[100], ebSi[100], ebAl[100], ebB[100], ebMy[100], ebP[100], ebCsI[100], ebSiO2[100];	
Double_t dedxbH[100], dedxbSi[100], dedxbAl[100], dedxbB[100], dedxbMy[100], dedxbP[100], dedxbCsI[100], dedxbSiO2[100];	

Double_t TargdE[2];

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

Bool_t detHitsL(PTrack tr, nucleus b, Double_t targetTh, TVector3 targetPos) // Calculate energy loss in detectors for light ejectile
{
	const Double_t SdThickness1=60.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t SdThickness2=1000.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t YYThickness[8]={ 112.*2.3212*0.1, 109.*2.3212*0.1, 110.*2.3212*0.1, 106.*2.3212*0.1, 101.*2.3212*0.1, 109.*2.3212*0.1, 111.*2.3212*0.1, 103.*2.3212*0.1 };
	const Double_t CsIThickness=12000.*4.51*0.1; //mu*g/cm^3*0.1

	Bool_t mask = maskClear(tr.T,tr.P);
	Bool_t shield = shieldClear(tr.T,tr.P);
	Bool_t YYHit = yd.Calculate(tr.T,tr.P,prm.DYY,targetPos) ;
	Bool_t CsIHit = csi.Calculate(tr.T,tr.P,prm.DYY+11.6,targetPos) ;
	Bool_t Sd1Hit = sd1.Calculate(tr.T,tr.P,prm.DS3,targetPos);
	Bool_t Sd2Hit = sd2.Calculate(tr.T,tr.P,prm.DS3+14.8,targetPos);

	Double_t ETmp = 0.;

	TRandom3 *rndm = new TRandom3(0);
	Int_t mul;

	TargdE[0] = simEloss(b,1.,tr.E,targetTh/2./Cos(tr.T),ebH,dedxbH);	
	ETmp = tr.E-TargdE[0];
	
	if(mask && shield && YYHit){
		mul = yd.mul-1;
  		ETmp -= simEloss(b,13./27.,ETmp,0.1*2.702*0.1/Cos(tr.T),ebAl, dedxbAl);
  		ETmp -= simEloss(b,5./10.,ETmp,0.05*2.3502*0.1/Cos(tr.T),ebB, dedxbB);
  		yd.dE[mul] = simEloss(b,14./28.,ETmp,YYThickness[yd.Seg[mul]]/Cos(tr.T),ebSi, dedxbSi);
		ETmp = ETmp-yd.dE[mul];
  		yd.dE[mul] = rndm->Gaus(yd.dE[mul],0.00225*5.73*Sqrt(5.73/yd.dE[mul]));
	}
	if(mask && shield && CsIHit){
		mul = csi.mul-1;
  		ETmp -= simEloss(b,15./31.,ETmp,0.1*1.8219*0.1/Cos(tr.T),ebP, dedxbP);
  		ETmp -= simEloss(b,13./27.,ETmp,0.3*2.702*0.1/Cos(tr.T),ebAl, dedxbAl);
  		ETmp -= simEloss(b,100./192.,ETmp,6.*1.4*0.1/Cos(tr.T),ebMy, dedxbMy);
		csi.dE[mul] = simEloss(b,108./260.,ETmp,CsIThickness/Cos(tr.T),ebCsI, dedxbCsI);
		csi.dE[mul] = rndm->Gaus(csi.dE[mul],0.031*14.1*Sqrt(14.1/csi.dE[mul]));
	}
	if(shield && mask && Sd1Hit){
		mul = sd1.mul-1;
		ETmp = ETmp - simEloss(b,13./27.,ETmp,0.1*2.702*1.5/Cos(tr.T),ebAl,dedxbAl); //first metal
		ETmp = ETmp - simEloss(b,30./60.,ETmp,0.1*2.65*3.5/Cos(tr.T),ebSiO2,dedxbSiO2); //SiO2
		ETmp = ETmp - simEloss(b,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ebAl,dedxbAl); //second metal
		ETmp = ETmp - simEloss(b,5./10.,ETmp,0.1*2.3502*0.5/Cos(tr.T),ebB,dedxbB); //boron junction implant 		
		sd1.dE[mul] = simEloss(b,14./28.,ETmp,SdThickness1/Cos(tr.T),ebSi,dedxbSi);
   		ETmp = ETmp - sd1.dE[mul];
		sd1.dE[mul] = rndm->Gaus(sd1.dE[mul],0.01*sd1.dE[mul]);
	}
	if(shield && mask && Sd2Hit){
		mul = sd2.mul-1;
		ETmp = ETmp - simEloss(b,15./31.,ETmp,0.1*1.8219*0.5/Cos(tr.T),ebP,dedxbP); //phosphorus implant
		ETmp = ETmp - simEloss(b,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ebAl,dedxbAl); //metal
		ETmp = ETmp - simEloss(b,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ebAl,dedxbAl); //metal
		ETmp = ETmp - simEloss(b,15./31.,ETmp,0.1*1.822*0.5/Cos(tr.T),ebP,dedxbP); //phosphorus implant
		sd2.dE[mul] = simEloss(b,14./28.,ETmp,SdThickness2/Cos(tr.T),ebSi,dedxbSi);
		sd2.dE[mul] = rndm->Gaus(sd2.dE[mul],0.01*sd2.dE[mul]);
	}

	return (mask && shield && YYHit && CsIHit);
}

Bool_t detHitsH(PTrack tr, nucleus B, Double_t targetTh, TVector3 targetPos) // Calculate energy loss in detectors for heavy ejectile
{
	const Double_t SdThickness1=60.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t SdThickness2=1000.*2.3212*0.1; //mu*g/cm^3*0.1
	const Double_t YYThickness[8]={ 112.*2.3212*0.1, 109.*2.3212*0.1, 110.*2.3212*0.1, 106.*2.3212*0.1, 101.*2.3212*0.1, 109.*2.3212*0.1, 111.*2.3212*0.1, 103.*2.3212*0.1 };
	const Double_t CsIThickness=12000.*4.51*0.1; //mu*g/cm^3*0.1
	
	Bool_t mask = maskClear(tr.T,tr.P);
	Bool_t shield = shieldClear(tr.T,tr.P);
	Bool_t YYHit = yd.Calculate(tr.T,tr.P,prm.DYY,targetPos) ;
	Bool_t CsIHit = csi.Calculate(tr.T,tr.P,prm.DYY+11.6,targetPos) ;
	Bool_t Sd1Hit = sd1.Calculate(tr.T,tr.P,prm.DS3,targetPos);
	Bool_t Sd2Hit = sd2.Calculate(tr.T,tr.P,prm.DS3+14.8,targetPos);

	Double_t ETmp = 0.;

	TRandom3 *rndm = new TRandom3(0);
	Int_t mul;

	TargdE[1] = simEloss(B,1.,tr.E,targetTh/2./Cos(tr.T),eBH,dedxBH);
	ETmp = tr.E - TargdE[1];
	
	if(mask && shield && YYHit){
		mul = yd.mul-1;
  		ETmp -= simEloss(B,13./27.,ETmp,0.1*2.702*0.1/Cos(tr.T),eBAl, dedxBAl);
  		ETmp -= simEloss(B,5./10.,ETmp,0.05*2.3502*0.1/Cos(tr.T),eBB, dedxBB);
  		yd.dE[mul] = simEloss(B,14./28.,ETmp,YYThickness[yd.Seg[mul]]/Cos(tr.T),eBSi, dedxBSi);
		ETmp = ETmp-yd.dE[mul];
  		yd.dE[mul] = rndm->Gaus(yd.dE[mul],0.00225*5.73*Sqrt(5.73/yd.dE[mul]));
	}
	if(mask && shield && CsIHit){
		mul = csi.mul-1;
  		ETmp -= simEloss(B,15./31.,ETmp,0.1*1.8219*0.1/Cos(tr.T),eBP, dedxBP);
  		ETmp -= simEloss(B,13./27.,ETmp,0.3*2.702*0.1/Cos(tr.T),eBAl, dedxBAl);
  		ETmp -= simEloss(B,100./192.,ETmp,6.*1.4*0.1/Cos(tr.T),eBMy, dedxBMy);
		csi.dE[mul] = simEloss(B,108./260.,ETmp,CsIThickness/Cos(tr.T),eBCsI, dedxBCsI);
		csi.dE[mul] = rndm->Gaus(csi.dE[mul],0.031*14.1*Sqrt(14.1/csi.dE[mul]));
	}
	if(shield && mask && Sd1Hit){
		mul = sd1.mul-1;
		ETmp = ETmp - simEloss(B,13./27.,ETmp,0.1*2.702*1.5/Cos(tr.T),eBAl,dedxBAl); //first metal
		ETmp = ETmp - simEloss(B,30./60.,ETmp,0.1*2.65*3.5/Cos(tr.T),eBSiO2,dedxBSiO2); //SiO2
		ETmp = ETmp - simEloss(B,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),eBAl,dedxBAl); //second metal
		ETmp = ETmp - simEloss(B,5./10.,ETmp,0.1*2.3502*0.5/Cos(tr.T),eBB,dedxBB); //boron junction implant 		
		printf("%s ",B.name.data());
		sd1.dE[mul] = simEloss(B,14./28.,ETmp,SdThickness1/Cos(tr.T),eBSi,dedxBSi);
   		ETmp = ETmp - sd1.dE[mul];
		sd1.dE[mul] = rndm->Gaus(sd1.dE[mul],0.01*sd1.dE[mul]);
	}
	if(shield && mask && Sd2Hit){
		mul = sd2.mul-1;
		ETmp = ETmp - simEloss(B,15./31.,ETmp,0.1*1.8219*0.5/Cos(tr.T),eBP,dedxBP); //phosphorus implant
		ETmp = ETmp - simEloss(B,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),eBAl,dedxBAl); //metal
		ETmp = ETmp - simEloss(B,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),eBAl,dedxBAl); //metal
		ETmp = ETmp - simEloss(B,15./31.,ETmp,0.1*1.822*0.5/Cos(tr.T),eBP,dedxBP); //phosphorus implant
		sd2.dE[mul] = simEloss(B,14./28.,ETmp,SdThickness2/Cos(tr.T),eBSi,dedxBSi);
		sd2.dE[mul] = rndm->Gaus(sd2.dE[mul],0.01*sd2.dE[mul]);
	}

	return (shield && mask && Sd1Hit && Sd2Hit);
}

