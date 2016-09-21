#include "header.h"

Bool_t detHits(PTrack tr, nucleus ncl, Double_t targetTh, TVector3 targetPos) // Calculate energy loss in detectors for light ejectile
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

	Double_t ETmp = tr.Ebt;

	TRandom3 *rndm = new TRandom3(0);
	Int_t mul;

	if(mask && shield && YYHit){
		mul = yd.mul-1;
  		ETmp -= eloss(ncl,13./27.,ETmp,0.1*2.702*0.1/Cos(tr.T),ncl.EL.eAl, ncl.EL.dedxAl);
  		ETmp -= eloss(ncl,5./10.,ETmp,0.05*2.3502*0.1/Cos(tr.T),ncl.EL.eB, ncl.EL.dedxB);
  		yd.dE[mul] = eloss(ncl,14./28.,ETmp,YYThickness[yd.Seg[mul]]/Cos(tr.T),ncl.EL.eSi, ncl.EL.dedxSi);
		ETmp = ETmp-yd.dE[mul];
  		yd.dE[mul] = rndm->Gaus(yd.dE[mul],0.00225*5.73*Sqrt(5.73/yd.dE[mul]));
	}
	if(mask && shield && CsIHit){
		mul = csi.mul-1;
  		ETmp -= eloss(ncl,15./31.,ETmp,0.1*1.8219*0.1/Cos(tr.T),ncl.EL.eP, ncl.EL.dedxP);
  		ETmp -= eloss(ncl,13./27.,ETmp,0.3*2.702*0.1/Cos(tr.T),ncl.EL.eAl, ncl.EL.dedxAl);
  		ETmp -= eloss(ncl,100./192.,ETmp,6.*1.4*0.1/Cos(tr.T),ncl.EL.eMy, ncl.EL.dedxMy);
		csi.dE[mul] = eloss(ncl,108./260.,ETmp,CsIThickness/Cos(tr.T),ncl.EL.eCsI, ncl.EL.dedxCsI);
		csi.dE[mul] = rndm->Gaus(csi.dE[mul],0.031*14.1*Sqrt(14.1/csi.dE[mul]));
	}
	if(shield && mask && Sd1Hit){
		mul = sd1.mul-1;
		ETmp = ETmp - eloss(ncl,13./27.,ETmp,0.1*2.702*1.5/Cos(tr.T),ncl.EL.eAl,ncl.EL.dedxAl); //first metal
		ETmp = ETmp - eloss(ncl,30./60.,ETmp,0.1*2.65*3.5/Cos(tr.T),ncl.EL.eSiO2,ncl.EL.dedxSiO2); //SiO2
		ETmp = ETmp - eloss(ncl,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ncl.EL.eAl,ncl.EL.dedxAl); //second metal
		ETmp = ETmp - eloss(ncl,5./10.,ETmp,0.1*2.3502*0.5/Cos(tr.T),ncl.EL.eB,ncl.EL.dedxB); //boron junction implant 		
		sd1.dE[mul] = eloss(ncl,14./28.,ETmp,SdThickness1/Cos(tr.T),ncl.EL.eSi,ncl.EL.dedxSi);
   		ETmp = ETmp - sd1.dE[mul];
		sd1.dE[mul] = rndm->Gaus(sd1.dE[mul],0.01*sd1.dE[mul]);
	}
	if(shield && mask && Sd2Hit){
		mul = sd2.mul-1;
		ETmp = ETmp - eloss(ncl,15./31.,ETmp,0.1*1.8219*0.5/Cos(tr.T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
		ETmp = ETmp - eloss(ncl,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
		ETmp = ETmp - eloss(ncl,13./27.,ETmp,0.1*2.702*0.3/Cos(tr.T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
		ETmp = ETmp - eloss(ncl,15./31.,ETmp,0.1*1.822*0.5/Cos(tr.T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
		sd2.dE[mul] = eloss(ncl,14./28.,ETmp,SdThickness2/Cos(tr.T),ncl.EL.eSi,ncl.EL.dedxSi);
		sd2.dE[mul] = rndm->Gaus(sd2.dE[mul],0.01*sd2.dE[mul]);
	}

	return (mask && shield && YYHit && CsIHit);
}

