#include "CsIHit.h"

CsIHit::CsIHit()
{
	Thickness =12000.*4.51*0.1; //um*g/cm^3*0.1;
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = -1;
		dE[i] = 0;
	}
}

CsIHit::CsIHit(Double_t th)
{
	Thickness = th*4.51*0.1;
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = sqrt(-1);
		dE[i] = sqrt(-1);
	}
}

void CsIHit::Clear()
{
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = sqrt(-1);
		dE[i] = sqrt(-1);
	}
}

Double_t CsIHit::ThetaMin(Double_t CsIDistance)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	Double_t theta_min = TMath::ATan2(RIn,CsIDistance)*TMath::RadToDeg();
	return theta_min;
}

Double_t CsIHit::ThetaMax(Double_t CsIDistance)
{	
	const Double_t ROut = 150.;	// Outer radius in mm
	Double_t theta_max = TMath::ATan2(ROut,CsIDistance)*TMath::RadToDeg();
	return theta_max;
}

Bool_t CsIHit::Hit(Double_t theta, Double_t phi, Double_t CsIDistance, TVector3 targetPos)
{
	const Double_t RIn = 50.;	// Inner radius in mm
	const Double_t ROut = 150.;	// Outer radius in mm

	TRandom3 fRandom(0);
	Double_t fX0, fY0;
	Double_t Seg0;
	Bool_t hitTheta = 0;
	Bool_t hitPhi = 0;
	Bool_t hitBool = 0; // return value
	//Double_t phiGap = 4.25*TMath::DegToRad(); //phi gap between CsI1s in rad 
	Double_t phiShift = -9.5*TMath::DegToRad();//shift from the vertical direction for the first CsI1
	//Double_t phiShift =  0.;//shift from the vertical direction for the first CsI1
	Double_t phiRel;//Relative phi after phishift
	//Double_t phiRange; //ring dependent phi range for each CsI1 

	phiRel = phi+phiShift;
	if (phiRel>TMath::Pi())  phiRel = phiRel - 2*TMath::Pi();
	
	// geometric efficiency
	fX0 = CsIDistance*tan(theta)*cos(phi);
	fY0 = CsIDistance*tan(theta)*sin(phi);
	
	TVector3 partVec(fX0,fY0,CsIDistance);
	
	partVec = partVec + targetPos; //taking into account the beam position at the target
	
	fX0 = partVec.X();
	fY0 = partVec.Y();
	theta = partVec.Theta();
	phi = partVec.Phi();
	
	hitTheta = ((CsIDistance*tan(theta)>RIn) && (CsIDistance*tan(theta)<ROut));
	hitPhi = kTRUE;
	Seg0 = int((phiRel+TMath::Pi())/(TMath::Pi()/8));
	//hitPhi = fabs(phiRel +TMath::Pi() - Seg*TMath::Pi()/4-TMath::Pi()/8) < phiRange/2;
  	hitBool = (hitPhi && hitTheta); 
	
	Seg0 = 11-Seg0;
	Seg0 = (Seg0<0) ? Seg0+16 : Seg0;
	
	if (hitBool){
		hit[mul] = 1;
	  	fX[mul] = fX0;
		fY[mul] = fY0;	
  		Seg[mul] = Seg0;
      	Double_t random = 0.5;//fRandom.Rndm();
 		fPhiCalc[mul] = (Seg[mul]-0.5+random)*22.5+phiShift*TMath::RadToDeg()-180.;
 		if (fPhiCalc[mul]<-180.)	fPhiCalc[mul] = fPhiCalc[mul]+360.;
		mul ++;
    }
  	else{
		hit[mul] = 0;
		fX[mul] = sqrt(-1);
		fY[mul] = sqrt(-1);
		fZ[mul] = sqrt(-1);
		fPhiCalc[mul] = sqrt(-1);
		Seg[mul] = -1;
	}
 	
	return hitBool;
}

Double_t CsIHit::ELoss(nucleus ncl, Double_t E, Double_t T)
{
	if(mul>0 && hit[mul-1]==1){
		TRandom3 *rndm = new TRandom3(0);
		E -= eloss(ncl,15./31.,E,0.1*1.8219*0.1/cos(T),ncl.EL.eP, ncl.EL.dedxP);
		E -= eloss(ncl,13./27.,E,0.3*2.702*0.1/cos(T),ncl.EL.eAl, ncl.EL.dedxAl);
		E -= eloss(ncl,100./192.,E,6.*1.4*0.1/cos(T),ncl.EL.eMy, ncl.EL.dedxMy);
		dE[mul-1] = eloss(ncl,108./260.,E,Thickness/cos(T),ncl.EL.eCsI, ncl.EL.dedxCsI);
		E -= dE[mul-1];
		dE[mul-1] = rndm->Gaus(dE[mul-1],0.031*14.1*sqrt(14.1/dE[mul-1]));
	}
	return E;
}
