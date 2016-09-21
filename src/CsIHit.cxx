#include <TMath.h>
#include "CsIHit.h"
#include "TRandom3.h"
#include "TVector3.h"

CsIHit::CsIHit()
{
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

void CsIHit::Clear()
{
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

Double_t CsIHit::ThetaMin(Double_t CsIDistance)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	Double_t theta_min = TMath::ATan2(RIn,CsIDistance)*TMath::RadToDeg();
	return theta_min;
}

Double_t CsIHit::ThetaMax(Double_t CsIDistance)
{	
	const Double_t ROut = 130.;	// Outer radius in mm
	Double_t theta_max = TMath::ATan2(ROut,CsIDistance)*TMath::RadToDeg();
	return theta_max;
}

Bool_t CsIHit::Calculate(Double_t theta, Double_t phi, Double_t CsIDistance, TVector3 targetPos)
{
	const Double_t RIn = 50.;	// Inner radius in mm
	const Double_t ROut = 130.;	// Outer radius in mm

	TRandom3 fRandom(0);
	Double_t fX0, fY0;
	Bool_t hitTheta = 0;
	Bool_t hitPhi = 0;
	Bool_t hitBool = 0; // return value
	//Double_t phiGap = 4.25*TMath::DegToRad(); //phi gap between CsI1s in rad 
	//Double_t phiShift =  9.5*TMath::DegToRad();//shift from the vertical direction for the first CsI1
	Double_t phiShift =  0.;//shift from the vertical direction for the first CsI1
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
	//hitPhi = fabs(phiRel +TMath::Pi() - Seg*TMath::Pi()/4-TMath::Pi()/8) < phiRange/2;
  	hitBool = (hitPhi && hitTheta); 
	if (hitBool){
		hit[mul] = 1;
	  	fX[mul] = fX0;
		fY[mul] = fY0;	
  		Seg[mul] = int((phiRel+TMath::Pi())/(TMath::Pi()/8));
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
