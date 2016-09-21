#include <TMath.h>
#include "YYHit.h"
#include "TRandom3.h"
#include "TVector3.h"

YYHit::YYHit()
{
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		fThetaCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = -1;
		Ring[i] = -1;
		dE[i] = 0;
	}
}

void YYHit::Clear()
{
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		fThetaCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = -1;
		Ring[i] = -1;
		dE[i] = 0;
	}
}

Double_t YYHit::ThetaMin(Double_t YdDistance)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	Double_t theta_min = TMath::ATan2(RIn,YdDistance)*TMath::RadToDeg();
	return theta_min;
}

Double_t YYHit::ThetaMax(Double_t YdDistance)
{	
	const Double_t ROut = 130.;	// Outer radius in mm
	Double_t theta_max = TMath::ATan2(ROut,YdDistance)*TMath::RadToDeg();
	return theta_max;
}

Bool_t YYHit::Calculate(Double_t theta, Double_t phi, Double_t YdDistance, TVector3 targetPos)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	const Double_t ROut = 130.;	// Outer radius in mm

	TRandom3 fRandom(0);
	Double_t fX0, fY0, Ring0, Seg0;
	
	Bool_t hitTheta=0;//does it hit in theta range?
	Bool_t hitPhi=0;//does it hit in phi range?
	Bool_t hitBool=0;//return value
	Double_t phiGap = 4.25*TMath::DegToRad(); //phi gap between YY1s in rad 
	Double_t phiShift =  9.5*TMath::DegToRad();//shift from the vertical direction for the first YY1
	Double_t phiRel;//Relative phi after phishift
	Double_t phiRange = TMath::Pi()/4-phiGap; //ring dependent phi range for each YY1 

	phiRel = phi+phiShift;
	if (phiRel>TMath::Pi())  phiRel = phiRel - 2*TMath::Pi();
	//ring number in YY1 and geometric efficiency
	fX0 = YdDistance*tan(theta)*cos(phi);
	fY0 = YdDistance*tan(theta)*sin(phi);
	
	TVector3 partVec(fX0,fY0,YdDistance);
	partVec = partVec + targetPos; //taking into account the beam position at the target
	
	fX0 = partVec.X();
	fY0 = partVec.Y();
	theta = partVec.Theta();
	phi = partVec.Phi();
  	
  	Ring0 = int((YdDistance*tan(theta) - RIn)/5.);
  	Seg0 = int((phiRel+TMath::Pi())/(TMath::Pi()/4));
	//if (Ring<13 && Ring>=0)  phiRange = TMath::Pi()/4-phiGap;
	if (Ring0==13) phiRange =  TMath::Pi()/4-phiGap*2; 
	else if (Ring0==14) phiRange =  TMath::Pi()/4-phiGap*3.5; 
	else if (Ring0==15) phiRange =  TMath::Pi()/4-phiGap*5.5; 
	hitPhi = fabs(phiRel +TMath::Pi() - Seg0*TMath::Pi()/4-TMath::Pi()/8) < phiRange/2;
	hitTheta = ((YdDistance*tan(theta)>RIn) && (YdDistance*tan(theta)<ROut));
	hitBool = (hitTheta && hitPhi);

	if (hitBool){
		hit[mul] = 1;
	  	fX[mul] = fX0;
		fY[mul] = fY0;	
  		Seg[mul] = Seg0;
  		Ring[mul] = Ring0;
      	Double_t random = 0.5;//fRandom.Rndm();
  		Double_t randomTheta = fRandom.Rndm();
 		fPhiCalc[mul] = (Seg[mul]-0.5+random)*45.+phiShift*TMath::RadToDeg()-180.;
 		if (fPhiCalc[mul]<-180.)	fPhiCalc[mul] = fPhiCalc[mul]+360.;
 		fThetaCalc[mul] = TMath::RadToDeg()*atan((50.+(Ring[mul]*5.)+5.*randomTheta)/YdDistance);
   		mul++;
   	}
  	else{
		hit[mul] = 0;
		fX[mul] = sqrt(-1);
		fY[mul] = sqrt(-1);
  	   	Ring[mul] = -1;
  	   	Seg[mul] = -1;
		fPhiCalc[mul] = sqrt(-1);
		fThetaCalc[mul] = sqrt(-1);
	}
 	
	return hitBool;
}
