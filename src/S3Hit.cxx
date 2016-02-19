#include <TMath.h>
#include "S3Hit.h"
#include "TRandom3.h"
#include "TVector3.h"

S3Hit::S3Hit()
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

void S3Hit::Clear()
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

Bool_t S3Hit::Calculate(Double_t theta, Double_t phi, Double_t S3Distance, TVector3 targetPos)
{
	TRandom3 fRandom(0);
	Double_t fX0, fY0;
	
	//ring number in S3 and geometric efficiency
	fX0 = S3Distance*tan(theta)*cos(phi);
	fY0 = S3Distance*tan(theta)*sin(phi);
	
	TVector3 partVec(fX0,fY0,S3Distance);
	
	partVec = partVec + targetPos; //taking into account the beam position at the target
	
	fX0 = partVec.X();
	fY0 = partVec.Y();
	theta = partVec.Theta();
	phi = partVec.Phi();

   	// geometric efficiency
	Bool_t hitBool = ((S3Distance*tan(theta)>11.) && (S3Distance*tan(theta)<35.));

	if (hitBool){
	  	hit[mul] = 1;	
		fX[mul] = fX0;
		fY[mul] = fY0;
		Ring[mul] = int(S3Distance*tan(theta) - 11.);
		Seg[mul] = int(phi/(TMath::Pi()/16.));
      	Double_t random = 0.5;//fRandom.Rndm();
  		Double_t randomTheta = fRandom.Rndm();
 		fPhiCalc[mul] = (Seg[mul]-0.5+random)*5.625;
 		if (fPhiCalc[mul]<-180.)	fPhiCalc[mul] = fPhiCalc[mul]+360.;
 		fThetaCalc[mul] = TMath::RadToDeg()*atan((11.+Ring[mul]+randomTheta)/S3Distance);
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
