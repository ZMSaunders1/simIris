#include <TMath.h>
#include "S3Hit.h"
#include "TRandom3.h"
#include "TVector3.h"

S3Hit::S3Hit()
{
	Orientation = 0; // 0 = rings first, 1 = sectors first
	Thickness = 60.*2.3212*0.1; //um*g/cm^3*0.1
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

void S3Hit::Init(Bool_t o, Double_t th)
{
	Orientation = o; // 0 = rings first, 1 = sectors first
	Thickness = th*2.3212*0.1; //um*g/cm^3*0.1
	mul = 0;
	for(UInt_t i=0; i<4; i++){
		fX[i] = sqrt(-1);
		fY[i] = sqrt(-1);
		fZ[i] = sqrt(-1);
		fPhiCalc[i] = sqrt(-1);
		fThetaCalc[i] = sqrt(-1);
		hit[i] = 0;
		Seg[i] = sqrt(-1);
		Ring[i] = sqrt(-1);
		dE[i] = sqrt(-1);
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
		Seg[i] = sqrt(-1);
		Ring[i] = sqrt(-1);
		dE[i] = sqrt(-1);
	}
}

Double_t S3Hit::ThetaMin(Double_t S3Distance)
{	
	const Double_t RIn = 11.;	// Inner radius in mm
	Double_t theta_min = TMath::ATan2(RIn,S3Distance)*TMath::RadToDeg();
	return theta_min;
}

Double_t S3Hit::ThetaMax(Double_t S3Distance)
{	
	const Double_t ROut = 35.;	// Outer radius in mm
	Double_t theta_max = TMath::ATan2(ROut,S3Distance)*TMath::RadToDeg();
	return theta_max;
}

Bool_t S3Hit::Hit(Double_t theta, Double_t phi, Double_t S3Distance, TVector3 targetPos)
{
	const Double_t RIn = 11.;	// Inner radius in mm
	const Double_t ROut = 35.;	// Outer radius in mm

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
	Bool_t hitBool = ((S3Distance*tan(theta)>RIn) && (S3Distance*tan(theta)<ROut));

	if (hitBool){
	  	hit[mul] = 1;	
		fX[mul] = fX0;
		fY[mul] = fY0;
		Ring[mul] = int(S3Distance*tan(theta) - RIn);
		Seg[mul] = int(phi/(TMath::Pi()/16.));
      	Double_t random = 0.5;//fRandom.Rndm();
  		Double_t randomTheta = fRandom.Rndm();
 		fPhiCalc[mul] = (Seg[mul]-0.5+random)*5.625;
 		if (fPhiCalc[mul]<-180.)	fPhiCalc[mul] = fPhiCalc[mul]+360.;
 		fThetaCalc[mul] = TMath::RadToDeg()*atan((RIn+Ring[mul]+randomTheta)/S3Distance);
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

Double_t S3Hit::ELoss(nucleus ncl, Double_t E, Double_t T)
{
	if(mul>0 && hit[mul-1]==1){
		TRandom3 *rndm = new TRandom3(0);
		if(Orientation==0){ // rings first
			E -= eloss(ncl,13./27.,E,0.1*2.702*1.5/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //first metal
			E -= eloss(ncl,30./60.,E,0.1*2.65*3.5/cos(T),ncl.EL.eSiO2,ncl.EL.dedxSiO2); //SiO2
			E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //second metal
			E -= eloss(ncl,5./10.,E,0.1*2.3502*0.5/cos(T),ncl.EL.eB,ncl.EL.dedxB); //boron junction implant 		
			dE[mul-1] = eloss(ncl,14./28.,E,Thickness/cos(T),ncl.EL.eSi,ncl.EL.dedxSi);
   			E -= dE[mul-1];
			if(dE[mul-1]<0.) dE[mul-1] = -dE[mul-1];
			dE_ideal[mul-1] = dE[mul-1];
			if(dE[mul-1]!=0.) dE[mul-1] = rndm->Gaus(dE[mul-1],0.01*dE[mul-1]);
			if(dE[mul-1]<0.) dE[mul-1] = -dE[mul-1];
		}
		else{ // sectors first
			E -= eloss(ncl,15./31.,E,0.1*1.8219*0.5/cos(T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
			E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
			E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
			E -= eloss(ncl,15./31.,E,0.1*1.822*0.5/cos(T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
			dE[mul-1] = eloss(ncl,14./28.,E,Thickness/cos(T),ncl.EL.eSi,ncl.EL.dedxSi);
   			E -= dE[mul-1];
			if(dE[mul-1]<0.) dE[mul-1] = -dE[mul-1];
			dE_ideal[mul-1] = dE[mul-1];
			if(dE[mul-1]!=0.) dE[mul-1] = rndm->Gaus(dE[mul-1],0.01*dE[mul-1]);
			if(dE[mul-1]<0.) dE[mul-1] = -dE[mul-1];
		}
	}
	return E;
}
