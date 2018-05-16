#include <TMath.h>
#include "S3Hit.h"
#include "TRandom3.h"
#include "TVector3.h"

S3Hit::S3Hit()
{
	Orientation = 0; // 0 = rings first, 1 = sectors first
	Thickness = 60.*2.3212*0.1; //um*g/cm^3*0.1
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear;
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
}

void S3Hit::Init(Bool_t o, Double_t th)
{
	Orientation = o; // 0 = rings first, 1 = sectors first
	Thickness = th*2.3212*0.1; //um*g/cm^3*0.1
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear;
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
}

void S3Hit::Clear()
{
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fPhiRand.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear;
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
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
	Double_t fX0, fY0, fZ0;
  	Double_t fPhiCalc0, fPhiRand0;
  	Double_t fThetaCalc0, fThetaRand0;
  	Int_t Seg0, Ring0;

	//ring number in S3 and geometric efficiency
	fX0 = S3Distance*tan(theta)*cos(phi);
	fY0 = S3Distance*tan(theta)*sin(phi);
	
	TVector3 partVec(fX0,fY0,S3Distance);
	
	partVec = partVec + targetPos; //taking into account the beam position at the target
	
	fX0 = partVec.X();
	fY0 = partVec.Y();
	fZ0 = partVec.Z();
	theta = partVec.Theta();
	phi = partVec.Phi();

   	// geometric efficiency
	Bool_t hitBool = ((S3Distance*tan(theta)>RIn) && (S3Distance*tan(theta)<ROut));

	if (hitBool){
    	mul++;
	  	//hit[mul] = 1;	
		Ring0 = int(S3Distance*tan(theta) - RIn);
		Seg0 = int((TMath::Pi()+phi)*16./TMath::Pi());
 		fPhiCalc0 = -180.+(Seg0*11.25);
 		if (fPhiCalc0<-180.) fPhiCalc0 = fPhiCalc0+360.;
      	Double_t rndm = 0.99*fRandom.Rndm();
 		fPhiRand0 = -180.+(Seg0+rndm)*11.25;
 		if (fPhiRand0<-180.) fPhiRand0 = fPhiRand0+360.;
 		if(Orientation==1){
			Seg0 = 31-Seg0;
			// fPhiCalc[mul] = -fPhiCalc[mul];
			// fPhiRand[mul] = -fPhiRand[mul];
		}
		fThetaCalc0 = TMath::RadToDeg()*atan((RIn+Ring0+0.5)/S3Distance);
 		fThetaCalc0 = (fThetaCalc0>0) ? fThetaCalc0 : fThetaCalc0+180.;
      	rndm = 0.99*fRandom.Rndm();
 		fThetaRand0 = TMath::RadToDeg()*atan((RIn+Ring0+rndm)/S3Distance);
 		fThetaRand0 = (fThetaRand0>0) ? fThetaRand0 : fThetaRand0+180.;
		fX.push_back(fX0);
		fY.push_back(fY0);
		fZ.push_back(fZ0);
		fPhiCalc.push_back(fPhiCalc0);
		fPhiRand.push_back(fPhiRand0);
		fThetaCalc.push_back(fThetaCalc0);
		fThetaRand.push_back(fThetaRand0);
		//hit.push_back;
		Seg.push_back(Seg0);
		Ring.push_back(Ring0);
	}
//  	else{
//	  	hit[mul] = 0;	
//		fX[mul] = sqrt(-1);
//		fY[mul] = sqrt(-1);
//  	   	Ring[mul] = sqrt(-1);
//		Seg[mul] = sqrt(-1);
//		fPhiCalc[mul] = sqrt(-1);
//		fThetaCalc[mul] = sqrt(-1);
//	}
 	
	return hitBool;
}

Double_t S3Hit::ELoss(nucleus ncl, Double_t E, Double_t Theta)
{
  	Double_t dE0, dE_ideal0;
	Double_t T = (Theta<TMath::Pi()/2.) ? Theta : TMath::Pi()-Theta;
	TRandom3 *rndm = new TRandom3(0);
	
	if(Orientation==0){ // rings first
		E -= eloss(ncl,13./27.,E,0.1*2.702*1.5/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //first metal
		E -= eloss(ncl,30./60.,E,0.1*2.65*3.5/cos(T),ncl.EL.eSiO2,ncl.EL.dedxSiO2); //SiO2
		E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //second metal
		E -= eloss(ncl,5./10.,E,0.1*2.3502*0.5/cos(T),ncl.EL.eB,ncl.EL.dedxB); //boron junction implant 		
		dE0 = eloss(ncl,14./28.,E,Thickness/cos(T),ncl.EL.eSi,ncl.EL.dedxSi);
   		// E -= dE0;
		// if(dE0<0.) dE0 = -dE0;
		// dE_ideal0 = dE0;
		// if(dE0!=0.) dE0 = rndm->Gaus(dE0,0.01*dE0);
		// if(dE0<0.) dE0 = -dE0;
	}
	else{ // sectors first
		E -= eloss(ncl,15./31.,E,0.1*1.8219*0.5/cos(T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
		E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
		E -= eloss(ncl,13./27.,E,0.1*2.702*0.3/cos(T),ncl.EL.eAl,ncl.EL.dedxAl); //metal
		E -= eloss(ncl,15./31.,E,0.1*1.822*0.5/cos(T),ncl.EL.eP,ncl.EL.dedxP); //phosphorus implant
		dE0 = eloss(ncl,14./28.,E,Thickness/cos(T),ncl.EL.eSi,ncl.EL.dedxSi);
   		// E -= dE0;
		// if(dE0<0.) dE0 = -dE0;
		// dE_ideal0 = dE0;
		// if(dE0!=0.) dE0 = rndm->Gaus(dE0,0.01*dE0);
		// if(dE0<0.) dE0 = -dE0;
	}
	E -= dE0;
	if(dE0<0.) dE0 = -dE0;
	dE_ideal0 = dE0;
	if(dE0!=0.) dE0 = rndm->Gaus(dE0,0.0225*dE0);
	if(dE0<0.) dE0 = -dE0;
	dE.push_back(dE0);
	dE_ideal.push_back(dE_ideal0);
	rndm->Delete();
	return E;
}

void S3Hit::SortByEnergy()
{
	Bool_t have_swapped = true;
	while(have_swapped == true){
		for (Int_t x=0; x<mul; x++)
		{
			have_swapped = false;
			for(Int_t y=0; y<mul-1; y++){
				if(dE[y]<dE[y+1]){
					std::swap(dE[y],dE[y+1]);
					std::swap(dE_ideal[y],dE_ideal[y+1]);
					std::swap(fX[y],fX[y+1]);
					std::swap(fY[y],fY[y+1]);
					std::swap(fZ[y],fZ[y+1]);
					std::swap(fPhiCalc[y],fPhiCalc[y+1]);
					std::swap(fPhiRand[y],fPhiRand[y+1]);
					std::swap(fThetaCalc[y],fThetaCalc[y+1]);
					std::swap(fThetaRand[y],fThetaRand[y+1]);
					std::swap(Seg[y],Seg[y+1]);
					std::swap(Ring[y],Ring[y+1]);
					have_swapped = true;
				}
			}
		}
	}
}
