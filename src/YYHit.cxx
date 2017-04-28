#include <TMath.h>
#include "YYHit.h"
#include "TRandom3.h"
#include "TVector3.h"

YYHit::YYHit()
{
	Thickness[0] = 104.65 * 2.3212 * 0.1; 
	Thickness[1] = 101.15 * 2.3212 * 0.1; 
	Thickness[2] = 106.13 * 2.3212 * 0.1; 
	Thickness[3] = 101.75 * 2.3212 * 0.1; 
	Thickness[4] = 100.05 * 2.3212 * 0.1; 
	Thickness[5] = 105.65 * 2.3212 * 0.1; 
	Thickness[6] = 102.48 * 2.3212 * 0.1; 
	Thickness[7] = 105.84 * 2.3212 * 0.1; 
	Avg_Thickness = 103.46 * 2.3212 * 0.1; 
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear();
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
}

void YYHit::Init(Double_t th[8])
{
	Thickness[0] = th[0] * 2.3212 * 0.1; 
	Thickness[1] = th[1] * 2.3212 * 0.1; 
	Thickness[2] = th[2] * 2.3212 * 0.1; 
	Thickness[3] = th[3] * 2.3212 * 0.1; 
	Thickness[4] = th[4] * 2.3212 * 0.1; 
	Thickness[5] = th[5] * 2.3212 * 0.1; 
	Thickness[6] = th[6] * 2.3212 * 0.1; 
	Thickness[7] = th[7] * 2.3212 * 0.1; 
	Avg_Thickness = (th[0]+th[1]+th[2]+th[3]+th[4]+th[5]+th[6]+th[7])/8. * 2.3212 * 0.1; 
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear();
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
}

void YYHit::Clear()
{
	mul = 0;
	fX.clear();
	fY.clear();
	fZ.clear();
	fPhiCalc.clear();
	fThetaCalc.clear();
	fThetaRand.clear();
	//hit.clear();
	Seg.clear();
	Ring.clear();
	dE.clear();
	dE_ideal.clear();
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

Bool_t YYHit::Hit(Double_t theta, Double_t phi, Double_t YdDistance, TVector3 targetPos)
{	
	const Double_t RIn = 50.;	// Inner radius in mm
	const Double_t ROut = 130.;	// Outer radius in mm

	TRandom3 fRandom(0);
	Double_t fX0, fY0, fZ0;
	Double_t fPhiCalc0, fThetaCalc0, fThetaRand0;
   	Int_t Seg0, Ring0;
	
	Bool_t hitTheta=0;//does it hit in theta range?
	Bool_t hitPhi=0;//does it hit in phi range?
	Bool_t hitBool=0;//return value
	Double_t phiGap = 4.25*TMath::DegToRad(); //phi gap between YY1s in rad 
	Double_t phiShift =  -9.5*TMath::DegToRad();//shift from the vertical direction for the first YY1
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
	fZ0 = partVec.Z();
	theta = partVec.Theta();
	phi = partVec.Phi();
	
	Seg0 = int((phiRel+TMath::Pi())/(TMath::Pi()/4));
  	Ring0 = int((YdDistance*tan(theta) - RIn)/5.);
		//if (Ring<13 && Ring>=0)  phiRange = TMath::Pi()/4-phiGap;
	if (Ring0==13) phiRange =  TMath::Pi()/4-phiGap*2; 
	else if (Ring0==14) phiRange =  TMath::Pi()/4-phiGap*3.5; 
	else if (Ring0==15) phiRange =  TMath::Pi()/4-phiGap*5.5; 
	hitPhi = fabs(phiRel +TMath::Pi() - Seg0*TMath::Pi()/4-TMath::Pi()/8) < phiRange/2;
	hitTheta = ((YdDistance*tan(theta)>RIn) && (YdDistance*tan(theta)<ROut));
	hitBool = (hitTheta && hitPhi);
	
	if (hitBool){
   		mul++;
		//hit0 = 1;
	  	fX0 = fX0;
		fY0 = fY0;	
  		Seg0 = int((TMath::Pi()+22.5*TMath::DegToRad()-phiRel)*4./TMath::Pi());
  		Ring0 = Ring0;
  		Double_t randomTheta = fRandom.Uniform();
 		fPhiCalc0 = 180.-(Seg0*45.+phiShift*TMath::RadToDeg());
 		if (fPhiCalc0<-180.)	fPhiCalc0 = fPhiCalc0+360.;
		Seg0 = Seg0 +6;
		if(Seg0>7) Seg0 = Seg0 - 8;
 		fThetaCalc0 = TMath::RadToDeg()*atan((50.+(Ring0*5.)+2.5)/YdDistance);
 		fThetaRand0 = TMath::RadToDeg()*atan((50.+(Ring0*5.)+5.*randomTheta)/YdDistance);
   		fX.push_back(fX0);
		fY.push_back(fY0);
		fZ.push_back(fZ0);
		fPhiCalc.push_back(fPhiCalc0);
		fThetaCalc.push_back(fThetaCalc0);
		fThetaRand.push_back(fThetaRand0);
		//hit.push_back();
		Seg.push_back(Seg0);
		Ring.push_back(Ring0);
	}
  	//else{
	//	hit[mul] = 0;
	//	fX[mul] = sqrt(-1);
	//	fY[mul] = sqrt(-1);
  	//   	Ring[mul] = sqrt(-1);
  	//   	Seg[mul] = sqrt(-1);
	//	fPhiCalc[mul] = sqrt(-1);
	//	fThetaCalc[mul] = sqrt(-1);
	//	fThetaRand[mul] = sqrt(-1);
	//}
 	
	return hitBool;
}

Double_t YYHit::ELoss(nucleus ncl, Double_t E, Double_t T)
{
	//if(mul>0 && hit[mul-1]==1) {	
  		Double_t dE0, dE_ideal0;
		TRandom3 *rndm = new TRandom3(0);
		E -= eloss(ncl,13./27.,E,0.1*2.702*0.1/cos(T),ncl.EL.eAl, ncl.EL.dedxAl);
  		E -= eloss(ncl,5./10.,E,0.05*2.3502*0.1/cos(T),ncl.EL.eB, ncl.EL.dedxB);
  		dE0 = eloss(ncl,14./28.,E,Thickness[Seg.at(Seg.size()-1)]/cos(T),ncl.EL.eSi, ncl.EL.dedxSi);
  		dE_ideal0 = eloss(ncl,14./28.,E,Avg_Thickness/cos(T),ncl.EL.eSi, ncl.EL.dedxSi);
		E = E-dE0;
		if(dE0<0.) dE0 = -dE0;
  		if(dE0!=0.) dE0 = rndm->Gaus(dE0,0.00225*dE0*sqrt(5.73/dE0));
		if(dE0<0.) dE0 = 0.;
		dE.push_back(dE0);
		dE_ideal.push_back(dE_ideal0);
//}
	return E;
}

void YYHit::SortByEnergy()
{
	Bool_t have_swapped = true;
	while(have_swapped == true){
		for (Int_t x=0; x<mul; x++){
			have_swapped = false;
			for(Int_t y=0; y<mul-1; y++){
				if(dE[y]<dE[y+1]){
					std::swap(dE[y],dE[y+1]);
					std::swap(dE_ideal[y],dE_ideal[y+1]);
					std::swap(fX[y],fX[y+1]);
					std::swap(fY[y],fY[y+1]);
					std::swap(fZ[y],fZ[y+1]);
					std::swap(fPhiCalc[y],fPhiCalc[y+1]);
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
