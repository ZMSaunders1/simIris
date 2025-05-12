#include "header.h"

Bool_t detHits(PTrack tr, nucleus ncl, TVector3 reacPos, Bool_t maskIn, Bool_t shieldIn)
{
	Bool_t mask = maskClear(tr.T,tr.P) || !maskIn;
	Bool_t shield = shieldClear(tr.T,tr.P) || !shieldIn;
	Bool_t forward = (tr.T<TMath::Pi()/2.);
	Bool_t backward = (tr.T>TMath::Pi()/2.);
	Bool_t YYHit = 0;
	Bool_t CsIHit = 0;
	Bool_t Sd1Hit = 0;
	Bool_t Sd2Hit = 0;
	Bool_t YuHit = 0;
	Bool_t SuHit = 0;

	Double_t ETmp = tr.Ebt;
	Double_t ETmpU = tr.Ebt;

	if(backward){
		YuHit = yu.Hit(tr.T,tr.P,geoPrm.DYYU,reacPos) ;
		SuHit = su.Hit(tr.T,tr.P,geoPrm.DS3U,reacPos);
		if(YuHit) ETmpU = yu.ELoss(ncl,ETmpU,tr.T);
		if(SuHit) ETmpU = su.ELoss(ncl,ETmpU,tr.T);
	}	
	
	if(mask && shield&&forward){
		YYHit = yd.Hit(tr.T,tr.P,geoPrm.DYY,reacPos) ;
		if(YYHit) CsIHit = csi.Hit(tr.T,tr.P,geoPrm.DYY+11.6,reacPos) ;
		Sd1Hit = sd1.Hit(tr.T,tr.P,geoPrm.DS3,reacPos);
		if(Sd1Hit) Sd2Hit = sd2.Hit(tr.T,tr.P,geoPrm.DS3+14.8,reacPos);
		if(YYHit) ETmp = yd.ELoss(ncl,ETmp,tr.T);
		if(CsIHit) ETmp = csi.ELoss(ncl,ETmp,tr.T);
		if(Sd1Hit) ETmp = sd1.ELoss(ncl,ETmp,tr.T);
		if(Sd2Hit) ETmp = sd2.ELoss(ncl,ETmp,tr.T);
	}	
	
	return (mask && shield && YYHit && CsIHit);
}

// Calculate the energy loss of the scattered particles in Foil and SHT
PTrack TgtELoss(PTrack tr, nucleus ncl, geoParams g, Double_t reacZ, Bool_t isSHTReac)
{
	if(isSHTReac){ //Reaction in SHT
		if(g.Orientation==0&&tr.T<TMath::Pi()/2.){ // foil before target, theta<90 deg
	   		tr.FoildE = 0.;	
			tr.TrgtdE = eloss(ncl,1./g.AoZTgt,tr.E,(g.TTgt-reacZ)/TMath::Cos(tr.T),ncl.EL.eTgt,ncl.EL.dedxTgt);	
		}
		if(g.Orientation==0&&tr.T>TMath::Pi()/2.){ // foil before target, theta>90 deg
			tr.TrgtdE = eloss(ncl,1./g.AoZTgt,tr.E,reacZ/TMath::Cos(TMath::Pi()-tr.T),ncl.EL.eTgt,ncl.EL.dedxTgt);	
			tr.FoildE = eloss(ncl,1./g.AoZFoil,tr.E-tr.TrgtdE,g.TFoil/TMath::Cos(TMath::Pi()-tr.T),ncl.EL.eFoil,ncl.EL.dedxFoil);	
		}
		if(g.Orientation==1&&tr.T<TMath::Pi()/2.){ // foil after target, theta<90 deg
			tr.TrgtdE = eloss(ncl,1./g.AoZTgt,tr.E,(g.TTgt-reacZ)/TMath::Cos(tr.T),ncl.EL.eTgt,ncl.EL.dedxTgt);	
			tr.FoildE = eloss(ncl,1./g.AoZFoil,tr.E-tr.TrgtdE,g.TFoil/TMath::Cos(tr.T),ncl.EL.eFoil,ncl.EL.dedxFoil);	
		}
		if(g.Orientation==1&&tr.T>TMath::Pi()/2.){ // foil after target, theta>90 deg
			tr.TrgtdE = eloss(ncl,1./g.AoZTgt,tr.E,reacZ/TMath::Cos(TMath::Pi()-tr.T),ncl.EL.eTgt,ncl.EL.dedxTgt);	
	   		tr.FoildE = 0.;	
		}
	}
    else{ //Reaction in foil
        
        if(tr.T<TMath::Pi()/2.){
            tr.TrgtdE = 0.;
            tr.FoildE = eloss(ncl,1./g.AoZFoil,tr.E,(g.TFoil-reacZ)/TMath::Cos(tr.T),ncl.EL.eFoil,ncl.EL.dedxFoil);
        }
        if(tr.T>TMath::Pi()/2.){ // foil after target, theta>90 deg
            tr.FoildE = eloss(ncl,1./g.AoZFoil,tr.E,reacZ/TMath::Cos(TMath::Pi()-tr.T),ncl.EL.eFoil,ncl.EL.dedxFoil);
            tr.TrgtdE =0;
        }
    }
	
	tr.Ebt = tr.E-tr.FoildE-tr.TrgtdE; // calculate energy of particle after foil and target
	//printf("In: %f\tFoil: %f\tTarget: %f\tLeft: %f\n",tr.E,tr.FoildE,tr.TrgtdE,tr.Ebt);
	return tr;
}

void sortEnergies(){
	if(yd.mul>1) yd.SortByEnergy();
	if(yu.mul>1) yu.SortByEnergy();
	if(csi.mul>1) csi.SortByEnergy();
	if(sd1.mul>1) sd1.SortByEnergy();
	if(sd2.mul>1) sd2.SortByEnergy();
	if(su.mul>1) su.SortByEnergy();
}

void setIDet(Double_t ICdE, Double_t SSBdE)
{
	if(yd.mul>0)
	{
  		det.TYdMul = yd.mul;
		for(Int_t i=0; i<det.TYdMul; i++){
			det.TYdEnergy.push_back(yd.dE[i]);
  			det.TYdTheta.push_back(yd.fThetaRand[i]);// Yd theta angle                                                                       
			det.TYdChannel.push_back(yd.Seg[i]*16+yd.Ring[i]);
			det.TYdNo.push_back(yd.Seg[i]);
			det.TYdRing.push_back(yd.Ring[i]);
		}
	}

	if(yu.mul>0)
	{
  		det.TYuMul = yu.mul;
		for(Int_t i=0; i<det.TYuMul; i++){
			det.TYuEnergy.push_back(yu.dE[i]);
  			det.TYuTheta.push_back(yu.fThetaRand[i]);// Yu theta angle                                                                       
			det.TYuChannel.push_back(yu.Seg[i]*16+yu.Ring[i]);
			det.TYuNo.push_back(yu.Seg[i]);
			det.TYuRing.push_back(yu.Ring[i]);
		}
	}

	if(csi.mul>0)
	{
		det.TCsI1Mul = csi.mul;
		for(Int_t i=0; i<det.TCsI1Mul; i++){
  			det.TCsI1Energy.push_back(csi.dE[i]);
			det.TCsI1Channel.push_back(csi.Seg[i]);
			det.TCsI1Phi.push_back(csi.fPhiRand[i]);
		}
		det.TCsI2Mul = csi.mul;
		for(Int_t i=0; i<det.TCsI2Mul; i++){
  			det.TCsI2Energy.push_back(csi.dE[i]);
			det.TCsI2Channel.push_back(csi.Seg[i]);
			det.TCsI2Phi.push_back(csi.fPhiRand[i]);
		}
	}

	det.TSSBEnergy=SSBdE;
   	det.TICEnergy.push_back(ICdE);
	det.TICChannel.push_back(15);

	if(sd1.mul>0)
	{
		det.TSd1rMul=sd1.mul;
		for(Int_t i=0; i<det.TSd1rMul; i++){
  			det.TSd1rEnergy.push_back(sd1.dE[i]);
			det.TSd1rChannel.push_back(sd1.Ring[i]);
  			det.TSd1Theta.push_back(sd1.fThetaRand[i]);
		}
		det.TSd1sMul=sd1.mul;
		for(Int_t i=0; i<det.TSd1sMul; i++){
  			det.TSd1sEnergy.push_back(sd1.dE[i]);
			det.TSd1sChannel.push_back(sd1.Seg[i]);
  			det.TSd1Phi.push_back(sd1.fPhiRand[i]);
		}
	}
	if(sd2.mul>0)
	{
		det.TSd2rMul=sd2.mul;
		for(Int_t i=0; i<det.TSd1rMul; i++){
  			det.TSd2rEnergy.push_back(sd2.dE[i]);
			det.TSd2rChannel.push_back(sd2.Ring[i]);
  			det.TSd2Theta.push_back(sd2.fThetaRand[i]);
		}
		det.TSd2sMul=sd2.mul;
		for(Int_t i=0; i<det.TSd1sMul; i++){
  			det.TSd2sEnergy.push_back(sd2.dE[i]);
			det.TSd2sChannel.push_back(sd2.Seg[i]);
  			det.TSd2Phi.push_back(sd2.fPhiRand[i]);
		}
	}
	
	if(su.mul>0)
	{
		det.TSurMul=su.mul;
		for(Int_t i=0; i<det.TSurMul; i++){
  			det.TSurEnergy.push_back(su.dE[i]);
			det.TSurChannel.push_back(su.Ring[i]);
  			det.TSuTheta.push_back(su.fThetaRand[i]);
		}
		det.TSusMul=su.mul;
		for(Int_t i=0; i<det.TSusMul; i++){
  			det.TSusEnergy.push_back(su.dE[i]);
			det.TSusChannel.push_back(su.Seg[i]);
  			det.TSuPhi.push_back(su.fPhiRand[i]);
		}
	}
}
