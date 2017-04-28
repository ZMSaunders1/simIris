#include "header.h"

Bool_t detHits(PTrack tr, nucleus ncl, TVector3 reacPos)
{
	Bool_t mask = maskClear(tr.T,tr.P);
	Bool_t shield = shieldClear(tr.T,tr.P);
	Bool_t YYHit = 0;
	Bool_t CsIHit = 0;
	Bool_t Sd1Hit = 0;
	Bool_t Sd2Hit = 0;

	Double_t ETmp = tr.Ebt;

	if(mask && shield){
		YYHit = yd.Hit(tr.T,tr.P,geoPrm.DYY,reacPos) ;
		CsIHit = csi.Hit(tr.T,tr.P,geoPrm.DYY+11.6,reacPos) ;
		Sd1Hit = sd1.Hit(tr.T,tr.P,geoPrm.DS3,reacPos);
		Sd2Hit = sd2.Hit(tr.T,tr.P,geoPrm.DS3+14.8,reacPos);
		if(YYHit) ETmp = yd.ELoss(ncl,ETmp,tr.T);
		if(CsIHit) ETmp = csi.ELoss(ncl,ETmp,tr.T);
		if(Sd1Hit) ETmp = sd1.ELoss(ncl,ETmp,tr.T);
		if(Sd2Hit) ETmp = sd2.ELoss(ncl,ETmp,tr.T);
	}	
	
	return (mask && shield && YYHit && CsIHit);
}

void sortEnergies(){
	if(yd.mul>1) yd.SortByEnergy();
	if(csi.mul>1) csi.SortByEnergy();
	if(sd1.mul>1) sd1.SortByEnergy();
	if(sd2.mul>1) sd2.SortByEnergy();
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

	det.SSB=SSBdE;
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
}
