#include "header.h"

Bool_t detHits(PTrack tr, nucleus ncl, TVector3 reacPos)
{
	Bool_t mask = maskClear(tr.T,tr.P);
	Bool_t shield = shieldClear(tr.T,tr.P);
	Bool_t YYHit = yd.Hit(tr.T,tr.P,prm.DYY,reacPos) ;
	Bool_t CsIHit = csi.Hit(tr.T,tr.P,prm.DYY+11.6,reacPos) ;
	Bool_t Sd1Hit = sd1.Hit(tr.T,tr.P,prm.DS3,reacPos);
	Bool_t Sd2Hit = sd2.Hit(tr.T,tr.P,prm.DS3+14.8,reacPos);

	Double_t ETmp = tr.Ebt;

	if(mask && shield){
	//if(shield){
		if(YYHit) ETmp = yd.ELoss(ncl,ETmp,tr.T);
		if(CsIHit) ETmp = csi.ELoss(ncl,ETmp,tr.T);
		if(Sd1Hit) ETmp = sd1.ELoss(ncl,ETmp,tr.T);
		if(Sd2Hit) ETmp = sd2.ELoss(ncl,ETmp,tr.T);
	}	
	
	return (mask && shield && YYHit && CsIHit);
}
void setIDet(Double_t ICdE, Double_t SSBdE)
{
	if(yd.mul>0)
	{
  		det.TYdMul=1;  
		det.TYdEnergy.push_back(yd.dE[0]);
  		det.TYdTheta.push_back(yd.fThetaRand[0]);// Yd theta angle                                                                       
		det.TYdChannel.push_back(yd.Seg[0]*16+yd.Ring[0]);
		det.TYdNo.push_back(yd.Seg[0]);
		det.TYdRing.push_back(yd.Ring[0]);
	}

	if(csi.mul>0)
	{
		det.TCsI1Mul=1;
  		det.TCsI1Energy.push_back(csi.dE[0]);
		det.TCsI1Channel.push_back(csi.Seg[0]);
		det.TCsI2Mul=1;
  		det.TCsI2Energy.push_back(csi.dE[0]);
		det.TCsI2Channel.push_back(csi.Seg[0]);
	}

	det.SSB=SSBdE;
   	det.TICEnergy.push_back(ICdE);
	det.TICChannel.push_back(15);

	if(sd1.mul>0)
	{
		det.TSd1rMul=1;
  		det.TSd1rEnergy.push_back(sd1.dE[0]);
		det.TSd1rChannel.push_back(sd1.Ring[0]);
  		det.TSdTheta.push_back(sd1.fThetaCalc[0]);
		det.TSd1sMul=1;
  		det.TSd1sEnergy.push_back(sd1.dE[0]);
		det.TSd1sChannel.push_back(sd1.Seg[0]);
  		det.TSdPhi.push_back(sd1.fPhiCalc[0]);
	}
	if(sd2.mul>0)
	{
		det.TSd2rMul=1;
  		det.TSd2rEnergy.push_back(sd2.dE[0]);
		det.TSd2rChannel.push_back(sd2.Ring[0]);
		det.TSd2sMul=1;
  		det.TSd2sEnergy.push_back(sd2.dE[0]);
		det.TSd2sChannel.push_back(sd2.Seg[0]);
	}
}
