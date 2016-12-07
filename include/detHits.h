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

	//if(mask && shield){
	if(shield){
		if(YYHit) ETmp = yd.ELoss(ncl,ETmp,tr.T);
		if(CsIHit) ETmp = csi.ELoss(ncl,ETmp,tr.T);
		if(Sd1Hit) ETmp = sd1.ELoss(ncl,ETmp,tr.T);
		if(Sd2Hit) ETmp = sd2.ELoss(ncl,ETmp,tr.T);
	}	
	
	return (mask && shield && YYHit && CsIHit);
}

