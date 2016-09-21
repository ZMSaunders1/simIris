#ifndef __YYHIT_H
#define __YYHIT_H

#include "TObject.h"
#include "TClass.h"
#include "TVector3.h"
 /****YY1 hit***/

class YYHit : public TObject{
 	public:
  		Int_t mul;
  		Double_t fX[4];
  		Double_t fY[4];
  		Double_t fZ[4];//should be equal to distance to YY1
  		Double_t fPhiCalc[4];//Calculated phi (using Seg)
  		Double_t fThetaCalc[4];//Calculated theta (using YdRing)
  		Double_t dE[4];//Energy loss
  		Bool_t hit[4];//hits YY1
  		Int_t Seg[4];
  		Int_t Ring[4];
 	public:
  		YYHit();//! Create
  		virtual ~YYHit() {} //!

  		//YYHit(const YYHit &);                          // The copy constructor.
  		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Calculate(Double_t, Double_t, Double_t, TVector3);  //!
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(YYHit,1);
};

#endif
