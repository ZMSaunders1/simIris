#ifndef __S3HIT_H
#define __S3HIT_H

#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "nucleus.h"
#include "eloss.h"

class S3Hit : public TObject{
 	public:
		Bool_t Orientation; // 0 = rings first, 1 = sectors first
		Double_t Thickness; // in um
  		Int_t mul;
  		Double_t fX[4];
  		Double_t fY[4];
  		Double_t fZ[4];//should be equal to distance to S3
  		Double_t fPhiCalc[4];//Hitd phi (using Seg)
  		Double_t fThetaCalc[4];//Hitd theta (using YdRing)
  		Double_t dE[4];//Energy loss
  		Double_t dE_ideal[4];//Energy loss
  		Bool_t hit[4];//hits S3
  		Int_t Seg[4];
  		Int_t Ring[4];
 	public:
  		S3Hit();//! Create
  		S3Hit(Bool_t, Double_t);//! Create
  		virtual ~S3Hit() {} //!

  		//S3Hit(const S3Hit &);                          // The copy constructor.
   		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Hit(Double_t, Double_t, Double_t, TVector3);  //!
		Double_t ELoss(nucleus, Double_t, Double_t);
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(S3Hit,1);
};

#endif
