#ifndef __CsI1PARTICLE_H
#define __CsI1PARTICLE_H

#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "nucleus.h"
#include "eloss.h"
 /****CsI hit***/

class CsIHit : public TObject{
	public:
		Double_t Thickness;
  		Int_t mul;
  		Double_t fX[4];
  		Double_t fY[4];
  		Double_t fZ[4];//should be equal to distance to CsI
  		Double_t fPhiCalc[4];//Hitd phi (using Seg)
  		Double_t dE[4];//Energy loss
  		Bool_t hit[4];//hits CsI
  		Int_t Seg[4];
 	public:
  		CsIHit();//! Create
  		CsIHit(Double_t);//! Create
  		virtual ~CsIHit() {} //!

  		//CsIHit(const CsIHit &);                          // The copy constructor.
   		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Hit(Double_t, Double_t, Double_t, TVector3);  //!
		Double_t ELoss(nucleus, Double_t, Double_t);
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(CsIHit,1);
};

#endif
