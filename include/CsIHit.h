#ifndef __CsI1PARTICLE_H
#define __CsI1PARTICLE_H

#include "TObject.h"
#include "TClass.h"
#include "TVector3.h"
 /****CsI hit***/

class CsIHit : public TObject{
 	public:
  		Int_t mul;
  		Double_t fX[4];
  		Double_t fY[4];
  		Double_t fZ[4];//should be equal to distance to CsI
  		Double_t fPhiCalc[4];//Calculated phi (using Seg)
  		Double_t dE[4];//Energy loss
  		Bool_t hit[4];//hits CsI
  		Int_t Seg[4];
 	public:
  		CsIHit();//! Create
  		virtual ~CsIHit() {} //!

  		//CsIHit(const CsIHit &);                          // The copy constructor.
   		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Calculate(Double_t, Double_t, Double_t, TVector3);  //!
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(CsIHit,1);
};

#endif
