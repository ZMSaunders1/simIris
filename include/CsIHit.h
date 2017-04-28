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
  		std::vector<Double_t> fX;
  		std::vector<Double_t> fY;
  		std::vector<Double_t> fZ;//should be equal to distance to CsI
  		std::vector<Double_t> fPhiCalc;//Hitd phi (using Seg)
  		std::vector<Double_t> fPhiRand;//Hitd phi (using Seg)
  		std::vector<Double_t> dE;//Energy loss
  		std::vector<Double_t> dE_ideal;//Energy loss
  		//std::vector<Bool_t> hit;//hits CsI
  		std::vector<Int_t> Seg;
 	public:
  		CsIHit();//! Create
  		CsIHit(Double_t);//! Create
  		virtual ~CsIHit() {} //!

  		//CsIHit(const CsIHit &);                          // The copy constructor.
   		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Hit(Double_t, Double_t, Double_t, TVector3);  //!
		Double_t ELoss(nucleus, Double_t, Double_t);
		void SortByEnergy();
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(CsIHit,1);
};

#endif
