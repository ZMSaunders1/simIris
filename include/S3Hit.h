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
  		std::vector<Double_t> fX;
  		std::vector<Double_t> fY;
  		std::vector<Double_t> fZ;//should be equal to distance to S3
  		std::vector<Double_t> fPhiCalc;//Hitd phi (using Seg)
  		std::vector<Double_t> fPhiRand;//Hitd phi (using Seg)
  		std::vector<Double_t> fThetaCalc;//Hitd theta (using YdRing)
  		std::vector<Double_t> fThetaRand;//Hitd theta (using YdRing)
  		std::vector<Double_t> dE;//Energy loss
  		std::vector<Double_t> dE_ideal;//Energy loss
  		//std::vector<Bool_t> hit;//hits S3
  		std::vector<Int_t> Seg;
  		std::vector<Int_t> Ring;
 	public:
  		S3Hit();//! Create
  		virtual ~S3Hit() {} //!

  		//S3Hit(const S3Hit &);                          // The copy constructor.
  		void Init(Bool_t, Double_t);//! Create
   		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Hit(Double_t, Double_t, Double_t, TVector3);  //!
		Double_t ELoss(nucleus, Double_t, Double_t);
		void SortByEnergy();
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(S3Hit,1);
};

#endif
