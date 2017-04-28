#ifndef __YYHIT_H
#define __YYHIT_H

#include "TObject.h"
#include "TClass.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "nucleus.h"
#include "eloss.h"

class YYHit : public TObject{
 	public:
		Double_t Thickness[8]; // in um
		Double_t Avg_Thickness; // in um
  		Int_t mul;
  		std::vector<Double_t> fX;
  		std::vector<Double_t> fY;
  		std::vector<Double_t> fZ;//should be equal to distance to YY1
  		std::vector<Double_t> fPhiCalc;//Hitd phi (using Seg)
  		std::vector<Double_t> fThetaCalc;//Hitd theta (using YdRing)
  		std::vector<Double_t> fThetaRand;//Hitd theta (using YdRing)
  		std::vector<Double_t> dE;//Energy loss
  		std::vector<Double_t> dE_ideal;//Energy loss
  		//std::vector<Bool_t> hit;//hits YY1
  		std::vector<Int_t> Seg;
  		std::vector<Int_t> Ring;
 	public:
  		YYHit();//! Create
  		virtual ~YYHit() {} //!

  		//YYHit(const YYHit &);                          // The copy constructor.
  		void Init(Double_t[8]);//! Create
  		Double_t ThetaMin(Double_t);  //!
  		Double_t ThetaMax(Double_t);  //!
  		Bool_t Hit(Double_t, Double_t, Double_t, TVector3);  //!
		Double_t ELoss(nucleus, Double_t, Double_t);
		void SortByEnergy();
  		void Clear();  //!
	protected:

 	private:
  
	ClassDef(YYHit,1);
};

#endif
