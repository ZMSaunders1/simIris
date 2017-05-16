#ifndef __DEDX_H
#define __DEDX_H

#include "TObject.h"
#include "TClass.h"
#include <fstream>
#include <string>

class dedx : public TObject{
 	public: 	
		
		Double_t eFoil[100];	
		Double_t eAl[100];	
		Double_t eB[100];	
		Double_t eC4H10[100];	
		Double_t eCsI[100];
		Double_t eTgt[100];	
		Double_t eMy[100];	
		Double_t eP[100];	
		Double_t eSi[100];	
		Double_t eSi3N4[100];	
		Double_t eSiO2[100];	

		Double_t dedxFoil[100];	
		Double_t dedxAl[100];	
		Double_t dedxB[100];	
		Double_t dedxC4H10[100];	
		Double_t dedxCsI[100];
		Double_t dedxTgt[100];	
		Double_t dedxMy[100];	
		Double_t dedxP[100];	
		Double_t dedxSi[100];	
		Double_t dedxSi3N4[100];	
		Double_t dedxSiO2[100];	

 	public:
  		dedx();//! Create
  		virtual ~dedx() {} //!
		
		void loadIncomingELoss(std::string, std::string, std::string, std::string, Double_t);
		void loadOutgoingELoss(std::string, std::string, std::string, std::string, Double_t);
		void loadELoss(std::string, Double_t[100], Double_t[100], Double_t);
  		void Clear();  //!
	protected:

 	private:
  
	//ClassDef(dedx,1);
};

#endif
