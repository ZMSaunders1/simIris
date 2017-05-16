// geoParams.h

#ifndef geoParams_H
#define geoParams_H
#include <TObject.h>
#include <TClass.h>
#include <string>

//Extern
//extern int gMesytecnitems;
class geoParams : public TObject {
	public:
		geoParams(); 
		virtual ~geoParams() {} //! 
		
		Double_t Bs; 
		Double_t TFoil; 
		Double_t TTgt; 
		Double_t DYY; 
		Double_t TYY[8]; 
		Double_t DS3; 
		Double_t TS3[2];
		std::string MFoil;
		std::string MTgt;	
		Bool_t Mask;
		Bool_t Shield;

		//virtual void ReadCalibPar(char* line);
		virtual void ReadParams(char* line);
		virtual void Load(std::string filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(geoParams,1)
};

#endif
// end
