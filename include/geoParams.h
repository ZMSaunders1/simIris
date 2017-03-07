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
		
		Double_t Tt; 
		Double_t Bs; 
		Double_t DYY; 
		Double_t DS3; 

		//virtual void ReadCalibPar(char* line);
		virtual void ReadParams(char* line);
		virtual void Load(std::string filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(geoParams,1)
};

#endif
// end
