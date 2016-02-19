// params.h

#ifndef params_H
#define params_H
#include <TObject.h>
#include <TClass.h>
#include <string>

//Extern
//extern int gMesytecnitems;
class params : public TObject {
	public:
		params(); 
		virtual ~params() {} //! 

		std::string A; 
		std::string a;
		std::string B;
		std::string b;
		std::string c;
		std::string d;

		Double_t E;
		Double_t R;
		Double_t W;
		Double_t Tt; 
		Double_t Bs; 
		
		Double_t DYY; 
		Double_t DS3; 

		Int_t N;

		//virtual void ReadCalibPar(char* line);
		virtual void ReadParams(char* line);
		virtual void Load(std::string filename);
		virtual void Print();
		virtual void Clear();
//		ClassDef(params,1)
};

#endif
// end
