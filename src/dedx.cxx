#include "dedx.h"

dedx::dedx()
{
	Clear();	
}

void dedx::Clear()
{
	for(int i=0; i<100; i++)
	{
		eFoil[i] = 0.;
		eAl[i] = 0.;	
		eB[i] = 0.;	
		eC4H10[i] = 0.;	
		eCsI[i] = 0.;
		eTgt[i] = 0.;	
		eMy[i] = 0.;	
		eP[i] = 0.;	
		eSi[i] = 0.;	
		eSi3N4[i] = 0.;	
		eSiO2[i] = 0.;	

		dedxFoil[i] = 0.;	
		dedxAl[i] = 0.;	
		dedxB[i] = 0.;	
		dedxC4H10[i] = 0.;	
		dedxCsI[i] = 0.;
		dedxTgt[i] = 0.;	
		dedxMy[i] = 0.;	
		dedxP[i] = 0.;	
		dedxSi[i] = 0.;	
		dedxSi3N4[i] = 0.;	
		dedxSiO2[i] = 0.;	
	}
}

void dedx::loadIncomingELoss(std::string path, std::string name, std::string foil, std::string tgt, Double_t mass)
{
	printf("Loading %s/lise_%s_in_%s.txt\n",path.data(),name.data(),foil.data());	
	loadELoss(Form("%s/lise_%s_in_%s.txt",path.data(),name.data(),foil.data()),eFoil,dedxFoil,mass/1000.);	
	printf("Loading %s/lise_%s_in_C4H10.txt\n",path.data(),name.data());	
	loadELoss(Form("%s/lise_%s_in_C4H10.txt",path.data(),name.data()),eC4H10,dedxC4H10,mass/1000.);	
	printf("Loading %s/lise_%s_in_%s.txt\n",path.data(),name.data(),tgt.data());	
	loadELoss(Form("%s/lise_%s_in_%s.txt",path.data(),name.data(),tgt.data()),eTgt,dedxTgt,mass/1000.);	
	printf("Loading %s/lise_%s_in_Si3N4.txt\n\n",path.data(),name.data());	
	loadELoss(Form("%s/lise_%s_in_Si3N4.txt",path.data(),name.data()),eSi3N4,dedxSi3N4,mass/1000.);	
}

void dedx::loadOutgoingELoss(std::string path, std::string name, std::string foil, std::string tgt, Double_t mass)
{
	printf("Loading %s/lise_%s_in_%s.txt\n",path.data(),name.data(),foil.data());	
	loadELoss(Form("%s/lise_%s_in_%s.txt",path.data(),name.data(),foil.data()),eFoil,dedxFoil,mass/1000.);	
	printf("Loading %s/lise_%s_in_Al.txt\n",path.data(),name.data());
	loadELoss(Form("%s/lise_%s_in_Al.txt",path.data(),name.data()),eAl,dedxAl,mass/1000.);
	printf("Loading %s/lise_%s_in_B.txt\n",path.data(),name.data());
	loadELoss(Form("%s/lise_%s_in_B.txt",path.data(),name.data()),eB,dedxB,mass/1000.);
	printf("Loading %s/lise_%s_in_CsI.txt\n",path.data(),name.data());	
	loadELoss(Form("%s/lise_%s_in_CsI.txt",path.data(),name.data()),eCsI,dedxCsI,mass/1000.);	
	printf("Loading %s/lise_%s_in_%s.txt\n",path.data(),name.data(),tgt.data());
	loadELoss(Form("%s/lise_%s_in_%s.txt",path.data(),name.data(),tgt.data()),eTgt,dedxTgt,mass/1000.);
	printf("Loading %s/lise_%s_in_My.txt\n",path.data(),name.data());	
	loadELoss(Form("%s/lise_%s_in_My.txt",path.data(),name.data()),eMy,dedxMy,mass/1000.);	
	printf("Loading %s/lise_%s_in_P.txt\n",path.data(),name.data());
	loadELoss(Form("%s/lise_%s_in_P.txt",path.data(),name.data()),eP,dedxP,mass/1000.);
	printf("Loading %s/lise_%s_in_Si.txt\n",path.data(),name.data());	
	loadELoss(Form("%s/lise_%s_in_Si.txt",path.data(),name.data()),eSi,dedxSi,mass/1000.);	
	printf("Loading %s/lise_%s_in_SiO2.txt\n\n",path.data(),name.data());
	loadELoss(Form("%s/lise_%s_in_SiO2.txt",path.data(),name.data()),eSiO2,dedxSiO2,mass/1000.);
}

void dedx::loadELoss(std::string filename, Double_t e[100], Double_t dedx[100], Double_t m)
{
	std::ifstream infile;
	Char_t line[2000];
			
	infile.open(filename.data()); 
	if(!infile.is_open()){
		printf("Cannot open the file %s!!\n",filename.data());
		exit(0);
	}

  	infile.getline(line,2000);
	for(UInt_t i=0; i<100; i++){
    	infile.getline(line,2000);
    	sscanf(line,"%*f\t%*f\t%lf\t%lf\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\t\n",&e[i],&dedx[i]);
  		e[i] *= m/0.931494061; //stopping power for this mass
   	}    
	infile.close();
}
	
