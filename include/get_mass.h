#include <stdio.h>
#include <stdlib.h>
#include <string>

Double_t get_mass(std::string input)
{
	const Double_t Amu=931.494061;
	FILE *massfile;
	massfile=fopen("/home/iris/mholl/mass.mas12","r");
	if(massfile==NULL){
		printf("Cannot open mass table !!\n");
		exit(0);
	}

	char buffer[124];

	char inputEl[2];
	Int_t inputA;

	std::string cN,cZ,cA,cEl,cME,cAM;
	Int_t N,Z,A;
	Double_t ME,AM;
	Double_t mass1=0.;
	Double_t mass2=0.;
	//Double_t factor=1.;

	if(input=="p") input="1H";
	else if(input=="d") input="2H";
	else if(input=="t") input="3H";
	else if(input=="n") input="1n";

	sscanf(input.data(),"%d%s",&inputA,inputEl);
	
	std::string inpEl(inputEl);
	inpEl.erase( std::remove_if( inpEl.begin(), inpEl.end(), ::isspace ), inpEl.end() );

	printf("%d%s\t===================================================================\n",inputA,inpEl.data());

	while(!feof(massfile)){
		if(!fgets(buffer,124,massfile)) break;

		cN.assign(buffer+6,buffer+9);
		cN.erase( std::remove_if( cN.begin(), cN.end(), ::isspace ), cN.end() );
		cZ.assign(buffer+11,buffer+14);
		cZ.erase( std::remove_if( cZ.begin(), cZ.end(), ::isspace ), cZ.end() );
		cA.assign(buffer+16,buffer+19);
		cA.erase( std::remove_if( cA.begin(), cA.end(), ::isspace ), cA.end() );
		cEl.assign(buffer+20,buffer+22);
		cEl.erase( std::remove_if( cEl.begin(), cEl.end(), ::isspace ), cEl.end() );
		cME.assign(buffer+30,buffer+41);
		cME.erase( std::remove_if( cME.begin(), cME.end(), ::isspace ), cME.end() );
		cAM.assign(buffer+96,buffer+112);
		cAM.erase( std::remove_if( cAM.begin(), cAM.end(), ::isspace ), cAM.end() );
		sscanf(cN.data(),"%d",&N);
		sscanf(cZ.data(),"%d",&Z);
		sscanf(cA.data(),"%d",&A);
		sscanf(cME.data(),"%lf",&ME);
		sscanf(cAM.data(),"%lf",&AM);
		mass1=Amu*A+ME/1000.;
		mass2=Amu*AM/1000000.;
		if(A==inputA&&cEl==inpEl) {
			printf("\tN=%d\t Z=%d\t A=%d\t%s\t%lf\t%lf\t \n",N,Z,A,cEl.data(),mass1,mass2);
			break;
		}
	}
	fclose(massfile);
	return mass1;
}
