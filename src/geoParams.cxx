#include "geoParams.h"
#include <stdio.h>
#include <stdlib.h>

geoParams::geoParams(){
  //geoParams::Class()->IgnoreTObjectStreamer;
  geoParams::Clear();
}

void geoParams::ReadParams(char* line)
{
	bool expect_val=true;
	char *from=line;
	char *to=line;
	while (*from) {
		if (*from>32) {*to=*from;to++;}
		from++;
	}
	*to=0;
	if (*line==0) return; // line is empty
	
	char* val=strchr(line,'=');
	if (!val){ 
		val=strchr(line, '!');
		expect_val=false;
	}
	if (!val) printf("Missing = or ! in input file, line: '%s'\n",line);
	*val=0;
	
	// trim param name
	char* trm=val-1;
	while (*trm<=32) *(trm--)=0;
	
	val++;
	if (*val==0 && expect_val) printf("Value missing for parameter %s",line);

	char cval[256];	
	double fval;	
	std::string strval;
	sscanf(val,"%lf",&fval);
	sscanf(val,"%s",cval);
	strval=cval;
	
	//	parameter of type string:

	if (strcmp(line,"Bs")==0){
	   Bs=fval;
	}
	if (strcmp(line,"ICPressure")==0){
	   ICPressure=fval;
	}
	if (strcmp(line,"TFoil")==0){
	   TFoil=fval;
	}
	if (strcmp(line,"AoZFoil")==0){
	   AoZFoil=fval;
	}
	if (strcmp(line,"TTgt")==0){
	   TTgt=fval;
	}
	if (strcmp(line,"DYY")==0){
	   DYY=fval;
	}
	if (strcmp(line,"DYU")==0){
	   DYYU=fval;
	}
	if (strcmp(line,"TYY1")==0){
	   TYY[0]=fval;
	}
	if (strcmp(line,"TYY2")==0){
	   TYY[1]=fval;
	}
	if (strcmp(line,"TYY3")==0){
	   TYY[2]=fval;
	}
	if (strcmp(line,"TYY4")==0){
	   TYY[3]=fval;
	}
	if (strcmp(line,"TYY5")==0){
	   TYY[4]=fval;
	}
	if (strcmp(line,"TYY6")==0){
	   TYY[5]=fval;
	}
	if (strcmp(line,"TYY7")==0){
	   TYY[6]=fval;
	}
	if (strcmp(line,"TYY8")==0){
	   TYY[7]=fval;
	}
	if (strcmp(line,"TYU1")==0){
	   TYYU[0]=fval;
	}
	if (strcmp(line,"TYU2")==0){
	   TYYU[1]=fval;
	}
	if (strcmp(line,"TYU3")==0){
	   TYYU[2]=fval;
	}
	if (strcmp(line,"TYU4")==0){
	   TYYU[3]=fval;
	}
	if (strcmp(line,"TYU5")==0){
	   TYYU[4]=fval;
	}
	if (strcmp(line,"TYU6")==0){
	   TYYU[5]=fval;
	}
	if (strcmp(line,"TYU7")==0){
	   TYYU[6]=fval;
	}
	if (strcmp(line,"TYU8")==0){
	   TYYU[7]=fval;
	}
	if (strcmp(line,"DSD")==0){
	   DS3=fval;
	}
	if (strcmp(line,"DSU")==0){
	   DS3U=fval;
	}
	if (strcmp(line,"TSD1")==0){
	   TS3[0]=fval;
	}
	if (strcmp(line,"TSD2")==0){
	   TS3[1]=fval;
	}
	if (strcmp(line,"TSU")==0){
	   TS3U=fval;
	}
	if (strcmp(line,"MASK")==0){
	   Mask=(fval==1.);
	}
	if (strcmp(line,"SHIELD")==0){
	   Shield=(fval==1.);
	}
	if (strcmp(line,"ORIENTATION")==0){
	   Orientation=(fval==1.);
	}
	if (strcmp(line,"FOIL")==0){
	   MFoil=strval;	
	}
	if (strcmp(line,"TARGET")==0){
	   MTgt=strval;	
	}
}

void geoParams::Load(std::string filename){	

	char line[256];
	FILE* file=fopen(filename.data(),"rb");
	if (!file)
	{
		printf("Cannot open geoParams file '%s' for reading. Stop.\n",filename.data());
		exit(0);
	}
	
	printf("Reading geoParams file '%s'\n",filename.data());
	
	while (!feof(file))
	{
		if (!fgets(line,256,file)) break;
		printf("%s",line);
		// skip leading white spaces
		char* ptr=line;
		while ((*ptr>0) && (*ptr<32)) ptr++;
		//printf("%s\n",ptr[0]);
		switch (ptr[0])
		{
			case 0   :
			case '#' :
			case '/' :  continue;
			default  :  ReadParams(ptr);
		}
	}
	fclose(file);
	file=NULL;
}

void geoParams::Print(){

	printf("********************************\n\n");
	printf("IC:\t Pressure: %.2lf Torr\n",ICPressure);
	printf("Foil:\t Material: %s\tthickness: %.2lf mg/cm2\tA/Z: %.3f\n",MFoil.data(),TFoil,AoZFoil);
	printf("Target:\t Material: %s\tthickness: %.1lf mg/cm2\tA/Z: %.3f\n",MTgt.data(),TTgt,AoZTgt);
	printf("Size of beamspot: %.1lf mm\n",Bs);
	printf("Distance target-YY: %.1lf mm\n",DYY);
	printf("Distance target-upstream YY: %.1lf mm\n",DYYU);
	printf("Thickness YY1: %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\n",TYY[0],TYY[1],TYY[2],TYY[3],TYY[4],TYY[5],TYY[6],TYY[7]);
	printf("Thickness upstream YY1: %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\n",TYYU[0],TYYU[1],TYYU[2],TYYU[3],TYYU[4],TYYU[5],TYYU[6],TYYU[7]);
	printf("Distance target-S3: %.1lf mm\n",DS3);
	printf("Distance target-upstream S3: %.1lf mm\n",DS3U);
	printf("Thickness 1st S3: %.1lf um\n",TS3[0]);
	printf("Thickness 2nd S3: %.1lf um\n",TS3[1]);
	printf("Thickness upstream S3: %.1lf um\n",TS3U);
	printf("Using heat shield: %s\n",Shield ? "Yes" : "No");
	printf("Using target mask: %s\n",Mask ? "Yes" : "No");
	printf("Target %s from foil.\n",Orientation ? "upstream" : "downstream");
	printf("********************************\n\n");
}

void geoParams::Clear(){

	DYY=0.; 
	DS3=0.; 
	DYYU=0.; 
	DS3U=0.; 
	TS3U=0.; 
	Bs=0.; 
	TFoil=0.; 
	AoZFoil=0.; 
	TTgt=0.; 
	AoZTgt=0.; 
	ICPressure=0.; 
	Mask=kTRUE; 
	Shield=kTRUE; 
	for(Int_t i=0; i<8; i++) { TYY[i] = 0.; }
	for(Int_t i=0; i<8; i++) { TYYU[i] = 0.; }
	for(Int_t i=0; i<2; i++) { TS3[i] = 0.; }
}
