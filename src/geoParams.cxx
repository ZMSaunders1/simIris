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
	if (strcmp(line,"TAg")==0){
	   TAg=fval;
	}
	if (strcmp(line,"TTgt")==0){
	   TTgt=fval;
	}
	if (strcmp(line,"DYY")==0){
	   DYY=fval;
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
	if (strcmp(line,"DS3")==0){
	   DS3=fval;
	}
	if (strcmp(line,"TS31")==0){
	   TS3[0]=fval;
	}
	if (strcmp(line,"TS32")==0){
	   TS3[1]=fval;
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
	printf("Ag-Foil thickness: %.2lf um\n",TAg);
	printf("Target thickness: %.1lf um\n",TTgt);
	printf("Size of beamspot: %.1lf mm\n",Bs);
	printf("Distance target-YY: %.1lf mm\n",DYY);
	printf("Thickness YY1: %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\t %.1lf um\n",TYY[0],TYY[1],TYY[2],TYY[3],TYY[4],TYY[5],TYY[6],TYY[7]);
	printf("Distance target-S3: %.1lf mm\n",DS3);
	printf("Thickness 1st S3: %.1lf um\n",TS3[0]);
	printf("Thickness 2nd S3: %.1lf um\n",TS3[1]);
	printf("********************************\n\n");
}

void geoParams::Clear(){

	DYY=0.; 
	DS3=0.; 

}
