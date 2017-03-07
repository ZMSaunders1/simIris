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
	if (strcmp(line,"Tt")==0){
	   Tt=fval;
	}
	if (strcmp(line,"DYY")==0){
	   DYY=fval;
	}
	if (strcmp(line,"DS3")==0){
	   DS3=fval;
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
	printf("Target thickness: %.1lf mm\n",Tt);
	printf("Size of beamspot: %.1lf mm\n",Bs);
	printf("Distance target-YY: %.1lf mm\n",DYY);
	printf("Distance target-S3: %.1lf mm\n",DS3);
	printf("********************************\n\n");
}

void geoParams::Clear(){

	DYY=0.; 
	DS3=0.; 

}
