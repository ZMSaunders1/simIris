#include "params.h"
#include <stdio.h>
#include <stdlib.h>

params::params(){
  //params::Class()->IgnoreTObjectStreamer;
  params::Clear();
}

void params::ReadParams(char* line)
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
	if (strcmp(line,"N")==0){
	   N=int(fval);
	}
	if (strcmp(line,"A")==0){
	   A=strval;	
	}
	if (strcmp(line,"a")==0){
	   a=strval;
	}
	if (strcmp(line,"B")==0){
	   B=strval;
	}
	if (strcmp(line,"b")==0){
	   b=strval;
	}
	if (strcmp(line,"c")==0){
	   c=strval;
	}
	if (strcmp(line,"d")==0){
	   d=strval;
	}
	if (strcmp(line,"R")==0){
	   R=fval;
	}
	if (strcmp(line,"W")==0){
	   W=fval;
	}
	if (strcmp(line,"E")==0){
	   E=fval;
	}
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

void params::Load(std::string filename){	

	char line[256];
	FILE* file=fopen(filename.data(),"rb");
	if (!file)
	{
		printf("Cannot open params file '%s' for reading. Stop.\n",filename.data());
		exit(0);
	}
	
	printf("Reading params file '%s'\n",filename.data());
	
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

void params::Print(){
	printf("\n********************************\n");
	switch(N){
		case 2:
			printf("%s(%s,%s)%s(%.1lf) @ %.1lf MeV, %.1lf mu, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R,E,Tt,W);
			break;
		case 3:
			printf("%s(%s,%s)%s(%.1lf)+%s @ %.1lf MeV, %.1lf mu, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R,c.data(),E,Tt,W);
			break;
		case 4:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s @ %.1lf MeV, %.1lf mu, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R,c.data(),d.data(),E,Tt,W);
			break;
		default:
			printf("%s(%s,%s)%s(%.1lf)+%s+%s @ %.1lf MeV, %.1lf mu, Width=%.2lf\n",A.data(),a.data(),b.data(),B.data(),R,c.data(),d.data(),E,Tt,W);
			break;
	}
	printf("********************************\n\n");
	printf("Size of beamspot: %.1lf mm\n",Bs);
	printf("Distance target-YY: %.1lf mm\n",DYY);
	printf("Distance target-S3: %.1lf mm\n",DS3);
	printf("********************************\n\n");
}

void params::Clear(){

	A.clear(); 
	a.clear();
	B.clear();
	b.clear();
	c.clear();
	d.clear();

	E=0.;
	R=0.;
	W=0.;
	Tt= 0.; 
	Bs=0.; 
	
	DYY=0.; 
	DS3=0.; 

	N=0;

}
