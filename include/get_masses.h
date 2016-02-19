#include <stdio.h>
#include <string.h>

void get_masses(string input)
{
	const double Amu=931.494061;
	FILE *file;
	file=fopen("/home/iris/mholl/mass_amdc2013.txt","r");
	char buffer[124];

	//string input="6He";
	char inputEl[3];
	int inputA;

	sscanf(input.data(),"%d%s",&inputA,inputEl);
	//inputEl[2]='\0';
	printf("%d\t%s\n",inputA,inputEl);

	int N,Z,A;
	double ME,AM;
	double mass1, mass2;

	while(!feof(file)){
		if(!fgets(buffer,124,file)) break;

		//for(int i=0; i<100; i++){
		char cN[4],cZ[4],cA[4],cEl[3],cME[12],cAM[16];

		strncpy(cN,buffer+6,sizeof(cN));
		cN[3] ='\0';
		strncpy(cZ,buffer+11,sizeof(cZ));
		cZ[3] ='\0';
		strncpy(cA,buffer+16,sizeof(cA));
		cA[3] ='\0';
		strncpy(cEl,buffer+20,sizeof(cEl));
		cEl[2] ='\0';
		strncpy(cME,buffer+30,sizeof(cME));
		cME[11] ='\0';
		strncpy(cAM,buffer+96,sizeof(cAM));
		strncpy(cAM+3,buffer+100,sizeof(cAM)-3);
		cAM[15] ='\0';

		sscanf(cN,"%d",&N);
		sscanf(cZ,"%d",&Z);
		sscanf(cA,"%d",&A);
		sscanf(cME,"%lf",&ME);
		sscanf(cAM,"%lf",&AM);
		mass1=Amu*A+ME/1000.;
		mass2=Amu*AM/1000000.;
		if(A==inputA&&strcmp(cEl,inputEl)==0) printf("N=%d\t Z=%d\t A=%d\t%s\t%lf\t%lf\t \n",N,Z,A,cEl,mass1,mass2);
		//printf("N=%d\t Z=%d\t A=%d\t%s\tmass excess=%lf\t atomic mass=%lf\t \n",N,Z,A,cEl,ME,AM);
	}
	fclose(file);
}
