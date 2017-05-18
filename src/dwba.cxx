#include "dwba.h"

Double_t eval_theta(Double_t in, Double_t x[181], Double_t y[181])
{
	Double_t dxin=0., dx=0., dy=0., e=0.;
	if(in<=0.){
		e = 0.;
	}
	else if(in<x[0]){
		e = y[0]*in/x[0];
	}
	else if(in>x[180]){
		dxin = in-x[180];
		dx = x[180]-x[179];
		dy = y[180]-y[179];
		e = y[180]+dy*dxin/dx;
	}
	else{
		for(int i=1; i<181;i++){   
			if(in>x[i-1]&&in<x[i]){
				dxin = in-x[i-1];
				dx = x[i]-x[i-1];
				dy = y[i]-y[i-1];
				e = y[i-1]+dy*dxin/dx;
				break;
			}
		}
	}
	return e;
}

Double_t load_dwba(std::string filename, Double_t x[181], Double_t y[181])
{
	Double_t ymax=0.;
	FILE *dFile;
	char line[32];
	std::string cTheta, cXsec;
	dFile = fopen(filename.data(),"r");
	if(dFile == NULL) printf("Could not open DWBA file!\n");
	else{
		for(Int_t i=0; i<181; i++){
			if (!fgets(line,32,dFile)) break;
			//printf("%s",line);
			cTheta.assign(line,line+7);
			cXsec.assign(line+8,line+31);
			sscanf(cTheta.data(),"%lf",&x[i]);
			sscanf(cXsec.data(),"%lf",&y[181-i]);
			printf("%lf\t%lf\n",x[i],y[i]);
			if(y[i]>ymax) ymax = y[i];
		}
	}
	return ymax;
}
