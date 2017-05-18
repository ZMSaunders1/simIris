#include "eloss.h"

Double_t eval(Double_t in, Double_t x[100], Double_t y[100])
{
	Double_t dxin=0., dx=0., dy=0., e=0.;
	if(in<=0.){
		e = 0.;
	}
	else if(in<x[0]){
		e = y[0]*in/x[0];
	}
	else if(in>x[99]){
		dxin = in-x[99];
		dx = x[99]-x[98];
		dy = y[99]-y[98];
		e = y[99]+dy*dxin/dx;
	}
	else{
		for(int i=1; i<100;i++){   
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

//Make it a method for a particle class.
Double_t eloss(nucleus P, Double_t TZoverA, Double_t ein, Double_t th , Double_t x[100], Double_t y[100])//initial energy and thickness are given as arguments 
{

	Double_t k;
	Double_t Bohr;
	Double_t sgm;
	Double_t strg;

	if(ein==0.) return 0;
	if(th==0.) return 0.;
	
	TRandom3 *rndm = new TRandom3(0);
	//Energy loss calculation including energy Straggling
	Double_t dx =th/100.; //in mg/cm2
	Bohr = TMath::Sqrt(0.157 * dx * P.Z*P.Z *TZoverA)/1000.;
	Double_t de = 0; //energy loss
	Double_t en= ein; //the energy variable
	for (int i=0; i<100; i++){
	  	de = (dx * eval(en,x,y));//energy loss in dx
		if(de>en){
		   	en=0.;	
			break;
		}
		k= 1.1+0.47*TMath::Log10(en/Double_t(P.A));
		sgm = k*Bohr;
		strg = rndm->Gaus(de,sgm);
		de = (strg>0.) ? strg : 0.;
	  	if(de>en){
		   	en=0.;	
			break;
		}
		en = en - de; // energy remaining after dx
	}
	rndm->Delete();
	return ein-en;
}

Double_t elossFi(Double_t efi, Double_t th, Double_t x[100], Double_t y[100])//final energy and thickness are given as arguments 
{
	if (th==0) return efi;
	//Energy loss calculation
	Double_t dx =th/100.; //in 
	Double_t de = 0; //energy loss
	Double_t en= efi; //the energy variable
	Double_t pos = 0.;// the position variable
	
	while (pos<= th){
		de =  (dx* eval(en,x,y))/2.; //
		en = en + de;
		de =(dx* eval(en,x,y))/2.;//energy loss in dx
		en = en +de; // energy remaining after dx
		pos = pos +dx;
	}
 	return en-efi;
}
