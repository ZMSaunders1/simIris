//header for simulation functions
#ifndef __SIMFUNCS_H
#define __SIMFUNCS_H

double dalitzCalc(double *,double*);//calculated 2 dim dalitz plot

double P3P1(double *, double*); //calculated P3 vs P1

double phaseSpace(double *, double*);//phase space function

double phaseSpace2b(double *, double*);//2-body phase space function

TGraph* calcKin(double, double, double, double, double);//calculate elastic kinematics  

Double_t calcQv(double ma, double mA, double mb, double bfT, double EBeam, double bTheta);//calculate Q-value       
Int_t Calc3body(Double_t, TLorentzVector *,TLorentzVector *,TLorentzVector *);
Int_t Calc2body(TLorentzVector *,TLorentzVector *,TLorentzVector *);
Int_t Calc3bodyPS(Double_t, TLorentzVector *,TLorentzVector *,TLorentzVector *);
Int_t Calc4bodyPS(Double_t ETotCM, TLorentzVector *part1,TLorentzVector *part2,TLorentzVector *part3,TLorentzVector *part4);
#endif
