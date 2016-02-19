#include <stdio.h>
#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include "TGraph.h"
#include "TVector3.h"
#include "TLorentzVector.h"

double dalitzCalc(double *x, double* par) //calculated 2 dim dalitz function, arguments are m13 and m23 invariant masses
{
	double m13 = x[0];//12 invariant mass squared
	double m23 = x[1];//23 invariant mass squared
	double E = par[0];// total energy in CM                  
	double m1 = par[1];
	double m2 = par[2];
	double m3 = par[3];
	double m12 = E*E+m1*m1+m2*m2+m3*m3-m13-m23;  
	double E1 = (m13+m12-m2*m2-m3*m3)/(2*E);
	double E2 = (m12+m23-m1*m1-m3*m3)/(2*E);
	double E3 = (m13+m23-m1*m1-m2*m2)/(2*E);
 	if ((E1<m1) || (E2<m2) || (E3<m3)) return 0;
   	double P1 = sqrt(E1*E1-m1*m1);
   	double P2 = sqrt(E2*E2-m2*m2);
 	double P3 = sqrt(E3*E3-m3*m3);
 	if ((P1+P2<P3) || (P1+P3<P2) || (P2+P3<P1)) return 0; //momentum conservation 

 	return 1;
}

double P3P1(double *x, double* par) //calculated P3 vs P1
{
  	double P3 = x[0];//12 invariant mass squared
  	double P1 = x[1];//23 invariant mass squared
  	double E = par[0];// total energy in CM                  
  	double m1 = par[1];
  	double m2 = par[2];
  	double m3 = par[3];

  	double E1 = sqrt(P1*P1+m1*m1);
  	double E3 = sqrt(P3*P3+m3*m3);
  	double E2 = E- E1 -E3;

 	if ((E1<m1) || (E2<m2) || (E3<m3)) return 0;

 	double P2 = sqrt(E2*E2-m2*m2);
 	if ((P1+P2<P3) || (P1+P3<P2) || (P2+P3<P1)) return 0; //momentum conservation 
 	double  eps = sqrt((E-E3)*(E-E3)-P3*P3);
 	double  p1 = sqrt((eps*eps-(m2-m1)*(m2-m1))*(eps*eps-(m2+m1)*(m2+m1)))/(2*eps);//momentum of particle 1 in the 1+2 CM  

 	//deltas
 	double dP1 = 0.01;
 	double dP3 = 0.01;
 	double dE1 = sqrt((P1+dP1)*(P1+dP1)+m1*m1)-sqrt(P1*P1+m1*m1);
 	double dE3 = sqrt((P3+dP3)*(P3+dP3)+m3*m3)-sqrt(P3*P3+m3*m3);
 	double dE2 = 0.-dE1-dE3;// conseravtion of energy
 	double dP2 =  sqrt((E2+dE2)*(E2+dE2)-m2*m2)-sqrt(E2*E2-m2*m2);
 	double theta = acos((P1*P1+P3*P3-P2*P2)/(2*P1*P3)); 
 	double theta2 = acos(((P1+dP1)*(P1+dP1)+(P3+dP3)*(P3+dP3)-(P2+dP2)*(P2+dP2))/(2*(P1+dP1)*(P3+dP3)));
 	double dtheta = theta2-theta; 

 	if (!(dP2*dP2>0)) return 0;
 	//  printf("parameters %lf %lf %lf %lf\n",dP2,P2,theta,theta2);
 	return P3*P3*p1*dP3/E3/E1/E2;
}

double phaseSpace(double *x, double* par)//phase space function
{
	//x is the B kinetic energy
	//par[0] is the total energy
	double t3 = x[0];//kinetic energy of particle 3, or B
	double E = par[0];// total energy in CM
	double m1 = par[1]; // b
	double m2 = par[2]; // c 
	double m3 = par[3]; // B
	double p3 = sqrt(t3*t3+2*t3*m3);//momentum
	double E3 = sqrt(p3*p3+m3*m3);
	double  eps = sqrt((E-E3)*(E-E3)-p3*p3);
	double norm = 1.;
	//double norm = 6E6;
	if (E-m1-m2-m3<0) return 0;
	
	if (eps-m1-m2<0) return 0;
	 double  p1 = sqrt((eps*eps-(m2-m1)*(m2-m1))*(eps*eps-(m2+m1)*(m2+m1)))/(2*eps);//momentum of particle 1 in the 1+2 CM
	 //  printf("parameters %lf %lf %lf %lf\n",par[0],par[1],par[2],par[3]);  
	 //  return p3*p3*p1*p1/norm;
	 return p3*p1/eps/norm;
}//end phase space function


double phaseSpace2b(double *x, double* par)//2-body phase space function       
{
	double eps = x[0];
	double m1 = par[0];
	double m2 = par[1];
	double  p1 = sqrt((eps*eps-(m2-m1)*(m2-m1))*(eps*eps-(m2+m1)*(m2+m1)))/(2*eps);//momentum of particle 1 in the 1+2 CM          
	if (p1 >0) return p1*sqrt(p1*p1+m1*m1);
	else return 0;
}//end phase space test function 

TGraph* calcKin(double ma, double mA, double mb, double mB, double TA)//calculate elastic kinematics
{
	TLorentzVector LVa = TLorentzVector(0,0,0,ma);
	TLorentzVector LVA = TLorentzVector(0,0,sqrt(TA*TA+2*mA*TA),mA+TA);
	TLorentzVector LVb = TLorentzVector(0,0,0,mb);
	
	const int NPoints =89;
	const Double_t thstep = 2.0;
	Double_t betacm, ETotCM, bfE,bfP;
	Double_t bfT[NPoints];
	Double_t thcm;
	Double_t theta[NPoints];//lab angle
	betacm = LVA.Z()/(ma+TA+mA);
	LVA.Boost(0,0,-1.0*betacm);
	LVa.Boost(0,0,-1.0*betacm);
	ETotCM = LVA.T()+LVa.T();
	bfE = (ETotCM*ETotCM+mb*mb-mB*mB)/(2.*ETotCM);
	bfP = sqrt(bfE*bfE-mb*mb);
	for (int thn = 0; thn <NPoints; thn++ ){//thetacm number
		thcm = thstep*(thn+1);
		LVb.SetXYZT(0,0,bfP,bfE);
		
		LVb.SetTheta(thcm*TMath::DegToRad());
		LVb.Boost(0,0,betacm);
		bfT[thn] = LVb.E()-mb;
		theta[thn]=LVb.Theta()*TMath::RadToDeg();
	}//for thn
	TGraph *gr = new TGraph(NPoints,theta,bfT);
	gr->SetMarkerStyle(21);
	gr->SetTitle("kinematics plot");
	gr->SetName("grKin");
	return gr;
}

Double_t calcQv(double ma, double mA, double mb, double Eb, double EBeam, double bTheta)//calculate Q-value , Eb is the kinetic energy of b
{
	Double_t thetaR = bTheta*TMath::DegToRad();
	Double_t PA = sqrt(EBeam*EBeam+2*EBeam*mA);//momentum of A
	Double_t Pb = sqrt(Eb*Eb+2*Eb*mb);//momentum of b
	Double_t  Q = mA+ma-mb- sqrt(mA*mA+mb*mb-ma*ma-2.*(mA+EBeam)*(mb+Eb)+2.*PA*Pb*cos(thetaR)+2.*(EBeam+mA+ma-Eb-mb)*ma);
	return Q;
}

Int_t Calc3body(Double_t ETotCM, TLorentzVector *part1,TLorentzVector *part2,TLorentzVector *part3)
// A function for calculating energies and momenta of 3 bodies given total energy in CM and energy of one particle (part1) (the directions of the other two will be chosen isotropically in their own CM)
{
	TRandom3 fRndm(0);
	// Double_t m1 = part1->M();
	Double_t m2 = part2->M();
	Double_t m3 = part3->M();
	if (ETotCM-part1->E()-m2-m3<0){
		std::cout << "energy not conserved! "<< std::endl;  
		return 0;
	}
	
	Double_t x1,y1,z1;
	fRndm.Sphere(x1,y1,z1,part1->Rho());//change P1 in an isotropically random direction

	part1->SetX(x1);
	part1->SetY(y1);
	part1->SetZ(z1);
	Double_t beta23x = -1*part1->Px()/(ETotCM-part1->E()); //CM betax of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E
	Double_t beta23y = -1*part1->Py()/(ETotCM-part1->E()); //CM betay of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E
	Double_t beta23z = -1*part1->Pz()/(ETotCM-part1->E()); //CM betaz of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E	
	
	TVector3 abeta = TVector3(beta23x,beta23y,beta23z);
	Double_t  E23 = sqrt((ETotCM-part1->E())*(ETotCM-part1->E())-part1->Rho()*part1->Rho()); //CM energy of d+n in their own CM  

	if (E23-m2-m3<-1E-6){
	  	printf("Error: E23 is smaller than m2+m3 (%lf)!\n",E23-m2-m3);
	}
	Double_t E2_23 = (E23*E23+m2*m2-m3*m3)/(2*E23); //Deuteron energy in the d+n CM	
	Double_t E3_23 = E23 - E2_23;
	Double_t P2_23 = sqrt(E2_23*E2_23-m2*m2);
	/**********2+3 CM isotropic simulation**********/
	//fRndm.SetSeed(11);
	fRndm.Sphere(x1,y1,z1,P2_23);//change it in an isotropically random direction
	
	part3->SetPxPyPzE(-1*x1,-1*y1,-1*z1,E3_23);
	
	part2->SetPxPyPzE(x1,y1,z1,E2_23);
	//printf("Rho = %lf %lf %lf\n",x1,y1,E23-m2-m3);                                                                                                                                             
	//going to general CM frame                                                       
	part2->Boost(abeta);
	part3->Boost(abeta);
	return 1;
}

Int_t Calc3bodyPS(Double_t ETotCM, TLorentzVector *part1,TLorentzVector *part2,TLorentzVector *part3)
// A function for calculating energies and momenta of 3 bodies given total energy in CM and energy of one particle (part1) (the directions of the other two will be chosen isotropically in their own CM)
{
	TRandom3 fRndm(0);
	Double_t m1 = part1->M();
	Double_t m2 = part2->M();
	Double_t m3 = part3->M();
	if (ETotCM-m1-m2-m3<0){
	  	std::cout << "Calc3bodyPS energy not conserved! "<< ETotCM << " " << m1+m2+m3<<std::endl;
		return 0;
	}
  	Double_t s12, s23, s31, E1, E2, E3, P1, P2, P3;
	Bool_t valid = 0;

	while (!valid){
   		s12 = fRndm.Uniform((m1+m2)*(m1+m2),(ETotCM-m3)*(ETotCM-m3));
   		s23 = fRndm.Uniform((m3+m2)*(m3+m2),(ETotCM-m1)*(ETotCM-m1));
   		s31 = ETotCM*ETotCM+m1*m1+m2*m2+m3*m3-s12-s23;

   		E1 = (ETotCM*ETotCM+m1*m1-s23)/(2*ETotCM);
   		E2 = (ETotCM*ETotCM+m2*m2-s31)/(2*ETotCM);
		//  E1 = fRndm.Uniform(m1,ETotCM-m2-m3);
		// E2 = fRndm.Uniform(m2,ETotCM-m1-m3);
  		E3 = ETotCM - E1 - E2;
  		P1 = sqrt(E1*E1-m1*m1);
  		P2 = sqrt(E2*E2-m2*m2);
  		if (E3>m3)	P3 = sqrt(E3*E3-m3*m3);
  		else {E3 =0;P3 = -10E6;}
  		if ((P1+P2<P3) || (P1+P3<P2) || (P2+P3<P1) || (E3 < m3)) valid = 0;
  		else valid = 1;
 	}

	part1->SetX(0);
	part1->SetY(0);
	part1->SetZ(P1);
	part1->SetE(E1);
	
	part2->SetX(0);
	part2->SetY(0);
	part2->SetZ(P2); 
	part2->SetE(E2);
	
	part3->SetX(0);
	part3->SetY(0);
	part3->SetZ(P3);
	part3->SetE(E3);
	// TVector3 abeta = TVector3(0,0,1);

 	Double_t theta12 = acos((P1*P1+P2*P2-P3*P3)/(-2*P1*P2));//angle between P1 and P2
 	part2->SetTheta(theta12);
 	Double_t phi12 = fRndm.Uniform(-1*TMath::Pi(),TMath::Pi());
 	part2->SetPhi(phi12);

 	Double_t theta13 = acos((P1*P1+P3*P3-P2*P2)/(-2*P1*P3));//angle between P1 and P2                                                 
 	part3->SetTheta(theta13);
 	part3->SetPhi(phi12+TMath::Pi());//opposite of part2 phi
 	// printf("Rhox =  %lf %lf %lf %lf\n",part1->X(),part2->X(),part3->X(), part1->X()+part2->X()+part3->X());         

	Double_t x1,y1,z1;
	fRndm.Sphere(x1,y1,z1,P1);//change P1 in an isotropically random direction                         
	part1->SetXYZT(x1,y1,z1,E1);                                             
	
	//rotate part2 and part3 with the same theta and phi as part1
 	part2->RotateY(part1->Theta());
	part2->RotateZ(part1->Phi());
 	part3->RotateY(part1->Theta());
 	part3->RotateZ(part1->Phi());

	return 1;
}

Int_t Calc2body(TLorentzVector *part1,TLorentzVector *part2,TLorentzVector *part3)
// A function for calculating energies and momenta of 2-body breakup of part1 into part2 and part3 (the directions of the other two will be chosen isotropically in their own CM) 
{
	TRandom3 fRndm(0);
	//Double_t m1 = part1->M();
	Double_t m2 = part2->M();
	Double_t m3 = part3->M();
	if (part1->E()-m2-m3<0){
		std::cout << "Calc2body energy not conserved! "<< std::endl;
		return 0;
	}
	
	Double_t beta23x = part1->Px()/(part1->E()); //CM betax of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E
	Double_t beta23y = part1->Py()/(part1->E()); //CM betay of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E  
	Double_t beta23z = part1->Pz()/(part1->E()); //CM betaz of d+n (or 9Li +n) in the general 9Li+d+n CM system //beta = P/E                                                                             
	
	TVector3 abeta = TVector3(beta23x,beta23y,beta23z);
	Double_t  E23 = sqrt(part1->E()*part1->E()-part1->Rho()*part1->Rho()); //CM energy of d+n in their own CM                                                                                     
	if (E23-m2-m3<0){
	  	printf("Error: E23 is smaller than m2+m3!\n");
	}
	Double_t E2_23 = (E23*E23+m2*m2-m3*m3)/(2*E23); //Deuteron energy in the d+n CM 
	Double_t E3_23 = E23 - E2_23;
	Double_t P2_23 = sqrt(E2_23*E2_23-m2*m2);
	/**********2+3 CM isotropic simulation**********/
	Double_t x1,y1,z1;
	fRndm.Sphere(x1,y1,z1,P2_23);//change it in an isotropically random direction	
	part3->SetPxPyPzE(-1*x1,-1*y1,-1*z1,E3_23);
	part2->SetPxPyPzE(x1,y1,z1,E2_23);
	 // printf("Rho = %lf %lf %lf\n",x1,y1,E23-m2-m3);
	//going to general CM frame	
	part2->Boost(abeta);
	part3->Boost(abeta);
	return 1;
}

Int_t Calc4bodyPS(Double_t ETotCM, TLorentzVector *part1,TLorentzVector *part2,TLorentzVector *part3,TLorentzVector *part4)//this function doesn't work right now
// A function for calculating energies and momenta of 4 bodies given total energy in CM and energy of one particle (part1) (the directions of the other two will be chosen isotropically in their own CM)
{
	TRandom3 fRndm(0);
	Double_t m1 = part1->M();
	Double_t m2 = part2->M();
	Double_t m3 = part3->M();
	Double_t m4 = part4->M();
	if (ETotCM-m1-m2-m3-m4<0){
	  	std::cout << "Calc3bodyPS energy not conserved! "<< ETotCM << " " << m1+m2+m3+m4<<std::endl;
		return 0;
	}
	Double_t  s12, s23, s34, s234,s341;
	Double_t  m34=0.;
	
	TLorentzVector part34(0,0,0,m34);
	// Calc3bodyPS(ETotCM,part1,part2,&part34);
	Double_t E1, E2, E34, P1, P2, P34;
	Bool_t valid = 0;
	while (!valid){
		s34 = fRndm.Uniform((m3+m4)*(m3+m4),(ETotCM-m1-m2)*(ETotCM-m1-m2));
		m34 = sqrt(s34);
		s12 = fRndm.Uniform((m1+m2)*(m1+m2),(ETotCM-m34)*(ETotCM-m34));
		s234 = fRndm.Uniform((m34+m2)*(m34+m2),(ETotCM-m1)*(ETotCM-m1));
		s341 = ETotCM*ETotCM+m1*m1+m2*m2+s34-s12-s234;
		
		E1 = (ETotCM*ETotCM+m1*m1-s234)/(2*ETotCM);
		E2 = (ETotCM*ETotCM+m2*m2-s341)/(2*ETotCM);
		
		E34 = ETotCM - E1 - E2;
		P1 = sqrt(E1*E1-m1*m1);
		P2 = sqrt(E2*E2-m2*m2);
		if (E34>m34) P34 = sqrt(E34*E34-m34*m34);
		else {E34 =0;P34 = -10E6;}
		if ((P1+P2<P34) || (P1+P34<P2) || (P2+P34<P1) || (E34 < m34)) valid = 0;
		else valid = 1;
	}
	
	part1->SetX(0);
	part1->SetY(0);
	part1->SetZ(P1);
	part1->SetE(E1);
	
	part2->SetX(0);
	part2->SetY(0);
	part2->SetZ(P2);
	part2->SetE(E2);
	
	part34.SetX(0);
	part34.SetY(0);
	part34.SetZ(P34);
	part34.SetE(E34);
	// TVector3 abeta = TVector3(0,0,1);	
	Double_t theta12 = acos((P1*P1+P2*P2-P34*P34)/(-2*P1*P2));//angle between P1 and P2
	part2->SetTheta(theta12);
	Double_t phi12 = fRndm.Uniform(-1*TMath::Pi(),TMath::Pi());
	part2->SetPhi(phi12);
	
	Double_t theta13 = acos((P1*P1+P34*P34-P2*P2)/(-2*P1*P34));//angle between P1 and P2
	part34.SetTheta(theta13);
	part34.SetPhi(phi12+TMath::Pi());//opposite of part2 phi
	// printf("Rhox =  %lf %lf %lf %lf\n",part1->X(),part2->X(),part3->X(), part1->X()+part2->X()+part3->X());
	
	Double_t x1,y1,z1;
	fRndm.Sphere(x1,y1,z1,P1);//change P1 in an isotropically random direction
	part1->SetXYZT(x1,y1,z1,E1);
	
	//rotate part2 and part3 with the same theta and phi as part1
	part2->RotateY(part1->Theta());
	part2->RotateZ(part1->Phi());
	part34.RotateY(part1->Theta());
	part34.RotateZ(part1->Phi());
	
	Calc2body(&part34,part3,part4);
	
	return 1;
}
                              
