//#include <stdio.h>
//#include <stdint.h>
#include <TMath.h>

// ****Heat shield clearance***
// The heat shield is a cylinder with Y-symmetry axis.
// The cylinder's outer diameter is 129 mm.
// The cylinder has a round corner rectangular window. 
// The dimension of the window cut is x= 2*67.4 mm and y = 2*50 mm, with corner radius of 20 mm. 
// Edit: A different heat shield was installed, where the vertical dimension of
// the cutout is not 2*50mm but instead 50mm + 111.7mm (above and below the y=0,
// respectively). Also, I don't think the 'clearSquare' is necessary any longer.
// (This size is before bending the window material to the cylindrical shape.)
//variables for the elliptic mask hole
Bool_t shieldClear(Double_t theta, Double_t phi)
{
  const Double_t horSide = 2*67.4;	//horizontal side
  const Double_t verSide = 2*50.0;	//vertical side
  const double verSideUp = 50.;
  const double verSideDown = 111.7;
  const Double_t cylRad = 129.0/2;	//cylinder radius
  //Double_t alphaR = 0;

  Bool_t clear = 0;	// does the particle clear the shield?
  Bool_t clearRect = 0;	// does the particle clear the rectangular shield (without the curves)?
  Bool_t clearCirc = 0;	// does the particle clear curve circle of the shield?
  Bool_t clearSquare = 0;	// does the particle clear square at the edge of the shield, with the opposite edge at the centre of the curve circle?
  Double_t cyCurve = 20;	//mm  

  Double_t polRad = cylRad/(sqrt(sin(theta)*cos(phi)*sin(theta)*cos(phi)+cos(theta)*cos(theta))); //polar radius, i.e. distance from the reaction point
// 	if(cos(phi)!=0) alphaR = acos(polRad*cos(theta)/cylRad)*fabs(cos(phi))/cos(phi); //alpha in cylindrical coordinates.
//
//	alphaR = asin(polRad*sin(theta)*cos(phi)/cylRad);	//alpha in cylindrical coordinates.
//	Double_t xOnShield = cylRad*alphaR;	//distance on the cylinder from 0.
  Double_t xOnShield = polRad*sin(theta)*cos(phi);	//distance on the cylinder from 0.
  Double_t yOnShield = polRad*sin(theta)*sin(phi);	//distance on the cylinder from 0.
	//clearRect = ((fabs(xOnShield)<horSide/2) && (fabs(yOnShield)<verSide/2));
  //clearCirc = pow((fabs(xOnShield)-horSide/2+cyCurve),2)+ pow((fabs(yOnShield)-verSide/2+cyCurve),2) < cyCurve*cyCurve;
  //clearSquare = fabs(fabs(xOnShield)-horSide/2+cyCurve/2)<cyCurve/2 &&  fabs(fabs(yOnShield)-verSide/2+cyCurve/2)<cyCurve/2;
  if(yOnShield < 0){
    clearRect = ((fabs(xOnShield)<horSide/2) && (fabs(yOnShield)<verSideDown));
    clearCirc = yOnShield > -(verSideDown-cyCurve)
      || fabs(xOnShield) < (horSide/2-cyCurve) ? true
      : pow((fabs(xOnShield)-horSide/2+cyCurve),2) 
      + pow((fabs(yOnShield)-verSideDown+cyCurve),2) < cyCurve*cyCurve;
  }
  else{
    clearRect = ((fabs(xOnShield)<horSide/2) && (fabs(yOnShield)<verSideUp));
    clearCirc = yOnShield < verSideUp - cyCurve
      || fabs(xOnShield) < (horSide/2-cyCurve) ? true
      : pow((fabs(xOnShield)-horSide/2+cyCurve),2) 
      + pow((fabs(yOnShield)-verSideUp+cyCurve),2) < cyCurve*cyCurve;
  }
	//  printf("xOnshield = %lf, yOnShield =%lf , xs= %lf \n",xOnShield, yOnShield,  cylRad*asin(polRad*sin(theta)*cos(phi)/cylRad));
	//clear = clearRect && (clearCirc || !clearSquare);
  clear = clearRect && clearCirc;
  return clear;
}
//****end Heat shield clearance***

//*** Mask clear******************
Bool_t maskClear(Double_t theta, Double_t phi)
{
	const Double_t radX = 54.0;
	const Double_t radY = 25.5;
	const Double_t maskDistance =129./4;

	Bool_t clear = 0;// does the particle clear the mask?
	Double_t xx = maskDistance*tan(theta)*cos(phi);  //x on the mask plane
	Double_t yy = maskDistance*tan(theta)*sin(phi);  //y on the mask plane
	clear = ((xx*xx/(radX*radX)+yy*yy/(radY*radY)) < 1.);
	
	return clear;
}
//********End Mask clear********

