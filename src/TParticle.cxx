//#include <TLorentzRotation.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TParticle.h"
//#include "PhysicalConstants.h"

ClassImp(TParticle);


//______________________________________________________________________________          
TParticle::TParticle(const TParticle &sp) : TObject(sp)
{
  // -- Copy constructor.                                                                 
  ((TParticle&)sp).Copy(*this);
}


//______________________________________________________________________________          
void TParticle::Copy(TObject &sp) const
{
  // -- Copy this method.                                                                 

  TObject::Copy((TObject&)sp);

  ((TParticle&)sp).fEulerPhi   = fEulerPhi;
  ((TParticle&)sp).fEulerTheta = fEulerTheta;
  ((TParticle&)sp).fEulerPsi   = fEulerPsi;

  //  ((TParticle&)sp).Clear(); //Why do we need this?                                
}



//______________________________________________________________________________          
void TParticle::Clear(Option_t *option)
{
  // --                                                                                   
  //                                                                                      

  SetA(0);
  SetZ(0);

 

     fT  = 0;                                                                             
     fE  = 0;                                                                             
     fPx = 0;                                                                             
     fPy = 0;                                                                             
     fPz = 0;                                                                             
 
}


// //______________________________________________________________________________       
// Double_t TParticle::GetBeta()                                                      
// {                                                                                      
//   TLorentzVector vec(fE,fPx,fPy,fPz);                                                  
//   return vec.Beta();                                                                   
// }            


//______________________________________________________________________________          
TVector3 TParticle::GetBeta()
{
  // -- Get beta for this particle.                                                       

  Get4Vec();
  fBeta.SetXYZ(fPx/fE,fPy/fE,fPz/fE);
  return fBeta;
}


//______________________________________________________________________________          
TLorentzVector TParticle::Get4Vec()
{
  // -- Get the 4-vector for this particle.                                               

  fPVec.SetPxPyPzE(fPx,fPy,fPz,fE);
  return fPVec;
}


// //______________________________________________________________________________       
// Double_t TParticle::GetRelVel()                                                    
// {                                                                                      
//   //                                                                                   

//   //  T05133Event *event = (T05133Event*)fEvent.GetObject();                           
//   TVector3 cmVec = fEvent->GetVelCM();                                                 

//   Double_t sMass = fEvent->fSp.fA*kAmu;                                                
//   TVector3 sVec(fEvent->fSp.fPx/sMass,fEvent->fSp.fPy/sMass,fEvent->fSp.fPz/sMass);    

//   Double_t mass = fA*kAmu;                                                             
//   TVector3 hVec(fPx/mass,fPy/mass,fPz/mass);                                           
//   TVector3 vec = hVec-cmVec;                                                           

//   return vec.Mag();                                                                    
// }                                                                                      


//______________________________________________________________________________          
Double_t TParticle::GetTCM(TVector3 boost)
{
  // -- Get relativistic kinetic energy for this particle.                                

  Double_t tCM=0;

  // Initialize the Lorentz Transformation with the given boost.                          
  // TLorentzRotation rotCM(boost);
  // Invert the transformation as we are going to the rod frame.                          
  // rotCM.Invert();
  // Get the particle's 4-vector.                                                         
  TLorentzVector vec = Get4Vec();
  // Transform 4-vector into COM system.                                                  
  // vec = rotCM * vec;
  // In the COM system the spatial momentum is 0.                                         
  // The relativistic K.E. is then Etot - Rest Energy.                                    
  tCM = vec.E() - GetA()*kAmu;

  return tCM;
}

//______________________________________________________________________________          
TVector3 TParticle::GetVelVec()
{
  //                                                                                      

  Double_t mass = fA*kAmu;
  fVelVec.SetXYZ(fPx/mass,fPy/mass,fPz/mass);

  return fVelVec;
}

//______________________________________________________________________________          
Double_t TParticle::GetVelMag()
{
  //                                                                                      

  TVector3 vec = GetVelVec();
  return vec.Mag();
}


//______________________________________________________________________________          
void TParticle::SetName(const Char_t *name)
{
  // -- Change (i.e. set) the name of the TParticle.                                  
  //                                                                                      

  fName = name;
}

//______________________________________________________________________________          
void TParticle::SetNameTitle(const Char_t *name, const Char_t *title)
{
  // -- Change (i.e. set) all the TParticle parameters (name and title).              
  //                                                                                      

  fName  = name;
  fTitle = title;
}

//______________________________________________________________________________          
void TParticle::SetTitle(const Char_t *title)
{
  // -- Change (i.e. set) the title of the TParticle.                                 
  //                                                                                      

  fTitle = title;
}


//Make it a method for a particle class.
//Double_t TParticle::eloss(Double_t ein, Double_t th , TF1 *func)//initial energy and thickness are given as arguments 
// {

//  // for (Int_t )

//  //Energy loss calculation
//  Double_t dx =0.1; //in microns
//  Double_t de = 0; //energy loss
//  Double_t en= ein; //the energy variable
//  Double_t pos = 0.;// the position variable


 
//  while (pos<= th){
//    de =  (dx*func->Eval(en))/2.; //
// en = en - de;
//    de =(dx*func->Eval(en))/2.;//energy loss in dx
//    en = en -de; // energy remaining after dx

//    pos = pos +dx;

//  }
//  return ein - en;
// }


//Make it a method for a particle class.
//Double_t TParticle::thickness(Double_t ein, Double_t efi , TF1 *func)//initial energy and final energy are given as arguments, calculates target thickness 
// {

//  // for (Int_t )

//  //Energy lossx calculation
//   Double_t dx = 0; 
//   Double_t de = 0.01; //energy loss step in MeV
//  Double_t en= ein; //the energy variable
//  Double_t th = 0.;// the thickness variable

//  //Integrate numerically 
 
//  while (en >= efi){
//    dx =  de/func->Eval(en)/2.; //
//    th = th +dx;
//    en = en-de;
//   dx =  de/func->Eval(en)/2.; //
//   th = th +dx;
//  }
//  return th;
// }


