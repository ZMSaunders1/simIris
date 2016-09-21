// Eloss.h

#ifndef EneLoss_H
#define EneLoss_H

#include "TGraph.h"
#include "nucleus.h"
//Extern

Double_t eval(Double_t, Double_t[100], Double_t[100]);
Double_t eloss(nucleus, Double_t, Double_t, Double_t, Double_t[100], Double_t[100]);
Double_t straggling(nucleus, Double_t, Double_t); //Energy straggling    
#endif
// end
