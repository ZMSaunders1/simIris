// Eloss.h
#ifndef EneLoss_H
#define EneLoss_H

#include <TMath.h>
#include <TRandom3.h>
#include "nucleus.h"

Double_t eval(Double_t, Double_t[100], Double_t[100]);
Double_t eloss(nucleus, Double_t, Double_t, Double_t, Double_t[100], Double_t[100]);
Double_t elossFi(Double_t, Double_t, Double_t[100], Double_t[100]);
#endif
// end
