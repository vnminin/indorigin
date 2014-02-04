#ifndef _cltctmc_PCLT_H
#define _cltctmc_PCLT_H

#include <RcppArmadillo.h>

double logspaceAdd(const double loga, const double logb);
double logspaceSubtract(const double loga, const double logb);
void hypgeoF11_k_2k(const double x, arma::subview_col<double> result, double * work_bi);
void hypgeoF11_k_2kp1(const double x, arma::subview_col<double> result, const double * work_bi);
void hypgeoF11_kp1_2kp1(const double x, arma::subview_col<double> result, const double * work_bi);
void pclt(const double t, const double lambda_1, const double lambda_2,
        arma::subview_col<double> result);
#endif
