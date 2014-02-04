#include "pclt.h"
using namespace Rcpp;

inline static double debyeU1(const double * tpow) {
    return (3.0*tpow[1] - 5.0*tpow[3])/24.0;
}

inline static double debyeU2(const double * tpow) {
    return (81.0*tpow[2] - 462.0*tpow[4] + 385.0*tpow[6])/1152.0;
}

inline static double debyeU3(const double * tpow) {
    return (30375.0*tpow[3] - 369603.0*tpow[5] + 765765.0*tpow[7] - 425425.0*tpow[9])/414720.0;
}

inline static double debyeU4(const double * tpow) {
    return (4465125.0*tpow[4] - 94121676.0*tpow[6] + 349922430.0*tpow[8] -
            446185740.0*tpow[10] + 185910725.0*tpow[12])/39813120.0;
}

inline static int getMaxK(const int n) {
    if (n % 2 == 0)
        return n / 2;
    return (n + 1) / 2;
}

double logspaceAdd(const double loga, const double logb) {
    if (!R_FINITE(loga))
        return logb;
    if (loga > logb)
        return logspaceAdd(logb, loga);
    return logb + log1p(exp(loga - logb));
}

double logspaceSubtract(const double loga, const double logb) {
    try {
        if (logb >= loga) {
            throw std::invalid_argument("Must have b < a for logspace subtraction.");
        }
        if(!R_FINITE(logb) && logb < 0.0)
            return loga;
        return loga + log1p(-exp(logb - loga));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    // not reached
    return NA_REAL;
}

static double scaledLogBesselLargeNu(const double nu, const double x) {
    //Use A&S 9.7.7 to compute log(exp(-x) * I_nu) for large values of nu.
    const double z = x / nu;
    const double root_term = sqrt(1.0 + z * z);
    const double t = 1.0 / root_term;
    const double eta = root_term + log(z / (1.0 + root_term));
    double sum = 0.0;
    double tpow[13];
    tpow[0] = 1.0;
    for(int i = 1; i < 13; i++) tpow[i] = tpow[i-1] * t;
    sum = debyeU1(tpow) / nu + debyeU2(tpow) / (nu * nu) +
            debyeU3(tpow) / (nu * nu * nu) + debyeU4(tpow) / (nu * nu * nu * nu);
    return (-.5 * log(2.0 * M_PI * nu * root_term) + nu * eta + log1p(sum) - x);
}

void hypgeoF11_k_2k(const double x, arma::subview_col<double> result, double *work_bi) {
// Compute 1F1(k, 2k)
    const int len = result.n_elem;
    if ((len == 1) || (x <= 0.0))
        return;
    const int nmax = len - 1;
    const int kmax = getMaxK(nmax);
    const double y = x / 2.0;
    const double lx = log(.25 * x);

    // compute log bessel function values.
    R::bessel_i_ex(y, (kmax + 1) - .5, 2.0, work_bi);
    int k = kmax;
    while((work_bi[k] == 0.0 || !R_FINITE(work_bi[k])) && k >= 0) {
        work_bi[k] = scaledLogBesselLargeNu((k + 1) - .5, y);
        k--;
    }
    while(k >= 0 && work_bi[k] > 0.0) {
        work_bi[k] = log(work_bi[k]);
        k--;
    }
    double lgamma_kp1 = R::lgammafn(kmax + 1 + .5);
    double lgamma_k = lgamma_kp1 - log(kmax + .5);
    double lgamma_km1 = 0.0;

    // initialize F11, then move down to the
    result[2 * kmax - 1] = work_bi[kmax - 1] + lgamma_k + (.5 - kmax) * lx + x;
    for(k = kmax; k > 1; k--) {
        lgamma_km1 = lgamma_k - log((k - 1) + .5);
        result[2 * (k - 1) - 1] = work_bi[k-2] + lgamma_km1 + (.5 - (k - 1)) * lx + x;
        lgamma_k = lgamma_km1;
    }
    return;
}

void hypgeoF11_kp1_2kp1(const double x, arma::subview_col<double> result, const double *work_bi) {
    const int len = result.n_elem;
    if ((len == 1) || (x <= 0.0))
        return;
    const int nmax = len - 1;
    const int kmax = getMaxK(nmax);

    if (nmax % 2 == 0) {
        const double lgkp1 = R::lgammafn(kmax + 1 + .5);
        result[2 * kmax] = logspaceAdd(result[2 * kmax - 1],
                    work_bi[kmax] + lgkp1 + (.5 - (kmax + 1))*log(x) +
                    x + log(x / (4.0 * kmax + 2.0)));
    }
    for (int k = kmax; k > 1; k--) {
        result[2 * (k - 1)] = logspaceAdd(result[2 * (k - 1) - 1],
                log(x / (4.0 * (k - 1) + 2.0)) + result[2 * k - 1]);
    }
    return;
}

void hypgeoF11_k_2kp1(const double x, arma::subview_col<double> result, const double *work_bi) {
    const int len = result.n_elem;
    if ((len == 1) || (x <= 0.0))
        return;
    const int nmax = len - 1;
    const int kmax = getMaxK(nmax);

    if (nmax % 2 == 0) {
        const double lgkp1 = R::lgammafn(kmax + 1 + .5);
        result[2 * kmax] = logspaceSubtract(result[2 * kmax - 1],
                    work_bi[kmax] + lgkp1 + (.5 - (kmax + 1)) * log(x) +
                    x + log(x / (4.0 * kmax + 2.0)));
    }
    for (int k = kmax; k > 1; k--) {
        result[2 * (k - 1)] = logspaceSubtract(result[2 * (k - 1) - 1],
                log(x / (4.0 * (k - 1) + 2.0)) + result[2 * k - 1]);
    }
    return;
}

void pclt(const double t, const double lambda_1, const double lambda_2,
        arma::subview_col<double> result){
    if (lambda_1 == lambda_2) {
        for(int k = 0; k < result.n_elem; k++)
            result[k] = R::dpois(k, t * lambda_1, 1);
        return;
    }
    const int len = result.n_elem;
    const int nmax = len - 1;
    const int kmax = getMaxK(nmax);
    double *work_bi = (double *) R_alloc((kmax + 1), sizeof(double));
    const double x = t * (lambda_2 - lambda_1);
    const double abs_x = fabs(x);
    const double log_l1 = log(lambda_1), log_l2 = log(lambda_2), log_t = log(t);
    const double max_rate = fmax(lambda_1, lambda_2);
    double log_gamma = 0.0;

    hypgeoF11_k_2k(abs_x, result, work_bi);
    if (x > 0.0) {
        hypgeoF11_kp1_2kp1(abs_x, result, work_bi);
    } else {
        hypgeoF11_k_2kp1(abs_x, result, work_bi);
    }
    result[0] = -lambda_1 * t;
    for (int k = 1; k <= kmax; k++) {
        result[2 * k - 1] = k * log_l1 + (k - 1) * log_l2 - log_gamma -
                max_rate * t + (2 * k - 1.0) * log_t + result[2 * k - 1];
        log_gamma = log(2 * k) + log_gamma;
        if (k < kmax || nmax % 2 == 0) {
            result[2 * k] = k * (log_l1 + log_l2) - log_gamma -
                    max_rate * t + (2 * k) * log_t + result[2 * k];
            log_gamma = log(2 * k + 1.0) + log_gamma;
        }
    }
    return ;
}

RcppExport SEXP hypgeo_test(SEXP rx, SEXP rout) {
    arma::vec x = Rcpp::as<arma::vec>(rx);
    arma::vec out = Rcpp::as<arma::vec>(rout);
    const int nmax = out.size() - 1;
    const int kmax = getMaxK(nmax);
    double * work = (double *) R_alloc((kmax + 1), sizeof(double));
    hypgeoF11_k_2k(x[0], out.col(0), work);
    return Rcpp::wrap(out);
}

RcppExport SEXP pclt_test(SEXP rpar, SEXP rout) {
    arma::vec out = Rcpp::as<arma::vec>(rout);
    NumericVector par(rpar);
    const double t = par[0], l1 = par[1], l2 = par[2];
    pclt(t, l1, l2, out.col(0));
    return Rcpp::wrap(out);
}



