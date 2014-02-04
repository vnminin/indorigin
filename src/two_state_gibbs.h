/*
 * two_state_gibbs.h
 *
 *  Created on: Oct 31, 2013
 *      Author: hoblitz
 */

#ifndef TWO_STATE_GIBBS_H_
#define TWO_STATE_GIBBS_H_

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

arma::mat twoStateRateMatrix(double lambda_01, double lambda_10);

arma::mat twoStateTransProb(double lambda_01, double lambda_10, double time);

int sampleOnce(arma::colvec weights, double rUnif);

NumericVector twoStateUnifSample(arma::mat& rateMatrix, int startState,
        int endState, double elapsedTime, double transProb);

arma::mat PartLikelihoods(const arma::Mat<int> & treeEdges, const IntegerVector & tipStates,
                                  const arma::cube & cubeProbMat);

double TwoStatePhyloLikelihood(arma::Mat<int>& treeEdges, IntegerVector& tipStates,
                               NumericVector& branchLengths, double lambda_01, double lambda_10,
                               NumericVector& rootDist);

double TwoStatePhyloLikelihood(const arma::Mat<int> & treeEdges, const IntegerVector & tipStates,
                               const arma::vec & branchLengths, const double & lambda_01,
                               const double & lambda_10, const arma::vec & armaRootDist);

NumericVector twoStateSufficientStatistics(arma::Mat<int>& treeEdges, IntegerVector& tipStates,
                                           NumericVector& branchLengths, double lambda_01, double lambda_10,
                                           NumericVector& rootDist);

double twoStateCompleteDataLogPosterior(NumericVector suffStat, double lambda_01, double lambda_10,
                                         double prior_alpha_01, double prior_beta_01,
                                         double prior_alpha_10, double prior_beta_10);

NumericVector twoStatePhyloGibbsSampler(IntegerVector treeEdges, IntegerVector cubeDims, NumericMatrix branchLengths,
                                    NumericVector rootDist, IntegerMatrix tipStates, double initial_lambda_01,
                                    double initial_lambda_10, double prior_alpha_01, double prior_beta_01,
                                    double prior_alpha_10, double prior_beta_10, int mcmcSize, int mcmcBurnin,
                                    int mcmcSubsample);
#endif /* TWO_STATE_GIBBS_H_ */
