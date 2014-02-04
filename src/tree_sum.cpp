/*
 * tree_sum.cpp
 *
 *  Created on: Nov 8, 2013
 *      Author: gimholte
 */

#include "tree_sum.h"

/* initialize the node jump distributions depending whether we
 * compute the prior or posterior. Currently, missing values are
 * supported for the posterior tip labels but this functionality is
 * not implemented in R, nor in the posterior likelihood calculator.
 */
void PhyConvolver::initializeTransdist(const PhyTree & tree,
        const std::string & dist_type) {
    /* choose between prior and posterior initialization patterns */
    if (dist_type == "prior") {
        transdist_0(arma::span(0, 0),
                arma::span(0, tree.getNumTips() - 1)).fill(1);
        transdist_1(arma::span(0, 0),
                arma::span(0, tree.getNumTips() - 1)).fill(1);
    }
    else {
        int state;
        for (int i = 0; i < tree.getNumTips(); i++) {
            state = tree.getTipState(i);
            if (state == 0)
                transdist_0(0, i) = 1;
            else if (state == 1)
                transdist_1(0, i) = 1;
            else if (state == -1) {
                transdist_0(0, i) = 1;
                transdist_1(0, i) = 1;
            } else
                Rcpp::stop("Invalid tip state in initialization.");
        }
    }
    return;
}

/* For a left and right branch index, this function convolves
 * their child node distribution and branch transition probabilities
 * combining to "n" jumps from a parent node in state "base_state".
 */
double PhyConvolver::convolveBelowNode(const int & left_br_idx,
        const int & right_br_idx, const int & n,
        const int & base_state, const PhyTree & tree) const
{
    const int left_ch_node = tree.getChildNode(left_br_idx);
    const int right_ch_node = tree.getChildNode(right_br_idx);
    const int par_node = tree.getParentNode(left_br_idx);
    if (par_node != tree.getParentNode(right_br_idx)) {
        Rcpp::stop("left/right branches must have same parent");
    }

    double sum_left(0.0);
    double sum_right(0.0);
    double sum_total(0.0);
    int j, k;

    for(int i = 0; i <= n; i++) {
        /* sum over combinations of left branch terminating
         * in state 0, starting in state "base" */
        for(j = 0, k = i; j <= i; j++, k--) {
            sum_left += transdist_0(j, left_ch_node) *
                    branchTransitionProb(k, left_br_idx, base_state, 0);
        }

        /* sum over combinations of left branch terminating in
         * state 1, starting in state "base" */
        for(j = i, k = 0; j >= 0; j--, k++) {
            sum_left += transdist_1(j, left_ch_node) *
                    branchTransitionProb(k, left_br_idx, base_state, 1);
        }

        /* right branch terminating in state 0, starting in state "base" */
        for(j = 0, k = n - i; j <= n - i; j++, k--) {
            sum_right += transdist_0(j, right_ch_node) *
                    branchTransitionProb(k, right_br_idx, base_state, 0);
        }

        /* right branch terminating in state 1, starting in state "base" */
        for(j = n - i, k = 0; j >= 0; j--, k++) {
            sum_right += transdist_1(j, right_ch_node) *
                    branchTransitionProb(k, right_br_idx, base_state, 1);
        }
        sum_total += sum_left * sum_right;
        sum_left = 0.0;
        sum_right = 0.0;
    }
    return sum_total;
}

/* divide a node's transition distribution by its maximum
 * probability value, to help prevent underflow
 */
void PhyConvolver::updateRescaleValue(const int & par_node_idx)
{
    const double m0 = arma::max(transdist_0.col(par_node_idx));
    const double m1 = arma::max(transdist_1.col(par_node_idx));
    const double m = std::max(m0, m1);
    if (m <= 0.0)
        return;
    transdist_0.col(par_node_idx) /= m;
    transdist_1.col(par_node_idx) /= m;
    rescale += log(m);
}

/* Fill the transition distribution for a parent node,
 * as identified via a post order index value.
 */
void PhyConvolver::fillNodeProbs(const int & post_order_idx,
        const PhyTree & tree)
{
    const int left_br_idx = post_order_idx;
    const int right_br_idx = post_order_idx + 1;
    const int par_node_idx = tree.getParentNode(post_order_idx);

    for(int n = 0; n <= n_max; n++) {
        transdist_0(n, par_node_idx) =
                convolveBelowNode(left_br_idx, right_br_idx,
                        n, 0, tree);

        transdist_1(n, par_node_idx) =
                convolveBelowNode(left_br_idx, right_br_idx,
                        n, 1, tree);
    }

    updateRescaleValue(par_node_idx);
    return;
}

/* "Q probs" are the probabilites of n jumps in time t, from
 * a chain state. We use the pclt function to compute these.
 * This function fills the q_probs_0/1 members with their proper
 * values according to branch length and initial state.
 */
void PhyConvolver::fillQProbs(const PhyTree & tree, const double & rate_0_to_1,
        const double & rate_1_to_0) {
    double cur_branch_len, prev_branch_len = -1.0;
    for(int edge_idx = 0; edge_idx < tree.getNumEdges(); edge_idx++) {
        cur_branch_len = tree.getBranchLength(edge_idx);
        if (cur_branch_len == prev_branch_len) {
            // Branch lengths are same. No need to re-compute, just copy.
            q_probs_0.col(edge_idx) = q_probs_0.col(edge_idx - 1);
            q_probs_1.col(edge_idx) = q_probs_1.col(edge_idx - 1);
        }
        else {
            pclt(cur_branch_len, rate_0_to_1, rate_1_to_0,
                    q_probs_0.col(edge_idx));
            pclt(cur_branch_len, rate_1_to_0, rate_0_to_1,
                    q_probs_1.col(edge_idx));
            q_probs_0.col(edge_idx).transform(ip_exp());
            q_probs_1.col(edge_idx).transform(ip_exp());
        }
        prev_branch_len = cur_branch_len;
    }
    return;
}

/* Primary function for computing the gains/losses of a trait over
 * an entire phylogeny.
 * INPUTS:
 * tree: a PhyTree object
 * n_max: max n for which to compute probability
 * rate_0_to_1, rate_1_to_0: rates of transition between states
 * root_node_prob_0: probability root node is in state 0
 * dist_type: either "prior" or "posterior" is computed
 * gains: if true, gains (0->1) distribution is computed, else losses (1->0)
 *  distribution is computed
 *
 * OUTPUTS:
 * arma::mat of dimension (n_max + 1, 1), whose k-th element
 * is the (prior/posterior) probability of k (gains/losses) over
 * the given phylogeny at the given rates and root distribution.
 */
arma::mat convolveTree(const PhyTree & tree, const int & n_max,
        const double & rate_0_to_1, const double & rate_1_to_0,
        const double & root_node_prob_0, const std::string & dist_type,
        const bool gains) {
    PhyConvolver conv(n_max, tree.getNumNodes(), gains);
    conv.initializeTransdist(tree, dist_type);
    conv.fillQProbs(tree, rate_0_to_1, rate_1_to_0);
    for(int i = 0; i < tree.getNumEdges(); i += 2) {
        conv.fillNodeProbs(i, tree);
    }

    // root node is now a rescaled probability
    const int root_node_idx = tree.getNumTips();
    IntegerVector tip_states_copy(tree.getNumTips());
    for(int i = 0; i < tree.getNumTips(); i++) {
        tip_states_copy[i] = tree.getTipState(i);
    }

    arma::vec arma_root_dist(2);
    arma_root_dist(0) = root_node_prob_0;
    arma_root_dist(1) = 1.0 - root_node_prob_0;

    // get likelihood of tip data if we are computing posterior
    if (dist_type == "posterior") {
        const double post_lik = TwoStatePhyloLikelihood(tree.getEdgeMatrix() + 1,
                tip_states_copy, tree.getBranchLengths(), rate_0_to_1,
                rate_1_to_0, arma_root_dist);
        conv.rescale -= log(post_lik);
    }
    // gather root node vectors
    arma::mat root_dists = conv.getRootDists(root_node_idx);

    // transform output to log scale, rescale, and combine
    root_dists.transform(ip_log());
    root_dists.col(0) += conv.rescale + log(root_node_prob_0);
    root_dists.col(1) += conv.rescale + log(1.0 - root_node_prob_0);

    arma::vec output(n_max + 1);
    for(int i = 0; i < output.n_rows; i++){
        output(i) = logspaceAdd(root_dists(i, 0), root_dists(i, 1));
    }

    return output;
}

RcppExport SEXP treeConvolve(SEXP r_edge_matrix, SEXP r_branch_lengths,
        SEXP r_tip_states, SEXP r_rate_0, SEXP r_rate_1, SEXP r_n_max, SEXP r_root_p_0,
        SEXP r_dist_type, SEXP r_gains)
{
    PhyTree tree(r_edge_matrix, r_branch_lengths, r_tip_states);
    const int n_max = Rcpp::as<int>(r_n_max);
    const double rate_0_to_1 = Rcpp::as<double>(r_rate_0);
    const double rate_1_to_0 = Rcpp::as<double>(r_rate_1);
    const double root_p_0 = Rcpp::as<double>(r_root_p_0);
    const std::string dist_type = Rcpp::as<std::string>(r_dist_type);
    bool gains = Rcpp::as<int>(r_gains);

    arma::mat output = convolveTree(tree, n_max, rate_0_to_1,
            rate_1_to_0, root_p_0, dist_type, gains);

    return wrap(output);
}

// [[Rcpp::export]]
arma::mat AverageTreesConvolve(IntegerVector r_edge_cube, IntegerVector cubeDims, 
                               arma::mat arma_branch_lengths, arma::imat arma_tip_states, 
                               double rate_0_to_1, double rate_1_to_0, int n_max, double root_p_0,
                               std::string dist_type, int gains){
  
  // Make a cube of tree edge matrices
  arma::Cube<int> cubeTreeEdges(r_edge_cube.begin(), cubeDims[0], cubeDims[1], cubeDims[2], false);


  int numTrees = arma_branch_lengths.n_cols; 

  PhyTree tree(cubeTreeEdges.slice(0), arma_branch_lengths.col(0), arma_tip_states.col(0));

  arma::mat output = convolveTree(tree, n_max, rate_0_to_1,
            rate_1_to_0, root_p_0, dist_type, gains);    
  

    return output;
}


