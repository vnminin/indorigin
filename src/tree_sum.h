/*
 * tree_sum_v2.h
 *
 *  Created on: Nov 8, 2013
 *      Author: gimholte
 */

#ifndef TREE_SUM_V2_H_
#define TREE_SUM_V2_H_

#include "pclt.h"
#include "two_state_gibbs.h"
using namespace Rcpp;

// in-place functors for arma::transform()
struct ip_exp {
    ip_exp() {};
    double operator()(const double x) {
        return std::exp(x);
    };
};

struct ip_log {
    ip_log() {};
    double operator()(const double x) {
        return std::log(x);
    };
};

/* PhyTree class doesn't do much, other than access tree elements
 * in an orderly fashion. The constructor takes three
 * arguments:
 * r_edge_matrix: edge index matrix as from a "phylo" object in
 *  post-order order (aka pruningwise)
 * r_branch_lengths: branch lengths from "phylo" object
 * r_tip_states: vector of 0/1 tip state integer values.
 *
 * The class's public interface is essentially immutable
 * and consists only of getters that help make code
 * more readable and the tree easier to work with.
 */
class PhyTree {
private:
    const arma::vec branch_lengths;
    const arma::ivec tip_states;
    const arma::imat edge_matrix;
    const int num_edges, num_tips;
public:
    PhyTree(SEXP r_edge_matrix, SEXP r_branch_lengths, SEXP r_tip_states) :
        branch_lengths(as<arma::vec>(r_branch_lengths)),
        tip_states(as<arma::ivec>(r_tip_states)),
        edge_matrix(as<arma::imat>(r_edge_matrix) - 1),
        num_edges(branch_lengths.n_elem),
        num_tips(tip_states.n_elem) {};
    PhyTree(arma::imat arma_edge_matrix, arma::vec arma_branch_lengths, arma::ivec arma_tip_states) :
        branch_lengths(arma_branch_lengths),
        tip_states(arma_tip_states),
        edge_matrix(arma_edge_matrix - 1),
        num_edges(branch_lengths.n_elem),
        num_tips(tip_states.n_elem) {};
    const arma::vec & getBranchLengths() const {
        return branch_lengths;
    }
    int getParentNode(const int & branch_idx) const {
        return edge_matrix(branch_idx, 0);
    }
    int getChildNode(const int & branch_idx) const {
        return edge_matrix(branch_idx, 1);
    }
    double getBranchLength(const int & idx) const {
        return branch_lengths(idx);
    }
    int getNumEdges() const {
        return num_edges;
    }
    int getNumNodes() const {
        return num_tips + num_edges / 2;
    }
    int getNumTips() const {
        return num_tips;
    }
    const arma::imat & getEdgeMatrix() const {
        return edge_matrix;
    }
    const arma::ivec & getTipStates() const {
        return tip_states;
    }
    int getTipState(const int & tip_idx) const {
        return tip_states(tip_idx);
    }
};

arma::mat convolveTree(const double & rate_0_to_1,
        const double & rate_1_to_0, const double & base_rate,
        const double & next_rate, const double & root_node_prob_0,
        const PhyTree & tree);

/* The PhyConvolver class manages objects used during PhyTree convolution.
 * In general, a single PhyConvolver should be constructed for each
 * tree convolution performed. The interface does not yet allow a PhyConvolver
 * object to be "reset" for use with, e.g., a new set of parameters.
 * This could be achieved with a simple new member function that basically
 * sets everything to zero again.
 */

class PhyConvolver {
private:
    int n_max;
    bool gains;
    arma::mat transdist_0, transdist_1;
    arma::mat q_probs_0, q_probs_1;
public:
    double rescale;

    const arma::mat & getTransdist(int id) const {
        if (id == 0)
            return transdist_0;
        return transdist_1;
    }

    const arma::mat & getQprobs(int id) const {
        if (id == 0)
            return q_probs_0;
        return q_probs_1;
    }

    arma::mat getRootDists(const int & root_node_idx) const {
        return arma::join_rows(transdist_0.col(root_node_idx),
                            transdist_1.col(root_node_idx));
    }

    void updateRescaleValue(const int & par_node_idx);

    double branchTransitionProb(const int & k, const int & branch_idx,
            const int base_state, const int end_state)
    const {
        const int sign = 2 * ((int) gains) - 1;
        const double d = (base_state - end_state) * sign;
        if (base_state == 0) {
            if (2 * k + d < 0)
                return 0.0;
            return q_probs_0(2 * k + d, branch_idx);
        } else {
            if (2 * k + d < 0)
                return 0.0;
            return q_probs_1(2 * k + d, branch_idx);
        }
    };

    double convolveBelowNode(const int & left_br_idx,
            const int & right_br_idx, const int & n,
            const int & base_state, const PhyTree & tree) const;

    void fillQProbs(const PhyTree & tree, const double & rate_0_to_1,
            const double & rate_1_to_zero);

    void fillNodeProbs(const int & post_order_idx, const PhyTree & tree);

    void initializeTransdist(const PhyTree & tree,
            const std::string & dist_type);

    PhyConvolver(const int & n_max, const int & num_nodes, const bool & gains) :
        n_max(n_max),
        gains(gains),
        rescale(0.0),
        transdist_0(n_max + 1, num_nodes, arma::fill::zeros),
        transdist_1(n_max + 1, num_nodes, arma::fill::zeros),
        q_probs_0(2 * (n_max + 1) + 2, num_nodes, arma::fill::zeros),
        q_probs_1(2 * (n_max + 1) + 2, num_nodes, arma::fill::zeros) {};
};
#endif /* TREE_SUM_V2_H_ */
