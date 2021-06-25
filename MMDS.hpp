#ifndef MMDS_hpp
#define MMDS_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "Armadillo\\armadillo"
#include <type_traits>
#include "lhs\\utilityLHS.h"
#include "lhs\\CRandom.h"
#include "lhs\\matrixlhs.h"


//using namespace std;
//using namespace arma;

enum class ClusteringType {
    kmeans = 0, /**< @brief kmeans.
                 * Removes one poll point from clusters with more than one poll point.
                 * Adds a new poll point for each cluster without a poll point. */
    notUsed = 1, /**< @brief NOT USED. */
    FSLA = 2, /**< @brief FSLA. Not tested. */
    QUAD = 4,
    LSQ = 5
};

template<typename T>
std::ostream& operator<<(typename std::enable_if<std::is_enum<T>::value, std::ostream>::type& stream, const T& e) {
    return stream << static_cast<typename std::underlying_type<T>::type> (e);
}

enum class SearchType {
    twoNRandom = 0, /**< @brief 2n random directions */
    singleRandom = 1, /**< @brief a single random direction */
    twoRandom = 2, /**< @brief two random directions */
    singleRandomSymmetric = 3 /**< @brief a single random direction and its symmetric */
};



class MMDSOption {
public:
    int maxObj; /**< @brief maximum number of objective function evaluations (_default_ **max(200*n,2000)**) */
    int maxIter; /**< @brief maximum number of iterations allowed */
    int maxIIter; /**< @brief maximum number of inner iterations allowed */
    SearchType type; /**< @brief selects type of search directions used (_default_ **3**) */
    int variables; /**< @brief number of problem variables */
    ClusteringType clusterType; /**< @brief Type of clustering used (_default_ **2**) */
    double alpha; /**< @brief Initial alpha for each point in the population */
    double alpha_tol; /**< @brief Alpha tolerance. Stopping criteria */
    double phi; /**< @brief Alpha factor on success >=1 (_default_ **1**).
                 * @details Used Phi is computed to comply with convergence theory `max(Phi,exp(log(Theta)/(1-2^(cardD-1))))`. */
    double theta; /**< @brief Alpha factor on insuccess <1 (_default_ **0.5**) */
    int pSize;
    int iSize;
    int proc;
    int randInit;
    int maxRuns;

    MMDSOption() {
        //cout << "Options create" << endl;
        this->maxObj = 10000;
        this->maxIter = 80000;
        this->maxIIter = 10;
        this->alpha = 1;
        this->alpha_tol = 1e-5;
        this->phi = 1;
        this->theta = 0.5;
        this->randInit = 1;
        this->type = SearchType::twoNRandom; //SearchType::twoRandom;
        this->pSize = 10;
        this->clusterType = ClusteringType::QUAD; //ClusteringType::QUAD ou LSQ
        this->proc = 0;
        this->iSize = 2;
        this->maxRuns = 20;
        this->variables = -1; //this means that is has not been initialized
    };

    MMDSOption(const MMDSOption& orig) {
        //cout << "Options clone" << endl;
        this->maxObj = orig.maxObj;
        this->maxIter = orig.maxIter;
        this->maxIIter = orig.maxIIter;
        this->alpha = orig.alpha;
        this->alpha_tol = orig.alpha_tol;
        this->phi = orig.phi;
        this->theta = orig.theta;
        this->randInit = orig.randInit;
        this->type = orig.type;
        this->pSize = orig.pSize;
        this->clusterType = orig.clusterType;
        this->proc = orig.proc;
        this->iSize = orig.iSize;
        this->maxRuns = orig.maxRuns;
        this->variables = orig.variables;
    };

    virtual ~MMDSOption() {
        //cout << "Going to free memory for Options" << endl;
    };

    void print() {
        std::cout << "Options values: " << std::endl;
        std::cout << "Max Obj: " << maxObj << std::endl;
        std::cout << "Max iter: " << maxIIter << std::endl;
        std::cout << "Max iiter" << this->maxIIter << std::endl;
        std::cout << "Alpha: " << this->alpha << std::endl;
        std::cout << "Alpha tol: " << this->alpha_tol << std::endl;
        std::cout << "Phi: " << this->phi << std::endl;
        std::cout << "Theta: " << this->theta << std::endl;
        std::cout << "Rand init: " << this->randInit << std::endl;
        std::cout << "Search type: " << this->type << std::endl;
        std::cout << "pSize: " << this->pSize << std::endl;
        std::cout << "Cluster Type: " << this->clusterType << std::endl;
        std::cout << "Proc: " << this->proc << std::endl;
        std::cout << "iSize: " << this->iSize << std::endl;
        std::cout << "Max Runs: " << this->maxRuns << std::endl;
        std::cout << "Variables: " << this->variables << std::endl;
    };
};


#pragma once

class Problem {
public:
    arma::colvec lb;
    arma::colvec ub;
    arma::mat x0;
    int n;


    Problem(arma::colvec lb, arma::colvec ub, arma::mat x0) {
        this->lb = lb;
        this->ub = ub;
        this->x0 = x0;
        this->n = max(arma::size(lb));
    };

    Problem(arma::colvec lb, arma::colvec ub) {
        this->lb = lb;
        this->ub = ub;
        this->x0 = arma::mat(0, 0);
        this->n = max(arma::size(lb));
    };

    Problem(const Problem& orig) {
        this->lb = orig.lb;
        this->ub = orig.ub;
        this->x0 = orig.x0;
        this->n = orig.n;
        //cout << "Problem clone" << endl;
    };

    virtual ~Problem() {
        //cout << "Problem Destroy" << endl;
    };

    void print() {
        std::cout << "Problem features: " << std::endl;
        lb.print("LB: ");
        ub.print("UB: ");
        x0.print("x0: ");
        std::cout << "n: " << n << std::endl;
    };

};



class Stat {
public:
    int func_eval; /**< @brief Function evaluation counter. */
    int iter_suc; /**< @brief Successful iteration counter. */
    int iter_uns; /**< @brief Unsuccessful iteration counter. */
    int iter; /**< @brief Iteration counter. */
    int random; /**< @brief number of random points */
    double time;
    double f;

    Stat(){
        func_eval = 0;
        iter_suc = 0;
        iter_uns = 0;
        iter = 0;
        random = 0;
        f = 0;
        time = 0;
        //cout << "Stats create" << endl;
    };
    Stat(const Stat& orig){
        func_eval = orig.func_eval;
        iter_suc = orig.iter_suc;
        iter_uns = orig.iter_uns;
        iter = orig.iter;
        random = orig.random;
        f = orig.f;
        time = orig.time;
        //cout << "Stats Cloning" << endl
    };
    virtual ~Stat(){
        //cout << "Going to free memory for Stats" << endl;
    };
};


class MMDS {
public:
    int maxIter;
    int maxIIter;
    int n;
    int maxEvals;
    SearchType type;
    ClusteringType clusterType;
    double alpha_tol;
    double alpha;
    double phi;
    double theta;
    int randInit;
    int pSize;
    int proc;
    int iSize;
    int maxRuns;

    arma::mat lb, ub;

    void (*objf)(int n, int m, double* vectorx, double* vectorfx); //objective function

    clock_t startClock, finishClock, processing_time;

    void executeMMDS(Problem& p, MMDSOption& o, arma::mat& X, arma::colvec& F, arma::colvec& Alphas, Stat& stats);
    arma::mat projection(arma::mat x0, arma::colvec lb, arma::colvec ub);
    arma::mat lhs(int n, int k);
    arma::mat pdist2(arma::mat X, arma::mat Y);
    arma::vec quadratic(arma::mat x, arma::vec p);
    int lsq_model(arma::mat history, arma::mat historyf, arma::vec & ps);
    int quad_model(arma::mat history, arma::colvec historyf, arma::vec& ps);
private:

};

class Run {
public:

    Run(arma::mat history, arma::colvec historyF, arma::mat x, double f, int convex) {
        this->history = history;
        this->historyF = historyF;
        this->x = x;
        this->f = f;
        this->convex = convex;
        this->count_pss = 0;
        this->sucess = 0;
    };

    virtual ~Run(){
        // free
    };

    arma::mat history;
    arma::colvec historyF;
    int count_pss;
    int sucess;
    arma::mat x;
    double f;
    int convex;

};


#endif