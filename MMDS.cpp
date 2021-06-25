#include "MMDS.hpp"



#include <lp_lib.h>

#include "blas/f2c.h"
#include "blas/clapack.h"
#include "pswarm.h"

using namespace std;
using namespace arma;

//#include "Armadillo/armadillo_bits/Mat_bones.hpp"

//extern int dgels_(char *trans, integer *m, integer *n, integer *
//	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
//	doublereal *work, integer *lwork, integer *info);

void MMDS::executeMMDS(Problem& problem, MMDSOption& options, arma::mat& X, arma::colvec& F, arma::colvec& Alphas, Stat& stats) {

    //Be sure to have a column vector
    int ULBSize = problem.lb.n_rows;
    if (ULBSize != problem.n) {
        problem.lb = problem.lb.t();
    }
    ULBSize = problem.ub.n_rows;
    if (ULBSize != problem.n) {
        problem.ub = problem.ub.t();
    }

    //%Initialize options
    maxIter = options.maxIter;
    maxIIter = options.maxIIter;
    n = problem.n;
    maxEvals = max(options.maxObj, 2000);
    type = options.type;
    clusterType = options.clusterType;
    alpha_tol = options.alpha_tol;
    phi = options.phi;
    theta = options.theta;
    randInit = options.randInit;
    pSize = options.pSize;
    proc = options.proc;
    iSize = options.iSize;
    maxRuns = options.maxRuns;

    if (iSize > pSize) {
        iSize = pSize;
    }

    lb = problem.lb;
    ub = problem.ub;

    //lb.print("lb: ");
    //ub.print("ub: ");

    //% start clock
    startClock = clock();


    //% Initialize random numbers
    //srand(time(NULL));

    arma_rng::set_seed(1234567891);
    srand(1234567891);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //% Initialization step.
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //% Compute f at the initial point.

    //% Multistart
    //% A Column vector of initial guesses
    if (problem.x0.n_rows == n && problem.x0.n_cols > pSize)
        pSize = problem.x0.n_cols; // history starts with a minimum of n_cols points
    mat history(n, pSize);
    history.fill(0.0);
    
    int start = 0;
    if (problem.x0.n_rows == n) {
        int m = problem.x0.n_cols;
        history(span(0, n - 1), span(0, m - 1)) = projection(problem.x0, problem.lb, problem.ub);
        //p.print("projection: ");
        //

        start = m;
    }
    //std::cout << "start history: " << history;

    double maxBound = 0.0;
    colvec scale(n);
    scale.fill(1.0);


    if (!problem.lb.is_finite() || !problem.ub.is_finite()) {
        for (int i = 0; i < n; i++) {
            if (is_finite(lb(i)) && is_finite(ub(i))) {
                history(i, span(start, history.n_cols - 1)) = (ub(i) - lb(i)) * randu(1, (pSize - start /* +1*/)) + lb(i);
                scale(i) = ub(i) - lb(i);

                if (maxBound < ub(i) - lb(i)) {
                    maxBound = ub(i) - lb(i);
                }
            } else if (is_finite(problem.lb(i))) {
                history(i, span(start, history.n_cols - 1)) = (lb(i)) + randu(1, (pSize - start /*+ 1*/));
            } else if (is_finite(problem.ub(i))) {
                history(i, span(start, history.n_cols - 1)) = (ub(i)) - randu(1, (pSize - start /*+ 1*/));

            }
        }
    } else {
        //% random start with latin hypercubes
        scale = ub - lb;
        mat lhsTemp = lhs(n, pSize - start /*+ 1*/);

        //TODO this is just to get the exact same solution to matlab
        //lhsTemp.reset();
        //lhsTemp << 0.3926 << 0.0503 << 0.7575 << 0.5999 << 0.8296 << 0.3553 << 0.5328 << 0.0132 << 0.6314 << 0.0551 << endr << 0.5299 << 0.7713 << 0.1098 << 0.1979 << 0.1252 << 0.7050 << 0.3839 << 0.9408 << 0.2582 << 0.9030 << endr;

        //lhsTemp.print("lhs: ");
        mat repTemp = repmat((ub - lb), 1, pSize - start /*+ 1*/);
        //repTemp.print("repTemp: ");

        mat multTemp = lhsTemp % repTemp;
        //multTemp.print("multTemp: ");

        mat rrepTemp = repmat(lb, 1, pSize - start /*+ 1*/);
        //rrepTemp.print("rrepTemp:");

        //(rrepTemp + multTemp).print("sum: ");
        history(span(0, history.n_rows - 1), span(start, history.n_cols - 1)) = rrepTemp + multTemp;

        for (int i = 0; i < n; i++) {
            if (is_finite(lb(i)) && is_finite(ub(i))) {
                if (maxBound < ub(i) - lb(i)) {
                    maxBound = ub(i) - lb(i);
                }
            }
        }
    }

    //std::cout << "History:" << history;
    //scale.print("scale: ");

    if (maxBound <= 0) {
        alpha = options.alpha;
    } else {
        alpha = maxBound / (iSize * 500);
    }

    double alpha_squared = 0;
    if (alpha > 1) {
        alpha_squared = alpha * alpha;
    } else if (alpha < 1) {
        alpha_squared = sqrt(alpha);
    } else {
        alpha_squared *= 2;
    }

    //cout << "alpha: " << alpha << endl;
    //cout << "alpha_squared: " << alpha_squared << endl;

    colvec historyF(pSize);
    historyF.fill(0.0); // = zeros<mat>(pSize, 1);

    //% compute objective function in parallel

#pragma omp parallel for
    for (int p = 0; p < pSize; p++) {

        //convert to double*
        //colvec temp = history.col(p);
        //double* temp2 = (double*)malloc(temp.n_rows * sizeof(double));
        //double* temp3;

        //for (int z = 0; z < temp.n_rows; z++) {
        //	temp2[z] = temp(z);
        //}

        double *temp2 = history.colptr(p);
        double *temp3 = historyF.memptr() + p;

        objf(2, 1, temp2, temp3);
        //cout << "insert obj value (historyf) for position " << (p + 1) << ":  ";
        //cin >> historyF(p);
        //historyF(p) = temp3;// temp3[0];

        //free(temp2);
    }

    stats.func_eval += pSize;
    
    //std::cout << "HistoryF: " << historyF;

    uvec idx = sort_index(historyF);
    //cout << "IDX" << idx;
    idx = idx(span(0, (iSize - 1)));
    //cout << "IDX" << idx;

    //% points started a local search
    colvec historyS(pSize);
    historyS.fill(false);

    for (int i = 0; i < idx.n_rows; i++) {
        historyS(idx(i)) = true;
    }

    //historyS.print("HistoryS: ");

    //% Random points
    colvec historyR(pSize);
    historyR.fill(true);

    stats.random += pSize;


    //% Initialize some auxiliary variables.
    //% Each cluster has its own alpha
    colvec alphas(iSize);
    alphas.fill(alpha);

    //% are points active?
    colvec flags(iSize);
    flags.fill(true);

    //alphas.print("alphas: ");
    //flags.print("flags: ");

    //% the number of directions used
    int cardD = 0;
    switch (type) {
        case SearchType::twoNRandom: //% 2n random
            cardD = 2 * n;
            break;
        case SearchType::singleRandom: //% one random
            cardD = 1;
            break;
        case SearchType::twoRandom:// % two random
            cardD = 2;
            break;
        case SearchType::singleRandomSymmetric:// % one random and its symmetric
            cardD = 2;
            break;
        default:
            std::cout << "Invalid local search type" << endl;
    }


    //% adapt Phi to the theory when we have random directions
    //% Corollary B.4. SIOPT
    phi = max(phi, exp(log(theta) / (1 - pow(2, cardD - 1))));
    //cout << "phi: " << phi << endl;


    //% shall we stop?
    bool halt = false;
    while (!halt && sum(flags) > 0) {

        //% Clustering based on distance and quadratic model fitting
        //% we could do this in parallel!

        //% Idx are the indexes for point where local search (polling) is happening
        int k = arma::max(arma::size(idx));
        int h = history.n_cols;

        //% compute distances need for clustering
        mat historyColumns = history.cols(idx);

        //pdist2 for armadillo...
        mat distances = pdist2(history.t(), historyColumns.t());
        
        //std::cout << "Distances: " << distances;


        //% we have a valid model?
        colvec valid_model(k);
        valid_model.fill(false);

        //% model convexity
        colvec convex(k);
        convex.fill(false);
        colvec concave(k);
        concave.fill(false);

        if (h > k) { //% enough points to cluster
            //% build k clusters, using the poll centers as centroids
            uvec idx_k = index_min(distances, 1);
            //idx_k.print("idx_k: ");

            //% some temporary variables
            //% how is the model adjusted to the quadratic
            colvec adjust(k);
            adjust.fill(0.0);

            //% index of maximum distance from poll center, i.e., index of point
            //% where next run should be started if model is not adjusted
            uvec model_max_diff_idx(k);
            model_max_diff_idx.fill(0);


            //% go for each cluster
            for (int i = 0; i < k; i++) {
                //% build quadratic for each cluster
                uvec cluster = find(idx_k == i);
                int e = -1; //exit flag
                colvec ps(2 * n + 1);

                if (clusterType == ClusteringType::QUAD) {
                    //% a underestimate quadratic model
                    //TODO quad_model
                    e = quad_model(history.cols(cluster), historyF(cluster), ps);
                } else {
                    //% a least squares quadratic model
                    e = lsq_model(history.cols(cluster), historyF(cluster), ps);
                }
                if (e >= 0) {
                    //% we have a quadratic model
                    valid_model(i) = true;

                    //% check for model convexity by looking into Hessian
                    //% diagonal values
                    colvec hess_ps = ps(span(n + 1, ps.n_elem - 1));
                    //cout << "Hess diag:" << hess_ps;
                    if (all(hess_ps > 0)) {
                        convex(i) = true;
                    }
                    if (all(hess_ps < 0)) {
                        concave(i) = true;
                    }

                    //% measure how good is the quadratic w.r.t. cluster
                    colvec model_diff = quadratic(history.cols(cluster), ps) - historyF(cluster);
                    adjust(i) = accu(model_diff % model_diff);

                    mat dTemp = distances.rows(cluster);
                    mat cTemp = dTemp.col(i);

                    uvec model_max_diff_idx_tmp = sort_index(cTemp, "descend");

                    //% By setting Idx we are not adding points from history
                    model_max_diff_idx(i) = idx(i);

                    int l=model_max_diff_idx_tmp.size();
                    bool notfoundit=true;
                    for(int j=0; notfoundit && j<l;j++){
                        int z=cluster(model_max_diff_idx_tmp(j));
                        if(distances(z,i)>alpha_squared && historyR(z) && !historyS(z)){
                            model_max_diff_idx(i) = z;
                            notfoundit=false;
                        }
                    }
                    
                    //std::cout << "Far away: " << model_max_diff_idx;
                        
                } else {
                    //% failed to build the quadratic model
                    //% local search will produce new points "close" to the poll
                    //% center
                }
            }

            //% check clusters to split
            if (k < maxRuns) {
                for (int i = 0; i < k; i++) {
                    if (valid_model(i) && adjust(i) > alpha) {
                        //% We have at leat one cluster to split
                        //cout << "We want to split cluster " << i << endl;
                        if (model_max_diff_idx(i) != idx(i)) {
                            uvec temp;
                            temp << model_max_diff_idx(i) << endr;
                            idx = join_vert(idx, temp);
                            flags = join_vert(flags, mat(1, 1, fill::ones));

                            //% assume convexity from the previous model
                            mat convTemp(1, 1);
                            convTemp.fill(convex(i));
                            mat concTemp(1, 1);
                            concTemp.fill(concave(i));

                            convex = join_vert(convex, convTemp);
                            concave = join_vert(concave, concTemp);

                            //% step size heritage
                            mat alphaTemp(1, 1);
                            alphaTemp.fill(alpha);
                            alphas = join_vert(alphas, alphaTemp);
                            //% we are starting a run, so do not use this point
                            //% for a new start
                            historyS(model_max_diff_idx(i)) = true;
                        } else {
                            //% add a random point ...
                            double away_dist = 0;
                            mat newPoint(n, 1, fill::zeros);
                            while (away_dist < alpha_squared) {
                                for (int j = 0; j < n; j++) {
                                    if (is_finite(lb(j)) && is_finite(ub(j))) {
                                        newPoint(j) = (ub(j) - lb(j)) * (double) ((double)rand() / (double) RAND_MAX) + lb(j);
                                    } else if (is_finite(problem.lb(j))) {
                                        newPoint(j) = lb(j) + (double) ((double)rand() / (double) RAND_MAX);
                                    } else if (is_finite(ub(j))) {
                                        newPoint(j) = ub(j) - (double) ((double)rand() / (double) RAND_MAX);
                                    }
                                }
                                away_dist = norm(history.col(idx(i)) - newPoint);
                            }
                            
                            //std::cout << "New random point" << newPoint;

                            //% no need to projet NewPoint to the feasible
                            //% region, so add it to history
                            history = join_horiz(history, newPoint);
                            
                            
                            
                            mat hft(1, 1);
                            objf(2, 1, newPoint.memptr(), hft.memptr());
                            

                            historyF = join_vert(historyF, hft);
                            historyS = join_vert(historyS, mat(1, 1, fill::zeros));
                            historyR = join_vert(historyR, mat(1, 1, fill::ones));
                            stats.random++;
                        }
                    }
                }
            }

            k = max(idx.n_rows, idx.n_cols);
            //cout << "We have " << accu(convex) << " convex model(s) of " << k << " models" << endl;
            //cout << "We have " << accu(concave) << " concave model(s) of " << k << " models" << endl;

        }

        //% Do not poll on points too "close" to each other
        //% points from one cluster may end up in other cluster
        k = max(idx.n_rows, idx.n_cols);
        distances = pdist2(history.cols(idx).t(), history.cols(idx).t());

        colvec beRemoved(k);
        beRemoved.fill(false);
        colvec beKept(k);
        beKept.fill(false);

        for (int i = 0; i < k; i++) {
            if (!beRemoved(i)) {
                for (int j = i + 1; j < k; j++) {
                    if (!beRemoved(j) && !beKept(j) && i != j) {
                        //% check for near poll centers

                        if (distances(j, i) < alpha_squared) {
                            //% j is reacheable by i, so remove j if worst
                            if (historyF(idx(i)) < historyF(idx(j))) {
                                beRemoved(j) = true;
                                beKept(i) = true;
                                //cout << "Removing point " << j << " due to proximity with " << i << endl;
                            } else {
                                beRemoved(i) = true;
                                beKept(j) = false;
                                //cout << "Removing point " << i << " due to proximity with " << j << endl;
                            }
                        }
                    }
                }
            }
        }

        if (any(beRemoved)) {
            flags.shed_rows(find(beRemoved == true));
            alphas.shed_rows(find(beRemoved == true));
            idx.shed_rows(find(beRemoved == true));
            convex.shed_rows(find(beRemoved == true));
            concave.shed_rows(find(beRemoved == true));
            //cout << "Removed " << accu(beRemoved) << " point(s) due to proximity" << endl;
        }

        //% just run for active point
        //uvec t(idx);
        uvec runs = idx((find(flags == true)));
        uvec runsIdx = find(flags == true);
        int nRuns = runs.size();
//        std::cout << "Idx=" << idx;
//        std::cout << "flags=" << flags;
//        std::cout << "Runs="<<runs;
//        std::cout << "RunsIdx="<<runsIdx;


        //% perform inner iterations
        for (int j = 0; nRuns >= 1 && j < maxIIter; j++) {
            //cout << "inner iteration: " << (j + 1) << endl;

            vector<Run> run;
            //%initialize memory for all runs
            for (int i = 0; i < nRuns; i++) {
                run.push_back(
                        Run(
                        mat(n, cardD, fill::zeros),
                        colvec(cardD),
                        history.col(runs(i)),
                        historyF(runs(i)),
                        convex(runsIdx(i))
                        ));
            }
            
//            std::cout << "Run:" << run.at(0).x;
//            std::cout << "Run:" << run.at(1).x;
//            std::cout << "Run:" << run.at(2).x;

            // this could be run in parallel
            for (int tIdx = 0; tIdx < nRuns; tIdx++) {
                double rho = min(1e-5, (1e-5) * pow(alphas(runsIdx(tIdx)), 2));

                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                //% Generation of the polling set
                //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                //% look for active bound constraints w.r.t. the step size
                uvec Ilb = find(run.at(tIdx).x <= lb + sqrt(alpha_tol));
                uvec Iub = find(run.at(tIdx).x >= ub - sqrt(alpha_tol));

                mat D;
                switch (type) {
                    case SearchType::twoNRandom:
                        D = mat(n, n);
                    {
                        double TAU;
                        double *work = (double *) malloc(n * sizeof (double));
                        double *Awork = (double *) malloc(n * n * sizeof (double));
                        integer N = 1, M = n, info, K = 1;

                        memset(work,0,n * sizeof (double));
                        memset(Awork,0,n * n * sizeof (double));
                        
                        for (int j = 0; j < n; j++){ // fortran is column major order
                            *(Awork + j) = (double) ((double)rand() / (double) RAND_MAX);
                            //printf("A%d: %f\n",j, *(Awork + j));
                        }

                        dgeqr2_(&M, &N, Awork, &M, &TAU, work, &info);

                        dorgqr_(&M, &M, &K, Awork, &M, &TAU, work, &M, &info);

                        if (!info) {
                            for (int j = 0; j < n; j++)
                                for (int l = 0; l < n; l++)
                                    D(j, l) = *(Awork + j * n + l);
                        }
                        free(Awork);
                        free(work);
                    }

                        //D.print("D: ");

//                        for (int i = 0; i < D.n_cols; i++) {
//                            D.col(i) = D.col(i) / norm(D.col(i));
//                        }
                        // D and its symmetric
                        D = join_horiz(D, -D);
                        break;
                    case SearchType::singleRandom:
                        D = 2 * randu(n, cardD) - 1;
                        for (int i = 0; i < D.n_cols; i++) {
                            D.col(i) = D.col(i) / norm(D.col(i));
                        }
                        break;
                    case SearchType::twoRandom:
                        D = 2 * randu(n, cardD) - 1;
                        for (int i = 0; i < D.n_cols; i++) {
                            D.col(i) = D.col(i) / norm(D.col(i));
                        }
                        break;
                    case SearchType::singleRandomSymmetric:
                        D = 2 * randu(n, 1) - 1;
                        D = D / norm(D);
                        D = join_horiz(D, -1 * D);
                        break;
                }


                //% positive direction if active at lower bound
                D.rows(Ilb) = abs(D.rows(Ilb));

                //% negative direction if active at upper bound
                D.rows(Iub) = -abs(D.rows(Iub));

                D = D % repmat(scale, 1, D.n_cols);


                while (!run.at(tIdx).sucess && run.at(tIdx).count_pss < cardD) {
                    colvec d = D.col(run.at(tIdx).count_pss);
                    colvec temp = run.at(tIdx).x + alphas(runsIdx(tIdx)) * d;
                    run.at(tIdx).history.col(run.at(tIdx).count_pss) = projection(temp, lb, ub);


                    //% Compute the function value at the point.
                    objf(2, 1, run.at(tIdx).history.colptr(run.at(tIdx).count_pss), run.at(tIdx).historyF.memptr()+run.at(tIdx).count_pss);


                    //% Test for a better function value.
                    if (run.at(tIdx).historyF(run.at(tIdx).count_pss) < run.at(tIdx).f - rho) {
                        run.at(tIdx).sucess = 1;
                    }

                    run.at(tIdx).count_pss++;

                }

            }

            //% Syncronize all information
            for (int tIdx = 0; tIdx < nRuns; tIdx++) {
                //% account for function evaluations
                stats.func_eval += run.at(tIdx).count_pss;

                //% Add to History
                int sizeH = history.n_cols;
                history = join_horiz(history, (run.at(tIdx).history.cols(0, run.at(tIdx).count_pss - 1)));
                historyF = join_vert(historyF, (run.at(tIdx).historyF.rows(0, run.at(tIdx).count_pss - 1)));
                historyS = join_vert(historyS, (mat(run.at(tIdx).count_pss, 1, fill::zeros)));
                historyR = join_vert(historyR, (mat(run.at(tIdx).count_pss, 1, fill::zeros)));

                //history.print("history: ");
                //historyF.print("historyF: ");
                //historyS.print("historyS: ");
                //historyR.print("historyR: ");

                if (run.at(tIdx).sucess) {
                    stats.iter_suc++;
                    alphas(runsIdx(tIdx)) = min(maxBound, phi * alphas(runsIdx(tIdx)));
                    idx(runsIdx(tIdx)) = sizeH + run.at(tIdx).count_pss - 1;
                } else {
                    stats.iter_uns++;
                    alphas(runsIdx(tIdx)) = theta * alphas(runsIdx(tIdx));
                }

                //% Test for ending.
                if (alphas(runsIdx(tIdx)) < alpha_tol) {
                    flags(runsIdx(tIdx)) = false;
                    //cout << "Halting point due to Alpha tolerance" << endl;
                }
                if (stats.func_eval >= maxEvals - cardD + 1 && halt == 0) {
                    halt = 1;
                    cout << "Halting due to Maximum Function Evaluations" << endl;
                }
                if (stats.iter > maxIter && halt == 0) {
                    halt = 1;
                    cout << "Halting due to Maximum iterations" << endl;
                }

                if (!any(flags) && halt == 0) {//% All have reached minima
                    halt = 1;
                    cout << "Halting due to all minima/maxima" << endl;
                }

            }

            runs = idx(find(flags == true));
            runsIdx = find(flags == true);
            nRuns = runs.size();
            stats.iter++;
        }




        //% compute distances between poll centers
        mat ttemp = history.cols(idx).t();
        distances = pdist2(ttemp, ttemp);

        k = max(idx.n_rows, idx.n_cols);

        //% remove local minimizers too close to each other
        beRemoved = colvec(k);
        beRemoved.fill(false);
        beKept = colvec(k);
        beKept.fill(false);
        for (int i = 0; i < k; i++) {
            if (!beRemoved(i)) {
                for (int j = i + 1; j < k; j++) {
                    if (!beRemoved(j) && !beKept(j) && i != j) {
                        //% check for near poll centers
                        if (distances(j, i) < alpha_squared) {
                            //% j is reacheable by i, so remove j if worst
                            if (historyF(idx(i)) < historyF(idx(j))) {
                                beRemoved(j) = true;
                                beKept(i) = true;
                                //cout << "Removing point " << j << " due to proximity with " << i << endl;
                                break;
                            } else {
                                beRemoved(i) = true;
                                beKept(j) = true;
                                //cout << "Removing point " << i << " due to proximity with " << j << endl;
                            }
                        }
                    }
                }
            }
        }


        if (any(beRemoved)) {
            uvec idxRemove = find(beRemoved == true);
            flags.shed_rows(idxRemove);
            alphas.shed_rows(idxRemove);
            idx.shed_rows(idxRemove);
            convex.shed_rows(idxRemove);
            concave.shed_rows(idxRemove);
            //cout << "Removed " << idxRemove.size() << " poll point(s) after polling" << endl;
        }
        
        //std::cout << "History:" << history;


    }




    finishClock = clock();

    processing_time = finishClock - startClock;
    //printf ("It took me %d clicks (%f seconds).\n",processing_time,((float)processing_time)/CLOCKS_PER_SEC);
    double totalProcessingTime = processing_time / CLOCKS_PER_SEC;
    cout << "total time: " << totalProcessingTime << " seconds" << endl;
    //cout << "total number of cycles: " << (int)processing_time << endl;

    stats.time = totalProcessingTime;

    stats.f = min(historyF(idx));
    X = history.cols(idx);
    F = historyF(idx);

    std::cout << "Flags " << flags << "Max fun eval " << stats.func_eval << endl;
    
    /*Stats.f=min(Historyf(Idx));
    X = History(:,Idx);
    F = Historyf(Idx);*/

}

mat MMDS::projection(mat x0, colvec lb, colvec ub) {
    mat X = x0;

    for (int i = 0; i < X.n_rows; i++) {
        for (int j = 0; j < X.n_cols; j++) {
            X(i, j) = min(max(X(i, j), lb(i)), ub(i));
        }
    }

    return X;
}

mat MMDS::lhs(int n, int k) {
    bool bPreserveDraw = false;
    bclib::matrixlhs<double> result = bclib::matrixlhs<double>(n, k);

    bclib::CRandomStandardUniform oRandom = bclib::CRandomStandardUniform();
    oRandom.setSeed(1234567891, 1234567891);
    lhslib::randomLHS(n, k, bPreserveDraw, result, oRandom);

    mat r(n, k);
    for (int i = 0; i < result.rowsize(); i++) {
        for (int j = 0; j < result.colsize(); j++) {
            r(i, j) = result(i, j);
        }
    }

    return r;
}

mat MMDS::pdist2(mat X, mat Y) {
    int N = X.n_rows, M = X.n_cols;
    int K = Y.n_rows, M2 = Y.n_cols;

    //TODO make sure that M == M2 otherwise error!
    //matlab equivalent to pdist2: sqrt(sum(X.^2,2)*ones(1,K) + ones(N,1)*sum( Y.^2, 2 )' - 2.*X*Y')

    mat ret = sqrt(sum(X % X, 1) * ones(1, K) + ones(N, 1) * sum(Y % Y, 1).t() - 2 * X * Y.t());

    return ret;
}

colvec MMDS::quadratic(mat x, colvec p) {
    int n = x.n_rows;
    int m = x.n_cols;

    double alpha = p(0);
    colvec c = p.subvec(1, n);
    colvec h = p.subvec(n + 1, p.n_elem - 1);

    colvec q(m);
    q.fill(0.0);
    
    for (int i = 0; i < m; i++){
        mat temp = c.t() * x.col(i) + 0.5 * (x.col(i) % x.col(i)).t() * h;
        q(i) = alpha + temp(0);
    }
    //cout << "q: "<<q;
    return q;
}

int MMDS::lsq_model(mat history, mat historyf, colvec& ps) {
    //% return error by default
    int e = -1;

    //% get problem dimension from data
    int n = history.n_rows;
    int m = history.n_cols;

    if (m < 2 * n + 1) {
        //% we have an underestimation model, so quit
        return e;
    }

    mat A(m, 1, fill::ones);
    //cout << "history" << history;
    //cout << "historyf" << historyf;
    A = join_horiz(A, history.t());
    A = join_horiz(A, 0.5 * ((history % history).t()));
    //ps = solve(A, historyf);
    double *work = (double *) malloc(2 * m * n * sizeof (double));
    double *Awork = (double *) malloc(m*(2*n+1)*sizeof(double));
    double *Bwork= (double *) malloc(m*sizeof(double));
    char type = 'N';
    integer NRHS = 1, N2 = 2*n+1, M = m, MN = 2 * m*n, info;
    //cout << "A:" << A;
    //cout << "b:" << historyf;
    
    memset(work,0,2 * m * n * sizeof (double));
    
    for (int i=0;i<m;i++){
        *(Bwork+i)=historyf(i);
        for (int j=0;j<2*n+1;j++)
            *(Awork+i+j*m)=A(i,j);
    }
    dgels_(&type, &M, &N2, &NRHS, Awork, &M, Bwork, &M, work, &MN, &info);
    
    if (!info) {
        for (int i=0;i<2*n+1;i++)
            ps(i) = *(Bwork+i);
        //cout << "ps"<<ps;
        e = 0;
    }
    
    free(work);
    free(Awork);
    free(Bwork);
    return e;
}

int MMDS::quad_model(mat history, colvec historyf, vec& ps) {
    //% return error by default
    int e = -1;
    lprec *lp=NULL;
        

    //% get problem dimension from data
    int n = history.n_rows;
    int m = history.n_cols;

    if (m < 2 * n + 1) {
        //% at least 2n+1 point to form the model
        return e;
    }
    //% Hessian is diagonal, so number of variables is 2n+1
    int psize = 2 * n + 1;
         
    //% m linear constraints for the set S
    // allocate memory and set the alpha part
    mat A(m, psize, fill::ones);

    mat temp(m, n, fill::zeros);
    rowvec tempH(n);
    tempH.fill(0.0);

    for (int j = 0; j < m; j++) {
        //% H's coeficients
        temp.row(j) = 0.5 * (history.col(j) % history.col(j)).t();
        //Temp(j,:)=0.5*(History(:,j).^2);
        tempH = tempH + temp.row(j);
    }

    //% alphas
    //A(span(0, m - 1), 0) = mat(m, 1, fill::ones);

    //% c's
    A(span(0, m - 1), span(1, n)) = history.t();

    //% Hessian
    A(span(0, m - 1), span(n + 1, psize - 1)) = temp;

    //% Objective function
    //% alphas
    rowvec d(1);
    d(0) = (int) (-1) * m;

    //% c's
    colvec temppp = (-sum(history, 1));
    d = join_horiz(d, temppp.t());

    //% Hessian
    d = join_horiz(d, -tempH);


    //A.print("A: ");
    //b.print("b: ");
    //d.print("d: ");
    //% solve linear program to get s_hat

    /* We will build the model row by row
     So we start with creating a model with 0 rows and psize columns */
    lp = make_lp(0, psize);
    if(!lp)
        return(e); /* couldn't construct a new model... */
    
//    set_col_name(lp, 1, "alpha");
//    set_col_name(lp, 2, "g1");
//    set_col_name(lp, 3, "g2");
//    set_col_name(lp, 4, "h1");
//    set_col_name(lp, 5, "h2");
    
    set_minim(lp); // a minimization problem

    for(int j=0;j<psize;j++){
        if(!set_unbounded(lp, j+1)){
            delete_lp(lp);
            return(e);
        };
    }

    set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */ 
    // first the objective function
    if(!set_obj_fn(lp, d.memptr()-1)){
        delete_lp(lp);
        return(e);
    }
    
    // and then the constraints
    for(int j=0;j<m;j++){
        rowvec tmp=A.row(j);
        if(!add_constraint(lp, tmp.memptr()-1, LE, historyf(j))){
            delete_lp(lp);
            return(e);
        }
    }
    set_add_rowmode(lp, FALSE);  /* makes building the model faster if it is done rows by row */

 
    /* I only want to see important messages on screen while solving */
    set_verbose(lp, IMPORTANT);
    
    //write_lp(lp, "mmds.lp");
    
    //set_outputfile(lp, "pois.txt");
    
    /* Now let lpsolve calculate a solution */
    e = solve(lp);
    if(e == OPTIMAL)
      e = 0;
    else
      e = -1;
    
    /* a solution is calculated, now lets get some results */
    /* objective value */
    //printf("Objective value: %f\n", (double)get_objective(lp));
    
    //print_objective(lp);

    /* variable values */
    get_variables(lp,ps.memptr());
//    for (int j = 0; j < psize; j++)
//        printf("%s: %f\n", get_col_name(lp, j + 1), ps(j));
    
   
    delete_lp(lp);
    return(e);
    
}