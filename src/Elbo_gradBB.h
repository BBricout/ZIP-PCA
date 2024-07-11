#pragma once

#include "RcppArmadillo.h"
#include <cmath>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"
#include "utilsBB.h"


//--------------------------------------------------------------------------------------------------------------------
// Calcul de l'Elbo et des gradients


std::tuple<
    arma::mat, double, double, double, double, double, double,
    arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat
>
Elbo_grad(const arma::mat & Y, const arma::mat & X, const arma::mat & R,
              const arma::mat & B, const arma::mat & D, const arma::mat & C, 
              const arma::mat & M, const arma::mat & S,
              double tolXi) {

    int n = Y.n_rows;
    int p = Y.n_cols;
    int q = M.n_cols;
    arma::vec XB = X * B;
    arma::vec XD = X * D;
    
    arma::mat mu = arma::mat(XB.memptr(), n, p, false, false);
    arma::mat nu = arma::mat(XD.memptr(), n, p, false, false);
    arma::vec vecY = arma::vectorise(Y);
    arma::vec vecR = arma::vectorise(R);
    arma::mat Z = mu + M * C.t();
    arma::vec vecZ = arma::vectorise(Z);
    arma::mat A = exp(Z + 0.5 * S * (C % C).t());
    arma::vec vecA = vectorise(A);
    arma::mat log_fact_Y = log_factorial_matrix(Y);
    arma::mat pi = 1./(1. + exp(-nu));
    arma::vec vecpi = vectorise(pi);
    arma::mat xi = ifelse_mat(Y, A, nu, R, tolXi);
    arma::vec vecxi = vectorise(xi);
    
   
    
    double objective = (accu(xi % nu - ifelse_exp(nu)) + 
                        accu(R % xi % (Y % (mu + M * C.t()) - A - log_fact_Y)) - 
                        0.5 * accu(M % M + S - log(S)) + 
                        entropie_logis(xi) + 0.5 * n * q);
                        
    arma::mat gradB = X.t() * (vecR % vecxi % (vecY - vecA));
    arma::mat gradD = X.t() * (vecR % (vecxi - vecpi));
    arma::mat gradC = (R % xi % (Y - A)).t() * M - (R % xi % A).t() * S % C;
    arma::mat gradM = (R % xi % (Y - A) * C - M);
    arma::mat gradS = 0.5 * (1. / S - 1. - R % xi % A * (C % C));
    
     
    
    double elbo1 = accu(xi % nu - ifelse_exp(nu));
    double elbo2 = - 0.5 * accu(M % M + S);
    double elbo3 = accu(R % xi % (Y % (mu + M * C.t()) - A - log_fact_Y));
    double elbo4 = entropie_logis(xi);
    double elbo5 = 0.5 * accu(0.5 * log(S % S)) + n * q * 0.5;
    
   
    
    return std::make_tuple(
        xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective,
        gradB, gradD, gradC, gradM, gradS, A
    );
}




std::tuple<
    arma::mat, double, double, double, double, double, double,
    arma::mat, arma::mat, arma::mat, arma::mat, arma::mat, arma::mat
>
Elbo_grad_LogS(const arma::mat & Y, const arma::mat & X, const arma::mat & R,
              const arma::mat & B, const arma::mat & D, const arma::mat & C, 
              const arma::mat & M, const arma::mat & logS,
              double tolXi
                ) {

    int n = Y.n_rows;
    int p = Y.n_cols;
    int q = M.n_cols;
    arma::mat S = exp(logS) ;
    auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
            Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);


    gradS = S % gradS;

     return std::make_tuple(
        xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective,
        gradB, gradD, gradC, gradM, gradS, A
    );
}

