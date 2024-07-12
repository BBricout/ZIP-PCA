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

//----------------------------------------------------------------------------------------
// From matrix to vector

// [[Rcpp::export]]
Rcpp::NumericVector MatrixToVector(const arma::mat & matrix) {
    int n = matrix.n_rows;
    int p = matrix.n_cols;

    // Check if the input matrix is not empty
    if (n == 0 || p == 0) {
        Rcpp::stop("Input matrix is empty");
    }

    // Reshape the matrix into a column vector
    arma::vec vectorized = arma::vectorise(matrix);

    // Convert the Armadillo vector to an Rcpp numeric vector
    Rcpp::NumericVector result(vectorized.begin(), vectorized.end());

    return result;
}

// From vector to matrix

// [[Rcpp::export]]

// Function to convert an arma::vec to an Rcpp NumericMatrix
Rcpp::NumericMatrix VectorToMatrix(const arma::vec & vector, int n, int p) {
  // Create an arma::mat with the specified dimensions
  arma::mat matrix = arma::reshape(vector, n, p);
  
  // Convert the arma::mat to an Rcpp NumericMatrix
   Rcpp::NumericMatrix result(n, p, matrix.memptr()); 
  
  return result;
}

// Calcul de log(Y!)


double log_factorial(double n) {
    return lgamma(n + 1);
}

arma::mat log_factorial_matrix(const arma::mat& Y) {
    arma::mat result = Y;
    result.transform([](double val) { return log_factorial(val); });
    return result;
}

arma::mat ifelse_mat(const arma::mat& Y, const arma::mat& A, const arma::mat& nu, const arma::mat& R, double tolXi){
    int n = Y.n_rows;
    int p = Y.n_cols;
    
    arma::Mat<double> xi(n, p, arma::fill::none); 
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            if (Y(i,j)*R(i, j) == 0){
                xi(i,j) = nu(i,j) - R(i,j) * A(i,j);
                if (xi(i,j) >= 0){
                	xi(i, j) = 1./(1. + exp(-xi(i, j)));
                } else {
                	xi(i,j) = exp(xi(i,j))/ (exp(xi(i,j)) + 1.);
                }
                xi(i,j) = (xi(i,j) + tolXi) / (1. + 2.*tolXi) ;
            } else {
                xi(i,j) = 1.;
            }
        }
    }
    
    
   return xi;
}

arma::mat ifelse_exp(const arma::mat& nu){
    int n = nu.n_rows;
    int p = nu.n_cols;
    
    arma::Mat<double> F(n, p, arma::fill::none); 
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            if (nu(i,j) <= 0) {
                F(i,j) = log(exp(nu(i,j)) + 1.);
            } else {
                F(i,j) = log((1. + exp(-nu(i,j)))/exp(-nu(i,j)));
            }
        }
    }
   
   return F;
}


double entropie_logis(arma::mat & xi){
   
    int n = xi.n_rows;
    int p = xi.n_cols;
    
    double H = 0.;
    
    for (size_t i = 0; i < n; ++i) {
    	for (size_t j = 0; j < p; ++j){
    		if (xi(i,j) == 0. || xi(i,j) == 1.){
    		H = H ;
    		} else {
    		H = H - (xi(i,j) * log(xi(i,j)) + (1 - xi(i,j))*log(1 - xi(i,j))) ; 
    		//H = H + xi(i,j) * log(xi(i,j)/(1. - xi(i,j))) + log(1. - xi(i,j));
    		}
    	}
    }
     return H ;
}
