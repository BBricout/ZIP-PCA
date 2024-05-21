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

mat log_factorial_matrix(const mat& Y) {
    mat result = Y;
    result.transform([](double val) { return log_factorial(val); });
    return result;
}
	




// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_rank_cov(
    const Rcpp::List & data  , // List(Y, R, X)
    const Rcpp::List & params, // List(B, C, M, S)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
    const auto init_D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
    const auto init_C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
    

    

    const auto metadata = tuple_metadata(init_B, init_D, init_C, init_M, init_S);
    enum { B_ID, D_ID, C_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<D_ID>(parameters.data()) = init_D;
    metadata.map<C_ID>(parameters.data()) = init_C;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<S_ID>(parameters.data()) = init_S;

    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    if(config.containsElementNamed("xtol_abs")) {
        SEXP value = config["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
            set_from_r_sexp(metadata.map<D_ID>(packed.data()), per_param_list["D"]);
            set_from_r_sexp(metadata.map<C_ID>(packed.data()), per_param_list["C"]);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }
    

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &R](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat D = metadata.map<D_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        
        
       
	
	int n = Y.n_rows;
	int p = Y.n_cols;
	int q = M.n_cols;
	arma::vec XB = X * B;
	arma::vec XD = X * D;
	arma::mat mu = arma::mat(XB.memptr(), n, p, false, false);
	arma::mat nu = arma::mat(XD.memptr(), n, p, false, false);
	arma::vec vecY = arma::vectorise(Y);
	arma::vec vecR = arma::vectorise(R);
        arma::mat Z =  mu + M * C.t();
        arma::vec vecZ = arma::vectorise(Z);
        arma::mat A = exp(Z + 0.5 * S * (C % C).t());
        arma::vec vecA = vectorise(A);
        arma::mat log_fact_Y = log_factorial_matrix(Y);
        arma::mat pi = exp(nu)/(1 + exp(nu));
        arma::vec vecpi = vectorise(pi);
        arma::mat E = nu - A + Y % (mu + M*C.t()) - log_fact_Y;
        arma::mat xi = log(E) - log(1. - E);
        arma::vec vecxi = vectorise(xi);
	
        
        
        double objective = -(accu(R % (xi % (nu - log(1 + exp(nu)) + Y % (nu + M*C.t()) - A - log_fact_Y + log(xi) - log (1. - xi)) + log(1 - xi)) - 0.5 * accu(diagmat(w) * (M % M + S - log(S))) + n*q/2);
     
     	std::cout << objective << std::endl;      
    
        

        metadata.map<B_ID>(grad) = - X.t() *  (vecR % vecxi % (vecY - vecA)) ;
        metadata.map<D_ID>(grad) = - X.t() *  (vecR % (vecxi - vecpi)) ;
	metadata.map<C_ID>(grad) = -((R % xi % (Y -A)).t() * M - (R % xi % A).t() * S * C);
        metadata.map<M_ID>(grad) = - (R % xi % (Y - A) * C - M);
        metadata.map<S_ID>(grad) =  - 1/2 * (1/S - 1. - R % xi % A * (C%C));
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat C = metadata.copy<C_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());



    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}
