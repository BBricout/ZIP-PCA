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

arma::mat ifelse_mat(const arma::mat & Y, arma::mat & A, arma::mat & nu){
    int n = Y.n_rows;
    int p = Y.n_cols;
    
    arma::Mat<double> E(n, p, arma::fill::none); 
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            if (Y(i,j) == 0) {
                E(i,j) = nu(i,j) - A(i,j);
            } else {
                E(i,j) = 1;
            }
        }
    }
   
   return E;
}

double entropie_logis(arma::mat & xi){
   
    int n = xi.n_rows;
    int p = xi.n_cols;
    
    double H = 0;
    
    for (size_t i = 0; i < n; ++i) {
    	for (size_t j = 0; j < p; ++j){
    		if (xi(i,j) == 0 || xi(i,j) == 1){
    		H = H ;
    		} else {
    		H = H + xi(i,j) * log(xi(i,j)/(1 - xi(i,j))) + log(1 - xi(i,j));
    		}
    	}
    }
     return H ;
}

	




// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_ZIP(
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
    
    std::vector<double> objective_values;
    

    // Optimize
    auto objective_and_grad = [&metadata, &X, &Y, &R, &objective_values](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat D = metadata.map<D_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        

        //std::cout << M(0,0) << std::endl;
        
	
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
        arma::mat pi = 1/(1 + exp(-nu));
        arma::vec vecpi = vectorise(pi);
    	arma::mat E = ifelse_mat(Y, A, nu);
        arma::mat xi = 1/(1 + exp(-E));
        arma::vec vecxi = vectorise(xi);

        
        
        double objective = -(accu(xi % nu - log((1 + exp(-nu))/exp(-nu))) + accu(R % xi % (Y % (mu + M*C.t()) - A - log_fact_Y)) - 1/2 * accu(M%M + S - log(S)) + entropie_logis(xi) + n*q/2);
        std::cout << objective << std::endl;
        objective_values.push_back(objective);
        

        metadata.map<B_ID>(grad) = - X.t() *  (vecR % vecxi % (vecY - vecA)) ;
        metadata.map<D_ID>(grad) = - X.t() *  (vecR % (vecxi - vecpi)) ;
	metadata.map<C_ID>(grad) = -((R % xi % (Y -A)).t() * M - (R % xi % A).t() * S % C);
        metadata.map<M_ID>(grad) = - (R % xi % (Y - A) * C - M);
        metadata.map<S_ID>(grad) =  - 1/2 * (1/S - 1. - R % xi % A * (C%C));

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat D = metadata.copy<D_ID>(parameters.data());
    arma::mat C = metadata.copy<C_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
  

    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("D", D),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("objective_values", objective_values),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}
