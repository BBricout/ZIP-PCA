#include "RcppArmadillo.h"


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
	




// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_rank_cov(
    const Rcpp::List & data  , // List(Y, R, X, O, w)
    const Rcpp::List & params, // List(B, C, M, S)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,p)
    const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (1,d)
    const auto init_C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
    const auto init_S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
    

    

    const auto metadata = tuple_metadata(init_B, init_C, init_M, init_S);
    enum { B_ID, C_ID, M_ID, S_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;
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
            set_from_r_sexp(metadata.map<C_ID>(packed.data()), per_param_list["C"]);
            set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
            set_from_r_sexp(metadata.map<S_ID>(packed.data()), per_param_list["S"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }
    
    std::vector<double> objective_values;
    

    // Optimize
    auto objective_and_grad = [&metadata, &O, &X, &Y, &w, &R, &objective_values](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        
        
       
	
	int numRowsY = Y.n_rows;
	int numColsY = Y.n_cols;
	arma::mat vecw = arma::repmat(w,1,numColsY).as_col();
	arma::vec XB = X * B;
	arma::mat matXB = arma::mat(XB.memptr(), numRowsY, numColsY, false, false);
	arma::vec vecY = arma::vectorise(Y);
	arma::vec vecR = arma::vectorise(R);
        arma::mat S2 = S % S;
        arma::mat Z = O + matXB + M * C.t();
        arma::vec vecZ = arma::vectorise(Z);
        arma::mat A = exp(Z + 0.5 * S2 * (C % C).t());
        arma::vec vecA = exp(vecZ);
        double objective = accu(diagmat(w) * (R % (A - Y % Z))) + 0.5 * accu(diagmat(w) * (M % M + S2 - log(S2) - 1.));
        
        objective_values.push_back(objective);        
    
        

        metadata.map<B_ID>(grad) = (X.each_col() % vecw).t() *  (vecR % (vecA - vecY)) ;
	metadata.map<C_ID>(grad) = (diagmat(w) * R % (A - Y)).t() * M + ((R % A).t() * (S2.each_col() % w)) % C;
        metadata.map<M_ID>(grad) = diagmat(w) * (R % (A - Y) * C + M);
        metadata.map<S_ID>(grad) = diagmat(w) * (S - 1. / S + R % A * (C % C) % S);
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat C = metadata.copy<C_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    arma::mat S2 = S % S;
    arma::mat Sigma = C * (M.t() * (M.each_col() % w) + diagmat(sum(S2.each_col() % w, 0))) * C.t() / accu(w);
    arma::mat Omega = C * inv_sympd((M.t() * (M.each_col() % w) + diagmat(sum(S2.each_col() % w, 0)))/accu(w))  * C.t() ;

    // Element-wise log-likelihood
    int numRowsY = Y.n_rows;
    int numColsY = Y.n_cols;
    arma::vec XB = X * B;
    arma::mat matXB = arma::mat(XB.memptr(), numRowsY, numColsY, false, false);
    arma::mat Z = O + matXB + M * C.t();
    arma::mat A = exp(Z + 0.5 * S2 * (C % C).t());
    arma::mat loglik = arma::sum(R % (Y % Z - A), 1) - 0.5 * sum(M % M + S2 - log(S2) - 1., 1) + ki(Y);

    Rcpp::NumericVector Ji = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(loglik));
    Ji.attr("weights") = w;
    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("Z", Z),
        Rcpp::Named("A", A),
        Rcpp::Named("Sigma", Sigma),
        Rcpp::Named("Omega", Omega),
        Rcpp::Named("Ji", Ji),
        Rcpp::Named("objective_values", objective_values),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}
