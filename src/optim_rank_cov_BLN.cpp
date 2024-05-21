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
Rcpp::List nlopt_optimize_rank_BLN(
    const Rcpp::List & data  , // List(Y, R, X)
    const Rcpp::List & params, // List(B, C, M, S)
    const Rcpp::List & config  // List of config values
    ) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
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
    

    // Optimize
    auto objective_and_grad = [&metadata, &X, &Y, &R](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        
        
     
	int n = Y.n_rows;
	int p = Y.n_cols;
	int q = C.n_cols;
	int d = X.n_cols;
	arma::vec XB = X * B;
	arma::mat Mu = arma::mat(XB.memptr(), n, p, false, false);
	arma::vec vecY = arma::vectorise(Y);
	arma::vec vecR = arma::vectorise(R);
	arma::mat A = sqrt(square(Mu + M*C.t()) + S*(C % C).t());
	arma::mat expA = exp(A);
	arma::mat expmA = exp(-A);
	arma::mat D = 1 / (1. + expmA);
	arma::mat E = (1 - expmA)/(4*A%(1 + expmA));
	arma::mat F = Y - D - 2*E%(Mu + M*C.t() - A);
	arma::vec vecF = arma::vectorise(F);
	arma::mat G = log(1. + expmA) + A - D%A;
	
	

	double objective = - (accu(R % ((Y - D) % (Mu + M * C.t()) - E % (S*(C%C).t() + square(Mu + M*C.t() - A)) - G))- 0.5 * accu((M%M + S - log(S)))) + (n*q)/2;
	
	std::cout << objective << std::endl;
	//std::cout << S << std::endl;
	

        metadata.map<B_ID>(grad) = -X.t() * (vecR % vecF) ;
        //metadata.map<B_ID>(grad) = arma::zeros<arma::mat>(1,d);
       	metadata.map<C_ID>(grad) = -((R % F).t() * M - 2 * (R%E).t() * S % C);
        metadata.map<M_ID>(grad) = - (R % F * C - M);
        metadata.map<S_ID>(grad) =  -(-R % E * (C%C)- 1/2 * (1. - 1/S));


        
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
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}
