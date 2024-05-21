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



// [[Rcpp::export]]
Rcpp::List nlopt_optimize_logit(
    const Rcpp::List & data  , // List(Y, R, X)
    const Rcpp::List & params, // List(B)
    const Rcpp::List & config  // List of config values
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (1,d)

    
    const auto metadata = tuple_metadata(init_B);
    enum { B_ID}; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;


    auto optimizer = new_nlopt_optimizer(config, parameters.size());
    if(config.containsElementNamed("xtol_abs")) {
        SEXP value = config["xtol_abs"];
        if(Rcpp::is<double>(value)) {
            set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
        } else {
            auto per_param_list = Rcpp::as<Rcpp::List>(value);
            auto packed = std::vector<double>(metadata.packed_size);
            set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }
    

    // Optimize
    auto objective_and_grad = [&metadata, &X, &Y, &R](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        
        int n = Y.n_rows;
	int p = Y.n_cols;
        arma::vec XB = X*B;
        arma::mat Mu = arma::mat(XB.memptr(), n, p, false, false);
        arma::mat expMu = exp(Mu);
        arma::mat A = 1+expMu ;
        arma::mat D = exp(Mu)/(1 + exp(Mu));
        arma::vec vecY = arma::vectorise(Y);
	arma::vec vecR = arma::vectorise(R);
	arma::mat F = D-Y;
	arma::vec vecF = arma::vectorise(F);
        
        
        
        double objective = accu(R % (log(A) - Y % Mu));
        
        

        metadata.map<B_ID>(grad) = X.t() * (vecR % vecF);
        std::cout << metadata.map<B_ID>(grad) << std::endl;
        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
    
        // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());


    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}
