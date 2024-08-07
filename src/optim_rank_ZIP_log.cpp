#include "RcppArmadillo.h"
#include <cmath>
#include <iostream>
#include <tuple>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

//-------------------------------------------------------------------------------------------------------------------

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
            if (nu(i,j) <= 0.) {
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
    
    double H = 0;
    
    for (size_t i = 0; i < n; ++i) {
    	for (size_t j = 0; j < p; ++j){
    		if (xi(i,j) == 0. || xi(i,j) == 1.){
    		H = H ;
    		} else {
    		H = H + xi(i,j) * log(xi(i,j)/(1. - xi(i,j))) + log(1. - xi(i,j));
    		}
    	}
    }
     return H ;
}

// [[Rcpp::export]]
Rcpp::List ElboB_neg(const Rcpp::List & data, // List(Y, R, X)
                 const Rcpp::List & params, // List(B, C, M, S)
                 double tolXi
                ) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
    const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
    const arma::mat & C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
    const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
    const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)
    
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
    

    
    double objective = - (accu(xi % nu - ifelse_exp(nu)) + 
                        accu(R % xi % (Y % (mu + M*C.t()) - A - log_fact_Y)) - 
                        0.5 * accu(M % M + S - log(S)) + 
                        entropie_logis(xi) + 0.5 * n * q);
                        
    
    arma::mat gradB = - X.t() * (vecR % vecxi % (vecY - vecA));
    arma::mat gradD = - X.t() * (vecR % (vecxi - vecpi));
    arma::mat gradC = - (R % xi % (Y - A)).t() * M - (R % xi % A).t() * S % C;
    arma::mat gradM = - (R % xi % (Y - A) * C - M);
    arma::mat gradS = - 0.5 * (1. / S - 1. - R % xi % A * (C % C));
    
    double elbo1 = accu(xi % nu - ifelse_exp(nu));
    double elbo2 = - 0.5 * accu(M % M + S);
    double elbo3 = accu(R % xi % (Y % (mu + M*C.t()) - A - log_fact_Y));
    double elbo4 = entropie_logis(xi);
    double elbo5 = 0.5 * accu(0.5*log(S%S)) + n * q * 0.5;
    
    return Rcpp::List::create(
    	Rcpp::Named("xi", xi),
        Rcpp::Named("elbo1", elbo1),
        Rcpp::Named("elbo2", elbo2),
        Rcpp::Named("elbo3", elbo3),
        Rcpp::Named("elbo4", elbo4),
        Rcpp::Named("elbo5", elbo5),
        Rcpp::Named("objective", objective),
        Rcpp::Named("gradB", gradB),
        Rcpp::Named("gradD", gradD),
        Rcpp::Named("gradC", gradC),
        Rcpp::Named("gradM", gradM),
        Rcpp::Named("gradS", gradS),
        Rcpp::Named("vecY", vecY),
        Rcpp::Named("vecxi", vecxi), 
        Rcpp::Named("A", A), 
       Rcpp::Named("nu", nu)
    );
}

std::tuple<
    arma::mat, double, double, double, double, double, double,
    arma::mat, arma::mat, arma::mat, arma::mat, arma::mat,
    arma::vec, arma::vec, arma::mat, arma::mat
>
ElboB_neg_cpp(const arma::mat & Y, const arma::mat & X, const arma::mat & R,
              const arma::mat & B, const arma::mat & D, const arma::mat & C, const arma::mat & M, const arma::mat & S,
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
    
    double objective = - (accu(xi % nu - ifelse_exp(nu)) + 
                        accu(R % xi % (Y % (mu + M * C.t()) - A - log_fact_Y)) - 
                        0.5 * accu(M % M + S - log(S)) + 
                        entropie_logis(xi) + 0.5 * n * q);
                        
    arma::mat gradB = - X.t() * (vecR % vecxi % (vecY - vecA));
    arma::mat gradD = - X.t() * (vecR % (vecxi - vecpi));
    arma::mat gradC = - (R % xi % (Y - A)).t() * M - (R % xi % A).t() * S % C;
    arma::mat gradM = - (R % xi % (Y - A) * C - M);
    arma::mat gradS = - 0.5 * (1. / S - 1. - R % xi % A * (C % C));
    
    double elbo1 = accu(xi % nu - ifelse_exp(nu));
    double elbo2 = - 0.5 * accu(M % M + S);
    double elbo3 = accu(R % xi % (Y % (mu + M * C.t()) - A - log_fact_Y));
    double elbo4 = entropie_logis(xi);
    double elbo5 = 0.5 * accu(0.5 * log(S % S)) + n * q * 0.5;
    
    return std::make_tuple(
        xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective,
        gradB, gradD, gradC, gradM, gradS,
        vecY, vecxi, A, nu
    );
}



// [[Rcpp::export]]
Rcpp::List ElboBLogS(const Rcpp::List & data, // List(Y, R, X)
                 const Rcpp::List & params, // List(B, C, M, logS)
                 double tolXi
                ) {
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
    const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
    const arma::mat & C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
    const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
    const arma::mat & logS = Rcpp::as<arma::mat>(params["logS"]); // (n,q)

    int n = Y.n_rows;
    int p = Y.n_cols;
    int q = M.n_cols;
    arma::mat S = exp(logS) ;
    Rcpp::List dataB = data;
    Rcpp::List paramsB = Rcpp::List::create(
    Rcpp::Named("B",B),
    Rcpp::Named("D",D),
    Rcpp::Named("C",C),
    Rcpp::Named("M",M),
    Rcpp::Named("S",S)
    );


    Rcpp::List Elbo_grad = ElboB_neg(dataB, paramsB, tolXi);


    const arma::mat & gradS = Rcpp::as<arma::mat>(Elbo_grad["gradS"]);
    arma::mat gradSfin = S % gradS;

    return Rcpp::List::create(
    Rcpp::Named("objective", Elbo_grad["objective"]),
    Rcpp::Named("gradB",Elbo_grad["gradB"]),
    Rcpp::Named("gradD",Elbo_grad["gradD"]),
    Rcpp::Named("gradC",Elbo_grad["gradC"]),
    Rcpp::Named("gradM",Elbo_grad["gradM"]),
    Rcpp::Named("gradS",gradSfin)
    );
}


// ---------------------------------------------------------------------------------------
// Rank-constrained covariance

// Rank (q) is already determined by param dimensions ; not passed anywhere

// [[Rcpp::export]]
Rcpp::List nlopt_optimize_ZIP_logS(
    const Rcpp::List & data  , // List(Y, R, X)
    const Rcpp::List & params, // List(B, C, M, S)
    const Rcpp::List & config,  // List of config values
    double tolXi
) {
    // Conversion from R, prepare optimization
    const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
    const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
    const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
    const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
    const auto init_D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
    const auto init_C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
    const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
    const auto init_logS = Rcpp::as<arma::mat>(params["logS"]); // (n,q)




    const auto metadata = tuple_metadata(init_B, init_D, init_C, init_M, init_logS);
    enum { B_ID, D_ID, C_ID, M_ID, logS_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<B_ID>(parameters.data()) = init_B;
    metadata.map<D_ID>(parameters.data()) = init_D;
    metadata.map<C_ID>(parameters.data()) = init_C;
    metadata.map<M_ID>(parameters.data()) = init_M;
    metadata.map<logS_ID>(parameters.data()) = init_logS;

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
            set_from_r_sexp(metadata.map<logS_ID>(packed.data()), per_param_list["logS"]);
            set_per_value_xtol_abs(optimizer.get(), packed);
        }
    }

    std::vector<double> objective_values;



 // Optimize
    auto objective_and_grad = [&metadata, &X, &Y, &R, &objective_values, tolXi](const double * params, double * grad) -> double {
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat D = metadata.map<D_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat logS = metadata.map<logS_ID>(params);
        
        arma::mat S = exp(logS);
        
        auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, vecY, vecxi, A, nu] = 
            ElboB_neg_cpp(Y, X, R, B, D, C, M, S, tolXi);
    
        objective_values.push_back(objective);
        std::cout << objective << std::endl;
        
        
        metadata.map<B_ID>(grad) = gradB;
        metadata.map<D_ID>(grad) = gradD;
        metadata.map<C_ID>(grad) = gradC;
        metadata.map<M_ID>(grad) = gradM;
        metadata.map<logS_ID>(grad) = S % gradS;
        

        return objective;
    };
    
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat D = metadata.copy<D_ID>(parameters.data());
    arma::mat C = metadata.copy<C_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<logS_ID>(parameters.data());

    	    int n = Y.n_rows;
	    int p = Y.n_cols;
	    int q = M.n_cols;
	    arma::vec XB = X * B;
	    arma::vec XD = X * D;
	    arma::mat S2 = S%S;
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
	    arma::mat xi = ifelse_mat(Y, A, nu, R, 1e-04);
	    arma::vec vecxi = vectorise(xi);


    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("D", D),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("A", A),
        Rcpp::Named("xi", xi),
        Rcpp::Named("objective_values", objective_values),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}























// // [[Rcpp::export]]
// Rcpp::List ElboBLogS(const Rcpp::List & data, // List(Y, R, X)
//                  const Rcpp::List & params, // List(B, C, M, logS)
//                  double tolXi
//                 ) {
//     const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//     const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
//     const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
//     const arma::mat & B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
//     const arma::mat & D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
//     const arma::mat & C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
//     const arma::mat & M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//     const arma::mat & logS = Rcpp::as<arma::mat>(params["logS"]); // (n,q)
// 
//     arma::mat S = exp(logS);
// 
//     int n = Y.n_rows;
//     int p = Y.n_cols;
//     int q = M.n_cols;
//     arma::vec XB = X * B;
//     arma::vec XD = X * D;
//     arma::mat mu = arma::mat(XB.memptr(), n, p, false, false);
//     arma::mat nu = arma::mat(XD.memptr(), n, p, false, false);
//     arma::vec vecY = arma::vectorise(Y);
//     arma::vec vecR = arma::vectorise(R);
//     arma::mat Z = mu + M * C.t();
//     arma::vec vecZ = arma::vectorise(Z);
//     arma::mat A = exp(Z + 0.5 * S * (C % C).t());
//     arma::vec vecA = vectorise(A);
//     arma::mat log_fact_Y = log_factorial_matrix(Y);
//     arma::mat pi = 1./(1. + exp(-nu));
//     arma::vec vecpi = vectorise(pi);
//     arma::mat xi = ifelse_mat(Y, A, nu, R, tolXi);
//     arma::vec vecxi = vectorise(xi);
// 
// 
// 
//     double objective = - (accu(xi % nu - ifelse_exp(nu)) +
//                           accu(R % xi % (Y % (mu + M*C.t()) - A - log_fact_Y)) -
//                           0.5 * accu(M % M + S - log(S)) +
//                           entropie_logis(xi) + 0.5 * n * q);
// 
// 
//     arma::mat gradB = - X.t() * (vecR % vecxi % (vecY - vecA));
//     arma::mat gradD = - X.t() * (vecR % (vecxi - vecpi));
//     arma::mat gradC = - (R % xi % (Y - A)).t() * M - (R % xi % A).t() * S % C;
//     arma::mat gradM = - (R % xi % (Y - A) * C - M);
//     arma::mat gradS = - 0.5 * (1. / S - 1. - R % xi % A * (C % C));
//     gradS = S%gradS ;
// 
// 
//     return Rcpp::List::create(
//     Rcpp::Named("objective", objective),
//     Rcpp::Named("gradB",gradB),
//     Rcpp::Named("gradD",gradD),
//     Rcpp::Named("gradC",gradC),
//     Rcpp::Named("gradM",gradM),
//     Rcpp::Named("gradS",gradS)
//     );
// }






// // ---------------------------------------------------------------------------------------
// // Rank-constrained covariance
// 
// // Rank (q) is already determined by param dimensions ; not passed anywhere
// 
// // [[Rcpp::export]]
// Rcpp::List nlopt_optimize_ZIP_logS(
//     const Rcpp::List & data  , // List(Y, R, X)
//     const Rcpp::List & params, // List(B, C, M, S)
//     const Rcpp::List & config,  // List of config values
//     double tolXi
// ) {
//     // Conversion from R, prepare optimization
//     const arma::mat & Y = Rcpp::as<arma::mat>(data["Y"]); // responses (n,p)
//     const arma::mat & R = Rcpp::as<arma::mat>(data["R"]); // missing data (n,p)
//     const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (np,d)
//     const auto init_B = Rcpp::as<arma::mat>(params["B"]); // (1,d) régresseurs pour la Poisson
//     const auto init_D = Rcpp::as<arma::mat>(params["D"]); // (1,d) régresseurs pour la logistique
//     const auto init_C = Rcpp::as<arma::mat>(params["C"]); // (p,q)
//     const auto init_M = Rcpp::as<arma::mat>(params["M"]); // (n,q)
//     const auto init_logS = Rcpp::as<arma::mat>(params["logS"]); // (n,q)
// 
// 
// 
// 
//     const auto metadata = tuple_metadata(init_B, init_D, init_C, init_M, init_logS);
//     enum { B_ID, D_ID, C_ID, M_ID, logS_ID }; // Names for metadata indexes
// 
//     auto parameters = std::vector<double>(metadata.packed_size);
//     metadata.map<B_ID>(parameters.data()) = init_B;
//     metadata.map<D_ID>(parameters.data()) = init_D;
//     metadata.map<C_ID>(parameters.data()) = init_C;
//     metadata.map<M_ID>(parameters.data()) = init_M;
//     metadata.map<logS_ID>(parameters.data()) = init_logS;
// 
//     auto optimizer = new_nlopt_optimizer(config, parameters.size());
// 
// 
//     if(config.containsElementNamed("xtol_abs")) {
//         SEXP value = config["xtol_abs"];
//         if(Rcpp::is<double>(value)) {
//             set_uniform_xtol_abs(optimizer.get(), Rcpp::as<double>(value));
//         } else {
//             auto per_param_list = Rcpp::as<Rcpp::List>(value);
//             auto packed = std::vector<double>(metadata.packed_size);
//             set_from_r_sexp(metadata.map<B_ID>(packed.data()), per_param_list["B"]);
//             set_from_r_sexp(metadata.map<D_ID>(packed.data()), per_param_list["D"]);
//             set_from_r_sexp(metadata.map<C_ID>(packed.data()), per_param_list["C"]);
//             set_from_r_sexp(metadata.map<M_ID>(packed.data()), per_param_list["M"]);
//             set_from_r_sexp(metadata.map<logS_ID>(packed.data()), per_param_list["logS"]);
//             set_per_value_xtol_abs(optimizer.get(), packed);
//         }
//     }
// 
//     std::vector<double> objective_values;
// 
// 
// 
//     // Optimize
//     auto objective_and_grad = [&metadata, &X, &Y, &R, &objective_values, tolXi](const double * params, double * grad) -> double {
//         const arma::mat B = metadata.map<B_ID>(params);
//         const arma::mat D = metadata.map<D_ID>(params);
//         const arma::mat C = metadata.map<C_ID>(params);
//         const arma::mat M = metadata.map<M_ID>(params);
//         const arma::mat logS = metadata.map<logS_ID>(params);
// 
//     arma::mat S = exp(logS);
// 
//     int n = Y.n_rows;
//     int p = Y.n_cols;
//     int q = M.n_cols;
//     arma::vec XB = X * B;
//     arma::vec XD = X * D;
//     arma::mat mu = arma::mat(XB.memptr(), n, p, false, false);
//     arma::mat nu = arma::mat(XD.memptr(), n, p, false, false);
//     arma::vec vecY = arma::vectorise(Y);
//     arma::vec vecR = arma::vectorise(R);
//     arma::mat Z = mu + M * C.t();
//     arma::vec vecZ = arma::vectorise(Z);
//     arma::mat A = exp(Z + 0.5 * S * (C % C).t());
//     arma::vec vecA = vectorise(A);
//     arma::mat log_fact_Y = log_factorial_matrix(Y);
//     arma::mat pi = 1./(1. + exp(-nu));
//     arma::vec vecpi = vectorise(pi);
//     arma::mat xi = ifelse_mat(Y, A, nu, R, tolXi);
//     arma::vec vecxi = vectorise(xi);
// 
// 
// 
//     double objective = - (accu(xi % nu - ifelse_exp(nu)) +
//                         accu(R % xi % (Y % (mu + M*C.t()) - A - log_fact_Y)) -
//                         0.5 * accu(M % M + S - log(S)) +
//                         entropie_logis(xi) + 0.5 * n * q);
// 
// 
//     arma::mat gradB = - X.t() * (vecR % vecxi % (vecY - vecA));
//     arma::mat gradD = - X.t() * (vecR % (vecxi - vecpi));
//     arma::mat gradC = - (R % xi % (Y - A)).t() * M - (R % xi % A).t() * S % C;
//     arma::mat gradM = - (R % xi % (Y - A) * C - M);
//     arma::mat gradS = - 0.5 * (1. / S - 1. - R % xi % A * (C % C));
// 
// 
//         objective_values.push_back(objective);
//         std::cout<< objective << std::endl;
// 
// 
//         metadata.map<B_ID>(grad) =  gradB;
//         metadata.map<D_ID>(grad) =  gradD;
// 	metadata.map<C_ID>(grad) = gradC;
//         metadata.map<M_ID>(grad) = gradM;
//         metadata.map<logS_ID>(grad) =  S % gradS;
// 
// 
//         return objective;
//     };
//     OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);
// 
//     // Model and variational parameters
//     arma::mat B = metadata.copy<B_ID>(parameters.data());
//     arma::mat D = metadata.copy<D_ID>(parameters.data());
//     arma::mat C = metadata.copy<C_ID>(parameters.data());
//     arma::mat M = metadata.copy<M_ID>(parameters.data());
//     arma::mat S = metadata.copy<logS_ID>(parameters.data());
// 
//     	    int n = Y.n_rows;
// 	    int p = Y.n_cols;
// 	    int q = M.n_cols;
// 	    arma::vec XB = X * B;
// 	    arma::vec XD = X * D;
// 	    arma::mat S2 = S%S;
// 	    arma::mat mu = arma::mat(XB.memptr(), n, p, false, false);
// 	    arma::mat nu = arma::mat(XD.memptr(), n, p, false, false);
// 	    arma::vec vecY = arma::vectorise(Y);
// 	    arma::vec vecR = arma::vectorise(R);
// 	    arma::mat Z = mu + M * C.t();
// 	    arma::vec vecZ = arma::vectorise(Z);
// 	    arma::mat A = exp(Z + 0.5 * S * (C % C).t());
// 	    arma::vec vecA = vectorise(A);
// 	    arma::mat log_fact_Y = log_factorial_matrix(Y);
// 	    arma::mat pi = 1./(1. + exp(-nu));
// 	    arma::vec vecpi = vectorise(pi);
// 	    arma::mat xi = ifelse_mat(Y, A, nu, R, 1e-04);
// 	    arma::vec vecxi = vectorise(xi);
// 
// 
//     return Rcpp::List::create(
//         Rcpp::Named("B", B),
//         Rcpp::Named("D", D),
//         Rcpp::Named("C", C),
//         Rcpp::Named("M", M),
//         Rcpp::Named("S", S),
//         Rcpp::Named("A", A),
//         Rcpp::Named("xi", xi),
//         Rcpp::Named("objective_values", objective_values),
//         Rcpp::Named("monitoring", Rcpp::List::create(
//             Rcpp::Named("status", static_cast<int>(result.status)),
//             Rcpp::Named("backend", "nlopt"),
//             Rcpp::Named("iterations", result.nb_iterations)
//         ))
//     );
// }






