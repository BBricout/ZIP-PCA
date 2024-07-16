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
#include "Elbo_gradBB.h"


//--------------------------------------------------------------------------------------------------------------------


// Sortie R de l'Elbo et des gradients


// [[Rcpp::export]]
Rcpp::List Elbo_grad_Rcpp(const Rcpp::List & data, // List(Y, R, X)
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
    const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)




    auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
            Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);


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
        Rcpp::Named("A", A)
    );
}

// [[Rcpp::export]]
Rcpp::List Elbo(const Rcpp::List & data, // List(Y, R, X)
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
    const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)




    auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
            Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);


    return Rcpp::List::create(Rcpp::Named("objective", objective));
}


// [[Rcpp::export]]
Rcpp::List Grad(const Rcpp::List & data, // List(Y, R, X)
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
    const arma::mat & S = Rcpp::as<arma::mat>(params["S"]); // (n,q)




    auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
            Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);


    return Rcpp::List::create(
        Rcpp::Named("gradB", gradB),
        Rcpp::Named("gradD", gradD),
        Rcpp::Named("gradC", gradC),
        Rcpp::Named("gradM", gradM),
        Rcpp::Named("gradS", gradS)
    );
}


//--------------------------------------------------------------------------------------------------------------------
// Optimisation



// [[Rcpp::export]]
Rcpp::List nlopt_optimize_ZIP(
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
    
    	  // Définition des bornes inférieures pour tous les paramètres
    if (config.containsElementNamed("lower_bounds")) {
        auto lower_bounds_r = Rcpp::as<std::vector<double>>(config["lower_bounds"]);
        if (lower_bounds_r.size() != metadata.packed_size) {
            Rcpp::stop("La taille du vecteur lower_bounds ne correspond pas à la taille totale des paramètres.");
        }
        std::vector<double> lower_bounds = lower_bounds_r;

        // Application des bornes à l'optimiseur
        nlopt_set_lower_bounds(optimizer.get(), lower_bounds.data());
    } else {
        std::vector<double> lower_bounds(metadata.packed_size, -HUGE_VAL);
        nlopt_set_lower_bounds(optimizer.get(), lower_bounds.data());
    }
    
 
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
    
    std::cout << "ok" << std::endl;
        
    

    // Optimize
    auto objective_and_grad = [&metadata, &X, &Y, &R, &objective_values, tolXi](const double * params, double * grad) -> double {
    
        const arma::mat B = metadata.map<B_ID>(params);
        const arma::mat D = metadata.map<D_ID>(params);
        const arma::mat C = metadata.map<C_ID>(params);
        const arma::mat M = metadata.map<M_ID>(params);
        const arma::mat S = metadata.map<S_ID>(params);
        
        
    auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
    Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);
    
    objective = -objective;

        objective_values.push_back(- objective);
        std::cout << objective << std::endl;

        metadata.map<B_ID>(grad) = - gradB;
        metadata.map<D_ID>(grad) = - gradD;
	metadata.map<C_ID>(grad) = - gradC;
        metadata.map<M_ID>(grad) = - gradM;
        metadata.map<S_ID>(grad) =  - gradS;
        

        return objective;
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    // Model and variational parameters
    arma::mat B = metadata.copy<B_ID>(parameters.data());
    arma::mat D = metadata.copy<D_ID>(parameters.data());
    arma::mat C = metadata.copy<C_ID>(parameters.data());
    arma::mat M = metadata.copy<M_ID>(parameters.data());
    arma::mat S = metadata.copy<S_ID>(parameters.data());
    
        auto [xi, elbo1, elbo2, elbo3, elbo4, elbo5, objective, gradB, gradD, gradC, gradM, gradS, A] = 
    Elbo_grad(Y, X, R, B, D, C, M, S, tolXi);
    
  	    

    return Rcpp::List::create(
        Rcpp::Named("B", B),
        Rcpp::Named("D", D),
        Rcpp::Named("C", C),
        Rcpp::Named("M", M),
        Rcpp::Named("S", S),
        Rcpp::Named("A", A),
        Rcpp::Named("xi", xi),
        Rcpp::Named("objective", objective),
        Rcpp::Named("objective_values", objective_values),
        Rcpp::Named("monitoring", Rcpp::List::create(
            Rcpp::Named("status", static_cast<int>(result.status)),
            Rcpp::Named("backend", "nlopt"),
            Rcpp::Named("iterations", result.nb_iterations)
        ))
    );
}




