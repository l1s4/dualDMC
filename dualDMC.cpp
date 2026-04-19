#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List simDDMCtrial(double mu_c, int b, int A1, int A2, int tau1, int tau2, 
    double dt, double sigma, int auto1, int auto2) {
        /**
         Simulate single trial of MDMC 
         Note: assuming a = 2
         @param mu_c drift rate of controlled process
         @param b boundary / threshold 
         @param A1, A2 amplitudes of automatic activations
         @param tau1, tau2 timepoint of max automatic activation 
         @param dt time discretization
         @param sigma diffusion constant of superimposed process
         @param auto1, auto2 type of automatic process (1 = congruent, -1 = incongruent)
         @return reaction time rt and decision of trial
         */

        double t = dt; 
        double X = 0.0; 
        int dec = 0;
        std::vector<int> x_traj;    // dynamic vector to store trajectory
        
        // with a1 = a2 = 2
        while (X > -b && X < b) {
            // drifts of automatic processes
            double mu_a1 = auto1 * A1 * exp(-t / tau1) * (exp(1) / tau1) * (1 - t / tau1); 
            double mu_a2 = auto2 * A2 * exp(-t / tau2) * (exp(1) / tau2) * (1 - t / tau2); 
            // superimposed drift
            double mu_t = mu_c + mu_a1 + mu_a2;

            // simulate noisy process
            double dX = mu_t * dt + sigma * sqrt(dt) * Rcpp::rnorm(1)[0];

            X += dX;     // X(t)
            t += dt;     // increase time
            x_traj.push_back(X); 
        }

        if (X < 0){dec = -1;}        // hit -b first (i.e. incorrect response)
        if (X > 0){dec = 1;}         // hit b first (i.e. correct response)

        return Rcpp::List::create(
            Rcpp::Named("rt") = t, 
            Rcpp::Named("dec") = dec, 
            Rcpp::Named("XTraj") = x_traj
        );
    }; 


// [[Rcpp::export]]
Rcpp::List simDDMC(Rcpp::DataFrame df, int N_sim) {

    Rcpp::NumericVector mu_c    = df["mu_c"]; 
    Rcpp::IntegerVector b       = df["b"];
    Rcpp::IntegerVector A1      = df["A1"]; 
    Rcpp::IntegerVector A2      = df["A2"]; 
    Rcpp::IntegerVector tau1    = df["tau1"]; 
    Rcpp::IntegerVector tau2    = df["tau2"]; 
    Rcpp::NumericVector dt      = df["dt"]; 
    Rcpp::NumericVector sigma   = df["sigma"];

    // pre-allocate vectors for parameters
    Rcpp::NumericVector mu_c_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector b_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector A1_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector A2_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector tau1_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector tau2_out(df.nrows()*N_sim); 
    Rcpp::NumericVector dt_out(df.nrows()*N_sim); 
    Rcpp::NumericVector sigma_out(df.nrows()*N_sim); 
    // pre-allocate vectors for simulation results
    Rcpp::IntegerVector cong1_out(df.nrows()*N_sim);
    Rcpp::IntegerVector cong2_out(df.nrows()*N_sim);
    Rcpp::NumericVector rt_out(df.nrows()*N_sim); 
    Rcpp::IntegerVector dec_out(df.nrows()*N_sim); 
    // pre-allocate vector for parameter set id
    Rcpp::IntegerVector set_id(df.nrows()*N_sim); 

    // to sample congruencies
    Rcpp:: IntegerVector in = Rcpp::IntegerVector::create(-1, 1);
    
    int row_out = 0;
    Rcpp::Rcout << "Number of rows: " << df.nrows() << "\n";

    // loop parameter sets
    for (int i = 0; i < df.nrows(); i++) {
        Rcpp::Rcout << "simulating parameter set " << i << "\n";

        // simulate N_sim trials
        for (int j = 0; j < N_sim; j++) {
            int auto1 = Rcpp::sample(in, 1, false)[0];
            int auto2 = Rcpp::sample(in, 1, false)[0];
            Rcpp::List sim = simDDMCtrial(mu_c[i], b[i], A1[i], A2[i], tau1[i], 
                tau2[i], dt[i], sigma[i], auto1, auto2);
            
            mu_c_out[row_out]   = mu_c[i]; 
            b_out[row_out]      = b[i]; 
            A1_out[row_out]     = A1[i]; 
            A2_out[row_out]     = A2[i]; 
            tau1_out[row_out]   = tau1[i]; 
            tau2_out[row_out]   = tau2[i]; 
            dt_out[row_out]     = dt[i]; 
            sigma_out[row_out]  = sigma[i]; 
            
            cong1_out[row_out]  = auto1;
            cong2_out[row_out]  = auto2;
            rt_out[row_out]     = sim["rt"];
            dec_out[row_out]    = sim["dec"]; 
            set_id[row_out]     = i;

            row_out++;
        };
    };

    return Rcpp::DataFrame::create(
        Rcpp::Named("mu_c")     = mu_c_out,
        Rcpp::Named("b")        = b_out,
        Rcpp::Named("A1")       = A1_out,
        Rcpp::Named("A2")       = A2_out,
        Rcpp::Named("tau1")     = tau1_out,
        Rcpp::Named("tau2")     = tau2_out,
        Rcpp::Named("dt")       = dt_out,
        Rcpp::Named("sigma")    = sigma_out, 
        Rcpp::Named("auto1")    = cong1_out,
        Rcpp::Named("auto2")    = cong2_out,
        Rcpp::Named("rt")       = rt_out,
        Rcpp::Named("dec")      = dec_out,
        Rcpp::Named("set_id")   = set_id
    ); 
}; 
