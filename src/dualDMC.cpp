#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List simDDMCtrial(double mu_c, int b, int A1, int A2, int tau1, int tau2, 
    double dt, double sigma, int auto1, int auto2) {
        /**
         @brief Simulate a single trial of dualDMC (DDMC). Note: assuming a1 = a2 = 2
         @param mu_c drift rate of controlled process
         @param b boundary / threshold 
         @param A1, A2 amplitudes of automatic activations
         @param tau1, tau2 timepoint of max automatic activations 
         @param dt time discretization
         @param sigma diffusion constant of superimposed process
         @param auto1, auto2 type of automatic process (1 = congruent, -1 = incongruent)
         @return reaction time rt (ms) and decision of trial (1 = correct, -1 = incorrect)
         */

        double t = dt; 
        double X = 0.0; 
        int dec = 0;                   // decision
        double t_max = 7000;            // set max time for safety
        std::vector<double> x_traj;    // to store trajectory
        x_traj.reserve(t_max / dt); 

        const double e = exp(1.0);
        const double sqrt_dt = sqrt(dt); 
        
        // with a1 = a2 = 2
        while (X > -b && X < b && t <= t_max) {
            // drifts of automatic processes
            double mu_a1 = auto1 * A1 * exp(-t / tau1) * (e / tau1) * (1 - t / tau1); 
            double mu_a2 = auto2 * A2 * exp(-t / tau2) * (e / tau2) * (1 - t / tau2); 
            double mu_t = mu_c + mu_a1 + mu_a2;     // superimposed drift

            // simulate noisy process
            double dX = mu_t * dt + sigma * sqrt_dt * R::rnorm(0.0, 1.0);

            X += dX;     // X(t)
            t += dt;     // increase time
            x_traj.push_back(X);
        }

        if (X <= -b){dec = -1;}     // hit -b first (i.e. incorrect response)
        else if (X >= b){dec = 1;}       // hit b first (i.e. correct response)
        else {dec = 0; Rcpp::Rcout << "No decision made!" << "\n";} // no decision

        return Rcpp::List::create(
            Rcpp::Named("rt") = t, 
            Rcpp::Named("dec") = dec,
            Rcpp::Named("XTraj") = x_traj
        );
    }; 


// [[Rcpp::export]]
Rcpp::List simDDMC(Rcpp::DataFrame df, int N_sim) {
    /**
     @brief Simulate multiple DDMC trials. Note: assuming a1 = a2 = 2
     @param df R DataFrame with names of DDMC parameters as columnames
     @param N_sim number of simulations per parameter set
     @return reaction time rt (ms) and decision of trial (1 = correct, -1 = incorrect)
     */

    Rcpp::NumericVector mu_c = df["mu_c"], dt = df["dt"], sigma = df["sigma"];
    Rcpp::IntegerVector b = df["b"], A1 = df["A1"], A2 = df["A2"], 
                        tau1 = df["tau1"], tau2 = df["tau2"]; 

    int n_rows = df.nrows(); 
    int n_out = n_rows * N_sim;  

    // pre-allocate vectors for parameters
    Rcpp::NumericVector mu_c_out(n_out), b_out(n_out), dt_out(n_out), 
                        sigma_out(n_out), rt_out(n_out); 
    Rcpp::IntegerVector A1_out(n_out), A2_out(n_out), tau1_out(n_out), 
                        tau2_out(n_out), cong1_out(n_out), cong2_out(n_out), 
                        dec_out(n_out), set_id(n_out); 

    // sample congruencies
    Rcpp::IntegerVector in = Rcpp::IntegerVector::create(-1, 1);
    Rcpp::IntegerVector auto1_vals = Rcpp::sample(in, n_out, true); 
    Rcpp::IntegerVector auto2_vals = Rcpp::sample(in, n_out, true); 
    
    int out_idx = 0;

    // loop parameter sets
    for (int i = 0; i < n_rows; i++) {
        Rcpp::Rcout << "simulating parameter set " << i << "\n";

        // simulate N_sim trials
        for (int j = 0; j < N_sim; j++) {
            int auto1 = auto1_vals[out_idx];
            int auto2 = auto2_vals[out_idx];

            Rcpp::List sim = simDDMCtrial(mu_c[i], b[i], A1[i], A2[i], tau1[i], 
                tau2[i], dt[i], sigma[i], auto1, auto2);
            
            mu_c_out[out_idx]   = mu_c[i]; 
            b_out[out_idx]      = b[i]; 
            A1_out[out_idx]     = A1[i]; 
            A2_out[out_idx]     = A2[i]; 
            tau1_out[out_idx]   = tau1[i]; 
            tau2_out[out_idx]   = tau2[i]; 
            dt_out[out_idx]     = dt[i]; 
            sigma_out[out_idx]  = sigma[i]; 
            
            cong1_out[out_idx]  = auto1;
            cong2_out[out_idx]  = auto2;
            rt_out[out_idx]     = sim["rt"];
            dec_out[out_idx]    = sim["dec"]; 
            set_id[out_idx]     = i;

            out_idx++;
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


// [[Rcpp::export]]
Rcpp::List simDDMCactivation(int N, double mu_c, int b, int A1, int A2, 
    int tau1, int tau2, int auto1, int auto2) {
        /**
         @brief Simulate activation functions of DDMC. Note: assuming a1 = a2 = 2
         @param N number of timepoints to be simulated [ms]
         @param mu_c drift rate of controlled process
         @param b boundary / threshold 
         @param A1, A2 amplitudes of automatic activations
         @param tau1, tau2 timepoint of max automatic activation 
         @param dt time discretization
         @param auto1, auto2 type of automatic process (1 = congruent, -1 = incongruent)
         @return reaction time rt and decision of trial (-1 = incorrect, 1 = correct)
         */

        double X = 0.0; 
        std::vector<double> auto1_traj;         // first automatic process
        std::vector<double> auto2_traj;         // second automatic process

        std::vector<double> cont_traj(N);       // controlled process
        for (size_t i = 0; i < cont_traj.size(); i++) {
            cont_traj[i] = (i + 1) * mu_c;
        }
        std::vector<double> super_traj;
        
        for (double t = 0; t < N; t++) {
            // activation of automatic processes
            double mu_a1 = auto1 * A1 * exp(-t / tau1) * ((t * exp(1)) / tau1); 
            double mu_a2 = auto2 * A2 * exp(-t / tau2) * ((t * exp(1)) / tau2); 
            // superimposed activation
            double mu_t = mu_c + mu_a1 + mu_a2;

            auto1_traj.push_back(mu_a1); 
            auto2_traj.push_back(mu_a2); 
            super_traj.push_back(mu_t); 
        }

        return Rcpp::List::create(
            Rcpp::Named("cont_traj") = cont_traj, 
            Rcpp::Named("auto1_traj") = auto1_traj,
            Rcpp::Named("auto2_traj") = auto2_traj, 
            Rcpp::Named("super_traj") = super_traj
        );
    }; 