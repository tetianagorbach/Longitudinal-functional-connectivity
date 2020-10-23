function data=generate_data_real(trace, cov_path, n_sub, n_pairs, seq_post, set_rng, ind_sel, phi, returnee) 
% GENERATE_DATA_REAL generates the data from the model fit to the
% Betula data. 
% Inputs:
%   trace : trace used for data generation
%   cov_path : path to the file with the matrix of baseline covariates
%       values of size [number of subjects  * 7] since 7 covariates are 
%       considered in this paper 
%   n_sub : number of subjects to generate. All have data at the baseline,
%       but only some have data at the follow-up.
%   n_pairs : number of node pairs for each subject
%   seq_post : sequence of MCMC iterations used to define the parameters 
%       for data generation
%   set_rng : seed for the MATLAB random number generator
%   ind_sel : subjects from the Betula sample to use the data from
%   phi : proportion of returnees
%   returnee : vector of size [1, n_sub] of 0s for dropouts and 1s for
%   returnees.
% Output:
%   data : structure, fields correspond to the values of 
%       parameters used in data generations, and the field gen_data
%       contains the generated data.
    rng(set_rng)
    %% Define values of the parameters for the data generation    
    % indicators of observed       
    ind_obs = repelem((1:n_sub)', returnee+1); % vector of subject id repeated 1 time for dropouts and 2 times for returnees
    n_obs = length(ind_obs);
    t = ones(n_obs, 1); % t =  1 if the first timepoint, t =  2 if the second timepoint
    for i = 2:n_obs
        if ind_obs(i) == ind_obs(i - 1)
            t(i) = 2;
        end
    end    
    
    % read cov for fixed effects
    cov = readtable(cov_path, 'sheet', 'Sheet1', 'ReadVariableNames', false);
    x_init_preproc = ones(length(ind_sel), 1);
    x_init_preproc(:, 2) = zscore(cov.Var1(ind_sel)); % standardize t1 age
    x_init_preproc(:, 3) = cov.Var2(ind_sel); % Gender
    x_init_preproc(:, 4) = zscore(cov.Var3(ind_sel));% FD
    x_init_preproc(:, 5) = zscore(cov.Var4(ind_sel));% EM
    x_init_preproc(:, 6) = zscore(cov.Var5(ind_sel));% Fl
    x_init_preproc(:, 7) = zscore(cov.Var6(ind_sel));% PS
    x_init_preproc(:, 8) = zscore(cov.Var7(ind_sel));% Bl
        
    x_init = x_init_preproc(ind_obs, :); % covariates for all observed subjects
    x_init(t == 2, 9) = 1;    
    x = repelem(x_init, n_pairs, 1);    
    
    % covariates for random effects 
    x_random = x(:, 1); 
    dim_cov_pop = size(x, 2);
    dim_cov_ind = size(x_random, 2);

    % mean component
    alpha = mean(trace.alpha(:, seq_post), 2);
   
    gamma_a_2 = mean(trace.gamma_a_2(1, seq_post));
    a = normrnd(0, sqrt(gamma_a_2(1)), 1, n_sub);
    
    gamma_0_2 = mean(trace.gamma_0_2(1, seq_post), 2);  
    mu_0i = normrnd(zeros(1, n_sub), sqrt(gamma_0_2), 1, n_sub);
       
    % variances
    sigma_0_2 = mean(trace.sigma_0_2(ind_sel, seq_post), 2);
    sigma_1_2 = mean(trace.sigma_1_2(ind_sel, seq_post), 2);

    % probabilities to be connected
    delta = mean(trace.delta(:, seq_post), 2);
    gamma_d_2 = mean(trace.gamma_d_2(1, seq_post));
    d = normrnd(0, sqrt(gamma_d_2), 1, n_sub);
   
    a_long = a(ind_obs)';
    d_long = d(ind_obs)';
    
    mu_1i = x_init*alpha+sum(x_init(:, 1).*a_long, 2);
    beta_i = x_init*delta+sum(x_init(:, 1).*d_long, 2);

    z_probit = normrnd(repelem(beta_i, n_pairs, 1), 1);
    omega = double(z_probit>0);
    
    %% Simulate the mixture
    con = zeros(n_obs, n_pairs);
    noncon = zeros(n_obs, n_pairs);
    sigma_1_2_ind_obs = sigma_1_2(ind_obs);
    mu_0i_ind_obs = mu_0i(ind_obs);
    sigma_0_2_ind_obs = sigma_0_2(ind_obs);
    for i = 1:n_obs
        con(i, :) = exp(normrnd(mu_1i(i), sqrt(sigma_1_2_ind_obs(i)), 1, n_pairs));
        noncon(i, :) = normrnd(mu_0i_ind_obs(i), sqrt(sigma_0_2_ind_obs(i)), 1, n_pairs);
    end
    
    gen_data = reshape(omega, n_pairs, n_obs)'.*con + (1 - reshape(omega, n_pairs, n_obs)').*noncon;
    
    %% Write to data
    data.alpha = alpha;
    data.delta = delta;
    data.a = a;
    data.d = d;
    data.mu_0 = mu_0i;
    data.gamma_0_2 = gamma_0_2;
    data.gamma_a_2 = gamma_a_2;
    data.gamma_d_2 = gamma_d_2;
    data.mu_i_c = mu_1i;
    data.sigma_0_2 = sigma_0_2;
    data.sigma_1_2 = sigma_1_2;
    data.probit = beta_i;
    data.gen_data = gen_data;
    data.ind_sel = ind_sel;
    data.ind_obs = ind_obs;
    data.x_init = x_init;
    data.phi = phi;
    data.n_sub = n_sub;
    data.n_pairs = n_pairs;
end