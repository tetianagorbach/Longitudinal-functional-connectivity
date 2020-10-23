function results=sim_study_2t(n_mcmc_it, n_sets_gen, trace_for_data_simulation_path, seq_post, cov_path, set_rng, name_results, phi, n_sub_vect, n_pairs_vect, separated)
% SIM_STUDY_2T  simulation study
% Input: 
    % n_mcm_it : number of MCMC iterations for each data replication;
    % n_sets_gen : number of data replications
    % trace_for_data_simulation_path : path to the trace used for data
    % generation
    % seq_post : sequence of MCMC iterations from the trace in
    % trace_for_data_simulation_path to be used in data generation
    % cov_path : path to the file with covariates information
    % set_rng : seed for the MATLAB random number generator
    % name_results : the character sequence to start the name of the output
    % file
    % phi : proportion of returnees
    % n_sub_vect : vector of number of subjects to generate with data at
    % the baseline
    % n_pairs_vect : vector of number of node pairs to generate for each subject
    % separated: 1 if generate separated mixture components (Design 1), otherwise 
    % generate from the model fit to the Betula data (Design 2)
% Output: 
%   For each number of subjects and number of node pairs the results are saved in the current directory in the
%   file named by concatenation of the strings name_results, "seed", set_rng
%   "sub", n_sub, "pairs", n_pairs ".mat".
%   Each file contains a structure with the fields:
%       gen_data_1 : example of the generated data
%       results : returnee: matrix with each row having of 1 for returnees
%                    and 0 otherwise for a data replication
%                 model parameters: matrices of size [9, 4, n_sets_gen] 
%                   with for each generated data replication and parameter 
%                   is a vector of 2.5%, 97.5% quantiles and the mean of 
%                   a posterior distribution, and the true value of the
%                   parameter
%                 alpha_diff_trace : columns are traces for theta_mu
%                 delta_diff_trace : columns are traces for theta_beta
%                 phi_trace : columns are traces for phi
%       set_rng_data_sim : seeds for the random number generator 
%                           used in data generation.   

rng(set_rng)
trace_for_data_gen=importdata(trace_for_data_simulation_path);
for n_sub = n_sub_vect
    ind_sel = datasample(1:310, n_sub, 'Replace', false); % select subjects to use the data from
    returnee = binornd(1, phi, n_sets_gen, n_sub);  % define returnees (returnee == 1)
    for n_pairs = n_pairs_vect
        % initialize
        results =  [];
        results_alpha =  [];
        results_delta =  [];
        results_gamma_a_2 =  [];
        results_gamma_0_2 =  [];
        results_gamma_d_2 =  [];
        results_sigma_1_2 =  [];
        results_sigma_0_2 =  [];
        results_mu_0 =  [];
        results_alpha_diff =  [];
        results_delta_diff =  [];
        results_phi =  [];
        set_rng_data_sim =  randi(10000, n_sets_gen,1); % set seeds for each parfor iteration
        if separated % dataset to save
            gen_data_1 =  generate_data_sep(trace_for_data_gen, cov_path, n_sub, n_pairs, seq_post, set_rng_data_sim(1), ind_sel, phi, returnee(1,:));
        else
            gen_data_1 =  generate_data_real(trace_for_data_gen, cov_path, n_sub, n_pairs, seq_post, set_rng_data_sim(1), ind_sel, phi, returnee(1,:));
        end
        % run parallel loop of simulations
        parfor k =  1:n_sets_gen 
            if separated
                gen_data =  generate_data_sep(trace_for_data_gen, cov_path, n_sub, n_pairs, seq_post, set_rng_data_sim(k), ind_sel, phi, returnee(k,:)); 
            else
                gen_data =  generate_data_real(trace_for_data_gen, cov_path, n_sub, n_pairs, seq_post, set_rng_data_sim(k), ind_sel, phi, returnee(k,:)); 
            end
            % starting value for the connected component parameters
            top_10 =  zeros(n_sub, n_pairs);
            for i =  1:size(gen_data.gen_data,1)
               [~,I] =  sort(gen_data.gen_data(i,:), 'descend');
               top_10(i, I(1:ceil(n_pairs/10))) =  gen_data.gen_data(i,I(1:ceil(n_pairs/10)));
            end
            start_values =  [mean(log(top_10(top_10~=0))), 0, -1.28, 0];  
            % run estimation and extract information
            [alpha, delta, gamma_a_2, gamma_0_2, gamma_d_2, sigma_1_2, sigma_0_2, mu_0, phi_est]...
                = extract_mcmc(mcmc_2t(n_mcmc_it, gen_data, set_rng_data_sim(k), start_values));
            % extract 2.5%, 97.5% quantiles and the mean of 
            % a posterior distribution, and the true value of the
            % parameter
            seq_post_sim =  n_mcmc_it/2:5:n_mcmc_it; % sequence for posterior inferences
            results_alpha(:,:,k) =  [quantile(alpha(:, seq_post_sim)', [0.025,0.975])',mean(alpha(:, seq_post_sim),2), gen_data.alpha];
            results_delta(:,:,k) =  [quantile(delta(:, seq_post_sim)', [0.025,0.975])',mean(delta(:, seq_post_sim),2), gen_data.delta];    
            results_gamma_a_2(:,:,k) =  [quantile(gamma_a_2(:,seq_post_sim)', [0.025,0.975]),mean(gamma_a_2(:,seq_post_sim),2),gen_data.gamma_a_2];
            results_gamma_d_2(:,:,k) =  [quantile(gamma_d_2(:,seq_post_sim)', [0.025,0.975]),mean(gamma_d_2(:,seq_post_sim),2),gen_data.gamma_d_2];
            results_gamma_0_2(:,:,k) =  [quantile(gamma_0_2(:,seq_post_sim)', [0.025,0.975]),mean(gamma_0_2(:,seq_post_sim),2),gen_data.gamma_0_2];
            results_sigma_1_2(:,:,k) =  [quantile(sigma_1_2(:,seq_post_sim)', [0.025,0.975])',mean(sigma_1_2(:, seq_post_sim),2), gen_data.sigma_1_2];
            results_sigma_0_2(:,:,k) =  [quantile(sigma_0_2(:,seq_post_sim)', [0.025,0.975])',mean(sigma_0_2(:, seq_post_sim),2), gen_data.sigma_0_2];
            results_mu_0(:,:,k) =  [quantile(mu_0(:, seq_post_sim)', [0.025,0.975])',mean(mu_0(:, seq_post_sim),2), gen_data.mu_0'];
            results_phi(:,:,k) =  [quantile(phi_est(:, seq_post_sim)', [0.025,0.975]),mean(phi_est(:, seq_post_sim),2), gen_data.phi];
            results_theta_mu(:,k) =  alpha(9,:);
            results_theta_beta(:,k) =  delta(9,:);
            results_phi_trace(:,k) =  phi_est;    
        end
        % gather outputs of all loops
        results.returnee =  returnee;
        results.alpha =  results_alpha;
        results.delta =  results_delta;
        results.gamma_a_2 =  results_gamma_a_2;
        results.gamma_0_2 =  results_gamma_0_2;
        results.gamma_d_2 =  results_gamma_d_2;
        results.sigma_1_2 =  results_sigma_1_2;
        results.sigma_0_2 =  results_sigma_0_2;
        results.mu_0 =  results_mu_0;
        results.phi =  results_phi;
        results.alpha_diff_trace =  results_theta_mu;
        results.delta_diff_trace =  results_theta_beta;
        results.phi_trace =  results_phi_trace;
        % save the output
        save(strcat(char(name_results), sprintf('seed%dsub%dpairs%d',0, n_sub, n_pairs)),'results', 'gen_data_1', 'set_rng_data_sim');
    end
end
end
