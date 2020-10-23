function coverage_mnar_mnar=coverage_theta1(sim_res, delta_mu, lower_bound_prior, random_values, seq)
% COVERAGE_THETA1 calculates the coverage of the credible intervals
% Input:
%   sim_res : cell with simulation results
%   delta_mu : true value of delta_mu
%   lower_bound_prior : lower bound of the prior for delta_mu
%   random_values : a vector of random U(0,1) values
%   seq : sequence for posteriors
% Output:
%   coverage_mnar_mnar : a vector coverages of every simulation setting
    theta1_mnar = sim_res{1}.gen_data_1.alpha(9) + (1 - sim_res{1}.gen_data_1.phi) * delta_mu;
    for i = 1:length(sim_res)
        theta1_mnar_l(i, :) = quantile(sim_res{i}.results.alpha_diff_trace(seq,:) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * lower_bound_prior, 0.025);
        theta1_mnar_u(i, :) = quantile(sim_res{i}.results.alpha_diff_trace(seq,:) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * lower_bound_prior, 0.975);
    end
    coverage_mnar_mnar = mean((theta1_mnar > theta1_mnar_l).* (theta1_mnar < theta1_mnar_u), 2);
end