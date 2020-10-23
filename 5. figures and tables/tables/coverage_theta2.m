function coverage_mnar_mnar=coverage_theta2(sim_res, delta_beta, bound_prior, random_values, seq)
% COVERAGE_THETA2 calculates the coverage of the credible intervals 
% Input:
%   sim_res : cell with simulation results
%   delta_beta : true value of delta_beta
%   bound_prior : bound of the prior for delta_beta
%   random_values : a vector of random U(0,1) values
%   seq : sequence for posteriors
% Output:
%   coverage_mnar_mnar : a vector coverages of every simulation setting
    theta2_mnar = sim_res{1}.gen_data_1.delta(9) + (1 - sim_res{1}.gen_data_1.phi) * delta_beta;
    for i = 1:length(sim_res)
        theta2_mnar_l(i, :) = quantile(sim_res{i}.results.delta_diff_trace(seq, :) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * bound_prior, 0.025);
        theta2_mnar_u(i, :) = quantile(sim_res{i}.results.delta_diff_trace(seq,:) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * bound_prior, 0.975);   
    end
    coverage_mnar_mnar = mean((theta2_mnar > theta2_mnar_l).* (theta2_mnar < theta2_mnar_u ), 2);
end