function length_cri=length_theta1(sim_res, lower_bound_prior, random_values, seq)
% LENGTH_THETA1 calculates the length of the credible intervals for theta_1
% Input:
%   sim_res : cell with simulation results
%   delta_mu : true value of delta_mu
%   lower_bound_prior : lower bound of the prior for delta_mu
%   random_values : a vector of random U(0,1) values
%   seq : sequence for posteriors
% Output:
%    length_cri : a vector coverages of every simulation setting
    for i = 1:length(sim_res)
        theta1_mnar_l(i, :) = quantile(sim_res{i}.results.alpha_diff_trace(seq, :) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * (lower_bound_prior), 0.025);
        theta1_mnar_u(i, :) = quantile(sim_res{i}.results.alpha_diff_trace(seq, :) +...
            (1 - sim_res{i}.results.phi_trace(seq, :)).*random_values{i} * (lower_bound_prior), 0.975);
    end
    length_cri = mean(theta1_mnar_u - theta1_mnar_l, 2);
end