function [alpha,delta, gamma_a_2, gamma_0_2, gamma_d_2,sigma_1_2, sigma_0_2, mu_0, phi]=extract_mcmc(trace)
% EXTRACT_MCMC is a help function for  that extracts fields of trace as separate
% variables. Needed for parallel loops.
% Input : 
%    trace
% Output : 
%   fields of the trace in separate variables.
        alpha = trace.alpha;
        delta = trace.delta;
        gamma_a_2 = trace.gamma_a_2;
        gamma_0_2 =  trace.gamma_0_2;
        gamma_d_2 =  trace.gamma_d_2;
        sigma_1_2 =  trace.sigma_1_2;
        sigma_0_2 =  trace.sigma_0_2;
        mu_0 =  trace.mu_0;
        phi =  trace.phi;
end