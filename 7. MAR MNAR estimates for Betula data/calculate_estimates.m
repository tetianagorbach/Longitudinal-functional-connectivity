% This file provides estimates from the model fitted to the Betula data.
% Note: theta is the last element of alphas and deltas.
% Input:
%   seq_post : MCMC iterations to use in the calculation of the scale 
%   reduction factor  

%% read the traces
% add the folder with "assess_convergence_and_calc_estimates" and 3 traces to search path
start_names_traces = 'trace_2t_with_missing_seed';
trace0=importdata(strcat(start_names_traces, '0.mat'));
seq_post =  1000:5:6000; % sequence of MCMC iterations for posterior inference

%% Parameters' estimates
mean(trace0.alpha(1:8, seq_post),2) % estimates of alpha 
quantile(trace0.alpha(1:8, seq_post)', 0.025)' % estimates of alpha 
quantile(trace0.alpha(1:8, seq_post)', 0.975)' % estimates of alpha 

mean(trace0.delta(1:8, seq_post),2) % estimates of delta 
quantile(trace0.delta(1:8, seq_post)', 0.025)' % estimates of delta 
quantile(trace0.delta(1:8, seq_post)', 0.975)' % estimates of delta

%% Estimates for thetas, MAR
mean(trace0.alpha(9, seq_post),2) 
% Credible interval, theta_1
quantile(trace0.alpha(9, seq_post)', 0.025)  
quantile(trace0.alpha(9, seq_post)', 0.975)  
% width of the credible interval
quantile(trace0.alpha(9, seq_post)', 0.975)  - quantile(trace0.alpha(9, seq_post)', 0.025)  

mean(trace0.delta(9, seq_post),2) 
% Credible interval, theta_2
quantile(trace0.delta(9, seq_post)', 0.025)  
quantile(trace0.delta(9, seq_post)', 0.975)  
% width of the credible interval
quantile(trace0.delta(9, seq_post)', 0.975) - quantile(trace0.delta(9, seq_post)', 0.025)

rng(0)
%% MNAR 
rand_theta1=(rand(1,length(seq_post))-0.5)*0.02; % draws of delta_mu
rand_theta2=(rand(1,length(seq_post))-0.5)*0.04; % draws of delta_beta
theta1_mnar1=trace0.alpha(9, seq_post)+rand_theta1*(1-197/310);
theta2_mnar1=trace0.delta(9, seq_post)+rand_theta2*(1-197/310);

mean(theta1_mnar1) 
% Credible interval, MNAR, theta_1
quantile(theta1_mnar1, 0.025) 
quantile(theta1_mnar1, 0.975)  
% width of the credible interval
quantile(theta1_mnar1, 0.975) - quantile(theta1_mnar1, 0.025)

mean(theta2_mnar1)
% Credible interval, MNAR, theta_2
quantile(theta2_mnar1, 0.025) 
quantile(theta2_mnar1, 0.975)
% width of the credible interval
quantile(theta2_mnar1, 0.975) - quantile(theta2_mnar1, 0.025)