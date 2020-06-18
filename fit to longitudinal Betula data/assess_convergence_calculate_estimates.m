% This file provides min and max scale reduction factor using seq_post
% iterations of MCMC as well estimates of the parameters using trace with seed 0.
% Note: theta is the last element of alphas and deltas.
% Input:
%   seq_post : MCMC iterations to use in te calculation of the scale reduction factor    
% Output: 
%   r: scale reduction factor for all the parameters in the model (min and max value and histogram)
% 
%% read the traces
% add the folder that contains "assess_convergence" and the folder "traces" with 3 traces to search path
% addpath('./traces')
start_names_traces = 'trace_2t202006_07583_standseed';
trace0=importdata(strcat(start_names_traces, '0.mat'));
trace1=importdata(strcat(start_names_traces, '30113.mat'));
trace2=importdata(strcat(start_names_traces, '104651.mat'));

%% transform the data into par_of_interest where each cell correspods to one parameter,
% and each cell contains 3*length(seq_post) matrix with seq_post MCMC iterations
% for 3 traces. 
seq_post=1000:5:5000;

fields=setdiff(fieldnames(trace0), {'a','d','N_c'});
par_of_interest=[];
names=[];
for k=1:length(fields)
    field0=trace0.(fields{k});
    field1=trace1.(fields{k});
    field2=trace2.(fields{k});
    for l=1:size(field0,1)
        q=length(par_of_interest);
        par_of_interest{q+1}=[field0(l,seq_post); field1(l,seq_post); field2(l,seq_post)];
        names{q+1}=strcat(fields{k}, sprintf( '_%d', l)); 
    end      
end

beta_zeta={'a','d'};
for k=1:length(beta_zeta)
    field0=reshape(cell2mat(trace0.(beta_zeta{k})(:,seq_post)),size(trace0.a{1},1),size(trace0.a{1},2), length(seq_post));
    field1=reshape(cell2mat(trace1.(beta_zeta{k})(:,seq_post)),size(trace0.a{1},1),size(trace0.a{1},2), length(seq_post));
    field2=reshape(cell2mat(trace2.(beta_zeta{k})(:,seq_post)),size(trace0.a{1},1),size(trace0.a{1},2), length(seq_post));
    
    for l=1:size(field0,1)
        for m=1:size(field0,2)
        q=length(par_of_interest);
        par_of_interest{q+1}=[reshape(field0(l,m,:),1, length(seq_post));...
            reshape(field1(l,m,:), 1,length(seq_post));...
            reshape(field2(l,m,:),1, length(seq_post))];
        names{q+1}=strcat(beta_zeta{k}, sprintf( 'cov%d_sub%d', l,m)); 
        end
    end      
end

%% Calculate scale reduction factor r
r=[]; 
for parameter=par_of_interest
    parameter_matrix=cell2mat(parameter);
    n=size(parameter_matrix,2);
    m=size(parameter_matrix,1);
    psi_dot_j=mean(parameter_matrix,2);
    b=n*var(psi_dot_j);
    w=1/m*sum(var(parameter_matrix,0,2));
    r=cat(1,r,sqrt(((n-1)/n*w+1/n*b)/w));
end

%% convergence
max(r)
min(r)
histogram(r); title('Scale reduction factor')

%% Parameters' estimates
mean(trace0.alpha(:, seq_post),2) % estimates of alpha and theta_mu
mean(trace0.delta(:, seq_post),2) % estimates of delta and theta_beta
mean(trace0.gamma_a_2(:, seq_post),2) % estimates of gamma_a^2
mean(trace0.gamma_d_2(:, seq_post),2) % estimates of gamma_d^2
mean(trace0.gamma_0_2(:, seq_post),2) % estimates of gamma_0^2
min(mean(trace0.sigma_0_2(:, seq_post),2)) %  min of the estimates of sigma_i0^2
max(mean(trace0.sigma_0_2(:, seq_post),2)) %  min of the estimates of sigma_i0^2
mean(mean(trace0.sigma_0_2(:, seq_post),2)) % mean of the estimates of sigma_i0^2
min(mean(trace0.sigma_1_2(:, seq_post),2)) %  min of the estimates of sigma_i1^2
max(mean(trace0.sigma_1_2(:, seq_post),2)) %  max of the estimates of sigma_i1^2
mean(mean(trace0.sigma_1_2(:, seq_post),2)) %  mean of the estimates of sigma_i1^2
