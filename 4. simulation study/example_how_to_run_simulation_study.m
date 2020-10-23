% Example on how to run MCMC locally and on High Performance Computing Center North
%% Locally
% Note: cd should include sim_study_2t.m, mcmc_2t.m, generate_data_real.m,
% and generate_data_sep.m
n_mcmc_it = 20;
n_sets_gen = 3;
trace_for_data_simulation_path = 'trace_2t_with_missing_seed0.mat';
seq_post = 1000:5:6000;
cov_path = 'Age_Sex_FD_stand_cov.xlsx';
set_rng = 0;
name_results = 'zzz';
phi = 0.6;
n_sub_vect = 50;
n_pairs_vect = 100;
bb=sim_study_2t(n_mcmc_it, n_sets_gen, trace_for_data_simulation_path, ...
    seq_post, cov_path, set_rng, name_results,phi,n_sub_vect, n_pairs_vect,0);

%% Terminal code to run line by line on High Performance Computing Center North (access rights are required)
ssh xxx
start_matlab % custom function that starts MATLAB on HPC2N
configCluster
c = parcluster('kebnekaise')
c.AdditionalProperties.AccountName = 'SNICyear-x-xxx'
c.AdditionalProperties.WallTime = '100:00:00'
c.saveProfile
% separated components
j=c.batch(@sim_study_2t, 1,{1000,1000, 'trace_2t_with_missing_seed0.mat', 1000:5:6000, 'Age_Sex_FD_stand_cov.xlsx',  0, 'sim_2t',0.6,[20,50,100],[100,400,800,1000,1500],1},'pool',11)
% less separated components
j=c.batch(@sim_study_2t, 1,{6000, 200, 'trace_2t_with_missing_seed0.mat', 1000:5:6000, 'Age_Sex_FD_stand_cov.xlsx',  0, 'sim_2treal',0.6,[20,50,100],[500, 1000, 1500, 2000, 3000, 5000, 10000],0},'pool',11)