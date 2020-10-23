% This file contains the code for Figure 1, Figure 2, and Figure 3
%% Figure1: histograms of the generated data
% histogram of separated data
h = figure('pos', [0 0 500 200]);
i = 1;
sim_res_gen = importdata('../4. simulation study/design 1 results/sim_2tseed0sub100pairs1000.mat');
edges = [-2 -1:0.25:10 10];
subplot('Position', [0.05 0.1 0.44 0.8]); 
histogram(sim_res_gen.gen_data_1.gen_data(1, :), 'FaceColor', 'none', 'BinEdges', edges); title('Design 1'); hold on; % generated data for the first subject, 
% that has generated both baseline and follow up ( 1st and 2nd element in sim_res_gen.gen_data_1.ind.obs are 1 meaning that 1st and 2nd orw of gen_data cooresponds to generated subject 1) 
histogram(sim_res_gen.gen_data_1.gen_data(2, :), 'FaceColor', 'black', 'FaceAlpha', 0.5 , 'BinEdges', edges); 
% histogram of real data
sim_res_real = importdata('../4. simulation study/design 2 results/sim_2trealseed0sub50pairs1000.mat');
edges = [-1 -1:0.1:1.5 1.5];
i = 1;
subplot('Position', [0.54 0.1 0.44 0.8]); histogram(sim_res_real.gen_data_1.gen_data(1, :), edges, 'FaceColor', 'none'); title('Design 2'); hold on;
histogram(sim_res_real.gen_data_1.gen_data(2, :), edges, 'FaceColor', 'black', 'FaceAlpha', 0.5); 
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(h, 'Figure1', '-dpdf', '-r500')

%% Figure 2: coverage for clearly separated mixtures 
addpath('../4. simulation study/design 1 results/')
start_names_traces = 'sim_2tseed0sub';
sim_res_path = {strcat(start_names_traces, '20pairs100.mat'), strcat(start_names_traces, '20pairs400.mat'),  strcat(start_names_traces, '20pairs800.mat'), strcat(start_names_traces, '20pairs1000.mat'), strcat(start_names_traces, '20pairs1500.mat'),...
                strcat(start_names_traces, '50pairs100.mat'), strcat(start_names_traces, '50pairs400.mat'),  strcat(start_names_traces, '50pairs800.mat'), strcat(start_names_traces, '50pairs1000.mat'), strcat(start_names_traces, '50pairs1500.mat'),...
                strcat(start_names_traces, '100pairs100.mat'), strcat(start_names_traces, '100pairs400.mat'), strcat(start_names_traces, '100pairs800.mat'), strcat(start_names_traces, '100pairs1000.mat'), strcat(start_names_traces, '100pairs1500.mat')};
seq = 100:1000;
delta_mu = -0.02;
delta_beta = 0.04;
rng(0);
for i = 1:length(sim_res_path)
    sim_res{i} =  importdata(sim_res_path{i});
    n_subjects(i) =  sim_res{i}.gen_data_1.n_sub;
    rand_theta1{i} =  rand(length(seq), size(sim_res{i}.results.phi_trace, 2));
    rand_theta2{i} =  rand(length(seq), size(sim_res{i}.results.phi_trace, 2));
end
for i =  1:length(sim_res_path)
    n_subjects(i) =  sim_res{i}.gen_data_1.n_sub;
    n_pairs(i) =  sim_res{i}.gen_data_1.n_pairs;
    theta1_mar_l(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :), 0.025);
    theta1_mar_u(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :), 0.975);
    theta1_mnar_l(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta1{i}*(delta_mu), 0.025);
    theta1_mnar_u(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta1{i}*(delta_mu), 0.975);
    theta2_mar_l(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :), 0.025);
    theta2_mar_u(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :), 0.975);
    theta2_mnar_l(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta2{i}*(delta_beta), 0.025);
    theta2_mnar_u(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta2{i}*(delta_beta), 0.975);
end
theta1 =  sim_res{1}.gen_data_1.alpha(9);
theta1_mnar =  theta1+(1 - sim_res{1}.gen_data_1.phi)*(delta_mu);
theta2 =  sim_res{1}.gen_data_1.delta(9);
theta2_mnar =  theta2+(1 - sim_res{1}.gen_data_1.phi)*(delta_beta);

h =  figure('pos', [0 0 500 500]);
% theta_1
coverage_mar =  mean((theta1 > theta1_mar_l).* (theta1 < theta1_mar_u ), 2);
coverage_data_mnar_est_mar =  mean((theta1_mnar > theta1_mar_l).* (theta1_mnar < theta1_mar_u ), 2);
coverage_data_mnar_est_mnar =  mean((theta1_mnar > theta1_mnar_l).* (theta1_mnar < theta1_mnar_u ), 2);

subplot(2, 3, 1); plot_coverage_sep(n_pairs, n_subjects, coverage_mar, '\theta_1', 1, [0.2 1]);  
title('MAR, MAR'); set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
subplot(2, 3, 2); plot_coverage_sep(n_pairs, n_subjects, coverage_data_mnar_est_mar, '', 0, [0.2 1]); 
title('MNAR, MAR'); set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
subplot(2, 3, 3); plot_coverage_sep(n_pairs, n_subjects, coverage_data_mnar_est_mnar, '', 0, [0.2 1]); 
title('MNAR, MNAR'); set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
% theta_2
coverage_mar =  mean((theta2 > theta2_mar_l).* (theta2 < theta2_mar_u ), 2);
coverage_data_mnar_est_mar =  mean((theta2_mnar > theta2_mar_l).* (theta2_mnar < theta2_mar_u ), 2); % data MNAR; but MAR in the estimation
coverage_data_mnar_est_mnar =  mean((theta2_mnar > theta2_mnar_l).* (theta2_mnar < theta2_mnar_u ), 2); % data MNAR; and MNAR in the estimation
% plot
subplot(2, 3, 4); plot_coverage_sep(n_pairs, n_subjects, coverage_mar, '\theta_2', 0, [0.2 1]); 
set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
subplot(2, 3, 5); plot_coverage_sep(n_pairs, n_subjects, coverage_data_mnar_est_mar, '', 0, [0.2 1]); 
set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
subplot(2, 3, 6); plot_coverage_sep(n_pairs, n_subjects, coverage_data_mnar_est_mnar, '', 0, [0.2 1]); 
set(gca, 'FontSize', 12); xlabel('Number of pairs'); xlim([100 1500]);
set(h,'Units', 'Inches');

pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)-0.5])
print(h, 'Figure2', '-dpdf', '-r500')

%% Figure 3: Coverage for Design 2.
addpath('../4. simulation study/design 2 results/')
start_names_traces = 'sim_2trealseed0sub';
sim_res_path =  {strcat(start_names_traces , '20pairs500.mat'), strcat(start_names_traces, '20pairs1000.mat'), strcat(start_names_traces, '20pairs1500.mat'), strcat(start_names_traces, '20pairs2000.mat'), strcat(start_names_traces, '20pairs3000.mat'), strcat(start_names_traces, '20pairs5000.mat'), strcat(start_names_traces, '20pairs10000.mat'),...
              strcat(start_names_traces , '50pairs500.mat'), strcat(start_names_traces, '50pairs1000.mat'), strcat(start_names_traces, '50pairs1500.mat'), strcat(start_names_traces, '50pairs2000.mat'), strcat(start_names_traces, '50pairs3000.mat'), strcat(start_names_traces, '50pairs5000.mat'), strcat(start_names_traces, '50pairs10000.mat'),...
              strcat(start_names_traces , '100pairs500.mat'), strcat(start_names_traces, '100pairs1000.mat'), strcat(start_names_traces, '100pairs1500.mat'), strcat(start_names_traces, '100pairs2000.mat'), strcat(start_names_traces, '100pairs3000.mat'), strcat(start_names_traces, '100pairs5000.mat'), strcat(start_names_traces, '100pairs10000.mat')};

seq =  3000:5:6000;
delta_mu =  -0.02;
delta_beta =  0.04;
rng(1077281);
for i =  1:length(sim_res_path)
    sim_res{i} =  importdata(sim_res_path{i});
    rand_theta1{i} =  rand(length(seq), size(sim_res{i}.results.phi_trace, 2));
    rand_theta2{i} =  rand(length(seq), size(sim_res{i}.results.phi_trace, 2));
end
for i =  1:length(sim_res_path)
    n_subjects(i) =  sim_res{i}.gen_data_1.n_sub;
    n_pairs(i) =  sim_res{i}.gen_data_1.n_pairs;
    theta1_mar_l(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :), 0.025);
    theta1_mar_u(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :), 0.975);
    theta1_mnar_l(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta1{i}*(delta_mu), 0.025);
    theta1_mnar_u(i, :) =  quantile(sim_res{i}.results.alpha_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta1{i}*(delta_mu), 0.975);
    theta2_mar_l(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :), 0.025);
    theta2_mar_u(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :), 0.975);
    theta2_mnar_l(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta2{i}*(delta_beta), 0.025);
    theta2_mnar_u(i, :) =  quantile(sim_res{i}.results.delta_diff_trace(seq, :)+...
        (1-sim_res{i}.results.phi_trace(seq, :)).*rand_theta2{i}*(delta_beta), 0.975);
end
theta1 =  sim_res{1}.gen_data_1.alpha(9);
theta1_mnar =  theta1+(1-sim_res{1}.gen_data_1.phi)*(delta_mu);
theta2 =  sim_res{1}.gen_data_1.delta(9);
theta2_mnar =  theta2+(1-sim_res{1}.gen_data_1.phi)*(delta_beta);
h =  figure('pos', [0 0 500 500]);

 q_l =  []; q_u =  [];
% theta1
coverage_mar =  mean((theta1>theta1_mar_l).* (theta1<theta1_mar_u ), 2);
coverage_data_mnar_est_mar =  mean((theta1_mnar>theta1_mar_l).* (theta1_mnar<theta1_mar_u ), 2);
coverage_data_mnar_est_mnar =  mean((theta1_mnar>theta1_mnar_l).* (theta1_mnar<theta1_mnar_u ), 2);

subplot(2, 3, 1); plot_coverage_real(n_pairs, n_subjects, coverage_mar, '\theta_1', 1, [0.2 1]);  
title('MAR, MAR');set(gca,'FontSize', 12); xlabel('Number of pairs'); xlim([500 10000]);
subplot(2, 3, 2); plot_coverage_real(n_pairs, n_subjects, coverage_data_mnar_est_mar, '', 0, [0.2 1]); 
title('MNAR, MAR');set(gca,'FontSize', 12);xlabel('Number of pairs');xlim([500 10000]);
subplot(2, 3, 3); plot_coverage_real(n_pairs, n_subjects, coverage_data_mnar_est_mnar, '', 0, [0.2 1]); 
title('MNAR, MNAR');set(gca,'FontSize', 12); xlabel('Number of pairs');xlim([500 10000]);
% theta2
coverage_mar =  mean((theta2 > theta2_mar_l).* (theta2 < theta2_mar_u ), 2);
coverage_data_mnar_est_mar =  mean((theta2_mnar > theta2_mar_l).* (theta2_mnar < theta2_mar_u ), 2);
coverage_data_mnar_est_mnar =  mean((theta2_mnar > theta2_mnar_l).* (theta2_mnar < theta2_mnar_u ), 2);

subplot(2, 3, 4); plot_coverage_real(n_pairs, n_subjects, coverage_mar, '\theta_2', 0, [0.2 1]); 
set(gca,'FontSize', 12); xlabel('Number of pairs'); xlim([500 10000]);
subplot(2, 3, 5); plot_coverage_real(n_pairs, n_subjects, coverage_data_mnar_est_mar, '', 0, [0.2 1]); 
set(gca,'FontSize', 12);xlabel('Number of pairs'); xlim([500 10000]);
subplot(2, 3, 6); plot_coverage_real(n_pairs, n_subjects, coverage_data_mnar_est_mnar, '', 0, [0.2 1]); 

set(gca, 'FontSize', 12); xlabel('Number of pairs');
set(h, 'Units', 'Inches');

pos = get(h,'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(h, 'Figure3', '-dpdf', '-r0')


