% This file plots Supplementary figures 1-4 and saves them in separate pdf
% files within specified current directory. The current dirrectory must
% contain the directory "traces" with the MCMC traces.

%% read data
start_names_traces = 'traces/trace_2t202006_07583_standseed';
trace0=importdata(strcat(start_names_traces, '0.mat'));
trace1=importdata(strcat(start_names_traces, '30113.mat'));
trace2=importdata(strcat(start_names_traces, '104651.mat'));


seq_post=1000:5:5000; % iterations used for posterior estimation
seq_autocorr=1000:5000; % iterations used for autocorrelation computation
%% Supplemenatry Figure 1: convergence of alphas
h=figure('pos',[0 0 500 600]);
for i=1:8
    subplot(3,3,i);
    plot(trace0.alpha(i,:)'); hold on; 
    plot(trace1.alpha(i,:)'); 
    plot(trace2.alpha(i,:)'); 
    title(sprintf('\\alpha_{%d}',i));
    set(gca,'FontSize',12)
end
subplot(3,3,9);
plot(trace0.alpha(9,:)'); hold on; 
plot(trace1.alpha(9,:)'); 
plot(trace2.alpha(9,:)'); title('\theta_{\mu}');
set(gca,'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)-0.5, pos(4)-0.5])
print(h,'Sfig1_conv_alpha','-dpdf','-r500')
%% Supplementary Figure 2: convergence of deltas
h=figure('pos',[0 0 500 600]);
for i=1:8
    subplot(3,3,i);
    plot(trace0.delta(i,:)'); hold on; 
    plot(trace1.delta(i,:)'); 
    plot(trace2.delta(i,:)'); 
    title(sprintf('\\delta_{%d}',i));  set(gca,'FontSize',12)
end
subplot(3,3,9);
plot(trace0.delta(9,:)'); hold on; 
plot(trace1.delta(9,:)'); 
plot(trace2.delta(9,:)'); title('\theta_{\beta}');
set(gca,'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)-0.5, pos(4)-0.5])
print(h,'Sfig2_conv_delta','-dpdf','-r500')
%% Supplementary Figure 3: autocorrelation of alphas
h=figure('pos',[0 0 500 600]);
for i=1:8
    subplot(3,3,i);
    autocorr(trace0.alpha(i,seq_autocorr)');  title(sprintf('\\alpha_{%d}',i)); ylabel('');
    set(gca,'FontSize',12)
end
subplot(3,3,9);
autocorr(trace0.alpha(9,seq_autocorr)');  title('\theta_{\mu}');ylabel('');
set(gca,'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)-0.5, pos(4)-0.5])
print(h,'Sfig3_autocorr_alpha','-dpdf','-r0')
%% Supplementary Figure 4: autocorrelation delta
h=figure('pos',[0 0 500 600]);
for i=1:8
    subplot(3,3,i);
    autocorr(trace0.delta(i,seq_autocorr)');  title(sprintf('\\delta_{%d}',i)); ylabel('');
    set(gca,'FontSize',12)
end
subplot(3,3,9);
autocorr(trace0.delta(9,seq_autocorr)');  title('\theta_{\beta}');ylabel('');
set(gca,'FontSize',12)
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3)-0.5, pos(4)-0.5])
print(h,'Sfig4_autocorr_delta','-dpdf','-r0')




