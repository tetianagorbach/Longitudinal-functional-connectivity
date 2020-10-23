function plot_coverage_real(n_pairs, n_subjects, mse,  title,i,ylimits)
% PLOT_COVERAGE_REAL plots coverage of credible intervals for design 2.
% Input:
%   n_pairs : vector of number of node pairs in each trace
%   n_sub :  vector of number of subjects in each trace
% Output:
%   plot:  coverage vs. number of node pairs, lines correspond to different
%   number of subjects.
plot( n_pairs(n_subjects==10), mse(n_subjects==10),'k->',...
      n_pairs(n_subjects==20), mse(n_subjects==20), 'k-o',...
      n_pairs(n_subjects==50), mse(n_subjects==50), 'k-*',...
      n_pairs(n_subjects==100), mse(n_subjects==100), 'k-s','MarkerFaceColor','auto', 'Linewidth', 0.5, 'MarkerSize',10);   ylabel(title);
  ylim(ylimits)
  if i==1
    legend({'20 sub', '50 sub', '100 sub'},'Location', 'south' )
  end
end
