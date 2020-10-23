function plot_coverage_sep(n_pairs, n_subjects, coverage,  title,i, ylimits)
% PLOT_COVERAGE_SEP plots coverage of credible intervals for design 1
% Input:
%   n_pairs : vector of number of node pairs in each trace
%   n_sub :  vector of number of subjects in each trace
% Output:
%   plot:  coverage vs. number of node pairs, lines correspond to different
%   number of subjects.
plot( n_pairs(n_subjects==10), coverage(n_subjects==10),'k->',...
      n_pairs(n_subjects==20), coverage(n_subjects==20), 'k-o',...
      n_pairs(n_subjects==50), coverage(n_subjects==50), 'k-*',...
      n_pairs(n_subjects==100), coverage(n_subjects==100), 'k-s','MarkerFaceColor','auto', 'Linewidth', 0.5, 'MarkerSize',7);   ylabel(title);
  if i==1 
  legend({'20 sub', '50 sub', '100 sub'},'Location', 'south' )
  end
  ylim(ylimits); 
end
