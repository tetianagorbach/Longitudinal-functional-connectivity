% This file compares the baseline characteristics of returnees and
% dropouts.
% The code reads the file "returnees.txt" that contains a vector of 1 for a returnee and 0
% for dropouts.
 
%% how "returnees.txt" were constructed
% "ID_inclusion_times.xlsx"  that has 3 columns:
%   1st column is subjects id's repeated once for dropouts or twice for
%   returnees
%   2nd column is Time with values 0 and 5 corresponding to baseline and
%   follow-up
%   3rd column is "Inclusion" with value 1 if a subjects satisfies health
%   inclusion criteria and 0 othewise.
% add path to the file 'ID_inclusion_times.xlsx'
% baseline_characteristic = readtable('ID_inclusion_times.xlsx'); 
% subjects = baseline_characteristic(baseline_characteristic.Inclusion==1,:); % subset to participants that satisfy the inclusion criterion
% n_sub = size(unique(subjects.MR_subject_ID), 1); % number of subjects at baseline
% unique_ids = table2array(subjects(subjects.Time==0, 'MR_subject_ID')); % 
% %% define who is returnee and who is dropout
% returnee = zeros(n_sub, 1); % vector of 1 for a returnee and 0 for a dropout
% for i=1:size(unique_ids,1) 
%     returnee(i) = size(subjects(strcmp(subjects.MR_subject_ID,  unique_ids{i}),:),1)-1;
% end
% 
% fileID = fopen('returnees.txt','w');
% fprintf(fileID,'%d\n',returnee);
% fclose(fileID);  
%% read data
% read vector of returnees
returnees=readtable('returnees.txt', 'ReadVariableNames',0);   
% read baseline characteristics
% add path to "Age_Sex_FD_stand_cov.xlsx"
cov = readtable('Age_Sex_FD_stand_cov.xlsx', 'sheet', 'Sheet1','ReadVariableNames' ,false);
cov(:,'returnees') = returnees(:,1);

% proportion of returnees
sum(cov.returnees)
sum(1-cov.returnees)

% test if proportion of males is different between dropouts and returnees.
[tbl,chi2,p,labels] = crosstab(cov.returnees, cov.Var2);

descr = array2table(zeros(0,7));
descr.Properties.VariableNames = {'mean_returnees','mean_dropout', ...
           'sd_returnees','sd_dropout',...
           't', 'df', 'p'};

for i=1:(size(cov, 2)-1)
    descr(2,i) = array2table(mean(table2array(cov(cov.returnees==1, i))));
    descr(3,i) = array2table(mean(table2array(cov(cov.returnees==0, i))));
    descr(4,i) = array2table(std(table2array(cov(cov.returnees==1, i))));
    descr(5,i) = array2table(std(table2array(cov(cov.returnees==0, i))));
    [h,p,ci,stats] = ttest2(table2array(cov(cov.returnees==1, i)),table2array(cov(cov.returnees==0, i)),'Vartype','unequal');
    descr(6,i) = array2table(stats.tstat);
    descr(7,i) = array2table(stats.df);
    descr(8,i) = array2table(p);
end


writetable(descr, 'descriptives.xlsx', 'Sheet', 'Sheet1')
