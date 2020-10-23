% subsets raw connectivity data to include subjects and timpeoints with Inclusion == 1 in 
% "ID_Inclusion_Age_Sex_FD_stand_cov.xlsx" constructed by
% "define_covariates.r".
FC = importdata('FC_Tanya.mat');
times=readtable('ID_Inclusion_Age_Sex_FD_stand_cov.xlsx');
data = FC(:,:,times.Inclusion==1);
% B0016.50 has only T5, B2502.65 has no cognition.
save("FC_Tanya_subbsetted_t5t6_inclusion1.mat",'data');