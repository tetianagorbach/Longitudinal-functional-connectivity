function trace = mcmc_2t_with_dropouts(n_it, seed)
    % MCMC_2T_WITH_DROPOUTS Runs MCMC for the Betula data
    % Input: 
    %   n_it: number of MCMC iterations
    %   seed for random number generation.
    %   The following files are also used:
    %   "FC_Tanya_subbsetted_t5t6_inclusion1.mat",
    %   "Age_Sex_FD_stand_cov.xlsx", and
    %   "ID_inclusion_times.xlsx" (columns MR.subject.ID, Time, Inclusion from 
    %   "ID_Inclusion_Age_Sex_FD_stand_cov.xlsx" constructed 
    %   during preprocessing to define what rows correspond to returnees and what to dropouts.)
    % Output:
    %   file "trace_2t_with_missing_seed(insert seed used here).mat" that contains  MCMC draws.
    
    times = readtable('ID_inclusion_times.xlsx');
    subjects = times(times.Inclusion==1,:);

    post_inv_gamma = @ (a_prior, n_used, b_prior, data_used) 1/gamrnd(a_prior+n_used/2, (1/b_prior+data_used'*data_used/2)^(-1));
    rng(seed);
    y = read_data_2t('FC_Tanya_subbsetted_t5t6_inclusion1.mat');       
    n_sub = size(unique(subjects.MR_subject_ID), 1); % number of subjects 
    n_pairs = size(y,2); % number of node pairs
    n_obs = size(y, 1); % number of observed subjects


    % define starting value for alpha
    top_10 = zeros(size(y, 1), ceil(n_pairs/10));
    for i = 1:size(y, 1)
        m = sort(y(i, :), 'descend');
        top_10(i, :) = log(m(1:ceil(n_pairs/10)));
    end
    start_for_alpha = mean(mean(top_10));

    y = reshape(y', size(y,1)*size(y,2), 1); % long vector of observations, sub1,t1; sub1,t2; sub2,t2, sub2,t2; ...

    returnee = zeros(n_sub, 1);
    unique_ids = table2array(subjects(subjects.Time==0, 'MR_subject_ID'));

    for i=1:size(unique_ids,1) 
        returnee(i) = size(subjects(strcmp(subjects.MR_subject_ID,  unique_ids{i}),:),1)-1;
    end
    ind_obs =  repelem((1:n_sub)', returnee+1); % vector of subject id repeated 1 time for dropouts and 2 times for returnees
    t =  ones(n_obs, 1); % t =  1 if the first timepoint, t =  2 if the second timepoint
        for i =  2:n_obs
            if ind_obs(i) == ind_obs(i - 1)
                t(i) = 2;
            end
        end    

    % read covariates for fixed effects
    cov = readtable('Age_Sex_FD_stand_cov.xlsx', 'sheet', 'Sheet1','ReadVariableNames' ,false);
    x_init_preproc = ones(n_sub,1);
    x_init_preproc(:,2) = zscore(cov.Var1); % standardize t1 age
    x_init_preproc(:,3) = cov.Var2; % Gender
    x_init_preproc(:,4) = zscore(cov.Var3);%  FD
    x_init_preproc(:,5) = zscore(cov.Var4);% EM
    x_init_preproc(:,6) = zscore(cov.Var5);% Fl
    x_init_preproc(:,7) = zscore(cov.Var6);% PS
    x_init_preproc(:,8) = zscore(cov.Var7);% Bl

    x_init = x_init_preproc(ind_obs,:); % covariates for all observed subjects
    x_init(t == 2, 9) = 1;    
    x = repelem(x_init, n_pairs,1);    

    % covariates for random effects 
    x_random = x(:,1); 
    dim_cov_pop = size(x,2);
    dim_cov_ind = size(x_random,2);

    a2_short =  zeros(n_sub, n_sub);
        for i =  1:n_obs
            a2_short(i, ind_obs(i)) =  1;
        end

    %% starting values of parameters
        % parameters for priors
        a_sigma_1 = 1.5;  b_sigma_1  = 1000;
        a_sigma_0 = 1.5;  b_sigma_0  = 1000;
        % parameters for hyperpriors
        a_gamma = 1.5;       b_gamma = 1000;

        trace =  [];
        % strength of connections
        trace.alpha(:, 1) = [start_for_alpha; zeros(dim_cov_pop-1, 1)];
        trace.gamma_a_2(:, 1) = 1./gamrnd(a_gamma*ones(dim_cov_ind, 1), b_gamma*ones(dim_cov_ind, 1));    
        trace.a{1} =  mvnrnd(zeros(dim_cov_ind, n_sub), reshape(repelem(trace.gamma_a_2(:, 1), 1, n_sub), 1, n_sub, dim_cov_ind));    
        % non-connected component
        trace.gamma_0_2(:, 1) = 1./gamrnd(a_gamma, b_gamma, 1, 1);
        trace.mu_0(:, 1) = mvnrnd(zeros(1, n_sub), trace.gamma_0_2(1)*ones(1, n_sub));
        % variances in normals
        trace.sigma_1_2(:, 1)  =  1./gamrnd(a_sigma_1, b_sigma_1, n_sub, 1);
        trace.sigma_0_2(:, 1) = 1./gamrnd(a_sigma_0, b_sigma_0, n_sub, 1);
        % number of connected edges
        trace.delta(:, 1) = [-1.28 ;zeros(dim_cov_pop-1, 1)];
        trace.gamma_d_2(:, 1) = 1./gamrnd(a_gamma*ones(dim_cov_ind, 1), b_gamma*ones(dim_cov_ind, 1));   
        trace.d{1} = mvnrnd(zeros(dim_cov_ind, n_sub), reshape(repelem(trace.gamma_d_2(:, 1), 1, n_sub), 1, n_sub, dim_cov_ind));
        trace.phi(1:n_it) =  betarnd(1+n_obs-n_sub, 1+2*n_sub-n_obs, 1, n_it);
        %% posterior estimation
        for l =  2:n_it
            %% initialize values to use
            % strength
            alpha =  trace.alpha(:, l-1);
            gamma_a_2 = trace.gamma_a_2(l-1);
            a = trace.a{l-1};
            a_rep =  repelem(a(ind_obs)', n_pairs, 1);
            gamma_0_2 = trace.gamma_0_2(l-1);
            mu_0 = trace.mu_0(:, l-1);
            sigma_1_2 = trace.sigma_1_2(:, l-1);
            sigma_1_2_ind_obs =  sigma_1_2(ind_obs, :);
            sigma_1_2_long =  repelem(sigma_1_2_ind_obs, n_pairs, 1);
            sigma_0_2 = trace.sigma_0_2(:, l-1);
            sigma_0_2_ind_obs =  sigma_0_2(ind_obs, :);
            sigma_0_2_long =  repelem(sigma_0_2_ind_obs, n_pairs, 1);
            % probability
            delta = trace.delta(:, l-1);
            gamma_d_2 =  trace.gamma_d_2(l-1);        
            d = trace.d{l-1};
            d_rep =  repelem(d(ind_obs)', n_pairs, 1);

            lambda =  normcdf(x*delta+sum(x_random.*d_rep, 2)); % probability for each of 2*nsub*npairs pairs being connected
            mu_ij = x*alpha+sum(x_random.*a_rep, 2); % strength of connections for connected component

            %% update mixture component indicators
            p_nom =  lambda.* lognpdf(y, mu_ij,  sqrt(sigma_1_2_long));
            omega =  binornd(1, p_nom./(p_nom+(1-lambda) .* normpdf(y, repelem(mu_0(ind_obs), n_pairs, 1),  sqrt(sigma_0_2_long))));
            %%  update hyperparameters 
            for k =  1:dim_cov_ind
              gamma_a_2(k) = post_inv_gamma(a_gamma, n_sub, b_gamma, a(k, :)');
              gamma_d_2(k) = post_inv_gamma(a_gamma, n_sub, b_gamma, d(k, :)');
            end         
            gamma_0_2 = post_inv_gamma(a_gamma, n_sub, b_gamma, mu_0);


            %% update parameters for means
            N_c = sum(reshape(omega, n_pairs, n_obs))'; % number of connected pairs by person
            ln_y_ij_c =  log(y.*omega); % ln_y_ij_c =  ln(y) for those pairs that are connected, -Inf otherwise
            ln_y_ij_c(ln_y_ij_c == -Inf) = 0; % ln_y_ij_c = ln(y) for those pairs that are connected, 0 otherwise

            N_nc =  n_pairs-N_c;
            y_nc =  reshape(y.*(1-omega), n_pairs, n_obs);
            psi_mu_0 =  zeros(n_sub, 1);
            mu_0 =  zeros(n_sub, 1);
            for i =  1:n_sub 
                psi_mu_0(i) =  (1/gamma_0_2+sum(N_nc(ind_obs == i))/sigma_0_2(i))^(-1);
                mu_0(i) =  normrnd(sum(sum(y_nc(:, ind_obs == i)))/sigma_0_2(i)*psi_mu_0(i), sqrt(psi_mu_0(i)));
            end

            % define z_probit
            mean_prob =  x*delta+sum(x_random.*d_rep, 2);
            f_limit =  normcdf(-mean_prob);
            z_probit =  rand(length(omega), 1);
            con = omega == 1;
            z_probit(con) =  f_limit(con)+z_probit(con).*(1-f_limit(con));
            z_probit(~con) =  z_probit(~con).*f_limit(~con);
            z_probit =  mean_prob+norminv(z_probit);

            % update alpha, a
            a1psia1 =  (x_init./sigma_1_2(ind_obs).*N_c)'*x_init+diag(repelem(1/10^3, dim_cov_pop, 1));
            a1psia2 =  (x_init./sigma_1_2(ind_obs).*N_c)'*a2_short;
            a2psia2 =  (a2_short./sigma_1_2(ind_obs).*N_c)'*a2_short+diag(repmat(1./gamma_a_2, n_sub, 1));
            dthetastar =  [a1psia1, a1psia2; a1psia2', a2psia2]^(-1);
            dthetastar =  (dthetastar+dthetastar')/2;
            etheta =  dthetastar*[(x_init./sigma_1_2(ind_obs))'*sum(reshape(ln_y_ij_c, n_pairs, n_obs), 1)';...
                                (a2_short./sigma_1_2(ind_obs))'*sum(reshape(ln_y_ij_c, n_pairs, n_obs), 1)'];
            theta =  mvnrnd(etheta', dthetastar);

            alpha =  theta(1:dim_cov_pop)';
            a =  theta(dim_cov_pop+1:length(theta));   
            a_rep =  repelem(a(ind_obs)', n_pairs, 1);

            % update delta, d
            a1psia1 =  (x_init*n_pairs)'*x_init+diag(repelem(1/10^3, dim_cov_pop, 1));
            a1psia2 =  (x_init*n_pairs)'*a2_short;
            a2psia2 =  (a2_short*n_pairs)'*a2_short+diag(repmat(1./gamma_d_2, n_sub, 1));
            dthetastar =  [a1psia1, a1psia2; a1psia2', a2psia2]^(-1);
            dthetastar =  (dthetastar+dthetastar')/2;
            etheta =  dthetastar*[(x_init)'*sum(reshape(z_probit, n_pairs, n_obs), 1)';...
                                a2_short'*sum(reshape(z_probit, n_pairs, n_obs), 1)' ];
            theta =  mvnrnd(etheta', dthetastar);

            delta =  theta(1:dim_cov_pop)';
            d =  theta(dim_cov_pop+1:length(theta));
            %% update variances
            mu_ij_c  =   reshape((x*alpha+sum(x_random.*a_rep, 2)).*omega, n_pairs, n_obs); % strength of connections for "connected pairs".
            ln_y_ij_c_resh =  reshape(ln_y_ij_c, n_pairs, n_obs);
            one_minus_omega_resh =  (1-reshape(omega, n_pairs, n_obs));
            for i =  1:n_sub 
                sigma_1_2(i) = post_inv_gamma(a_sigma_1, sum(N_c(ind_obs == i)), b_sigma_1, ...
                    reshape(ln_y_ij_c_resh(:, ind_obs == i) - mu_ij_c(:, ind_obs == i), sum(ind_obs == i)*n_pairs, 1)); 
                sigma_0_2(i) = post_inv_gamma(a_sigma_0, sum(N_nc(ind_obs == i)), ...
                    b_sigma_0,  reshape(y_nc(:, ind_obs == i)-one_minus_omega_resh(:, ind_obs == i)*mu_0(i), sum(ind_obs == i)*n_pairs, 1));          
            end

            %% store updates
            trace.alpha(:, l) =  alpha;
            trace.a{l} =  a;
            trace.gamma_0_2(:, l) = gamma_0_2;
            trace.mu_0(:, l) = mu_0;
            trace.sigma_1_2(:, l) =  sigma_1_2;
            trace.sigma_0_2(:, l) =  sigma_0_2;
            trace.gamma_a_2(:, l) =  gamma_a_2;
            trace.delta(:, l) = delta;
            trace.d{l} =  d;        
            trace.gamma_d_2(:, l) =  gamma_d_2;
            trace.N_c(:, l) =  N_c;  
            if floor(l/100)==l/100
                save(strcat(sprintf('trace_2t_with_missing_seed%d', seed)),'-struct','trace');
            end
        end
end
