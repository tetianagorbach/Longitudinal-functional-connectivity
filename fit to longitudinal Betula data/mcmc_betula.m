function trace = long_real(n_it,data_path, cov_path, set_rng,name_trace, start_values) 
    %read data    
    post_inv_gamma = @ (a_prior, n_used, b_prior, data_used) 1/gamrnd(a_prior+n_used/2,(1/b_prior+data_used'*data_used/2)^(-1));
    rng(set_rng);
    %% read data
    y=read_data_2_t(data_path);       
    n_sub=size(y,1)/2; % number of subjects 
    n_pairs=size(y,2); % number of node pairs
    y=reshape(y', n_pairs*n_sub*2,1); % long vector of observations, sub1,t1; sub1,t2; sub2,t2, sub2,t2; ...

    n_obs=2*n_sub; %number of observed subjects
    ind_obs=repelem(1:n_sub,1,2); %indeces of observed subjects
    %% read covariates for fixed effects
    cov=readtable(cov_path, 'sheet', 'Sheet1','ReadVariableNames' ,false);
    cov.Var1=zscore(str2double(cov.Var1)); % standardize t1 age
    cov.Var3=zscore(cov.Var3);% FD
    cov.Var4=zscore(cov.Var4);% EM
    cov.Var5=zscore(cov.Var5);% Fl
    cov.Var6=zscore(cov.Var6);% PS
    cov.Var7=zscore(cov.Var7);% Bl
        
    x_init=ones(n_obs,1);   
    x_init(:,2)=cov.Var1(ind_obs); 
    x_init(:,3)=cov.Var2(ind_obs); 
    x_init(:,4)=cov.Var3(ind_obs);
    x_init(:,5)=cov.Var4(ind_obs);
    x_init(:,6)=cov.Var5(ind_obs);
    x_init(:,7)=cov.Var6(ind_obs);
    x_init(:,8)=cov.Var7(ind_obs);
    x_init(2:2:size(x_init,1),9)=1;
    
    x=repelem(x_init, n_pairs,1);    
    % covariates for random effects 
    x_random = x(:,1); % just intercept
    dim_cov_pop=size(x,2);
    dim_cov_ind=size(x_random,2);      
   
    a2_short=zeros(n_sub, n_sub);
    for i=1:n_obs
        a2_short(i,ind_obs(i))=1;
    end
    
    %% initial values of parameters
    % parameters for priors
    a_sigma_1 = 1.5;  b_sigma_1  = 1000;
    a_sigma_0 = 1.5;  b_sigma_0  = 1000;
    % parameters for hyperpriors
    a_gamma = 1.5;    b_gamma = 1000;

    trace=[];
    % strength of connections
    trace.alpha(:,1) = [start_values(1);start_values(2);zeros(dim_cov_pop-2,1)];
    trace.gamma_a_2(:,1) = 1./gamrnd(a_gamma*ones(dim_cov_ind,1),b_gamma*ones(dim_cov_ind,1));    
    trace.a{1} =  mvnrnd(zeros(dim_cov_ind,n_sub),reshape(repelem(trace.gamma_a_2(1), 1,n_sub), 1,n_sub,dim_cov_ind));    
    % non-connected component
    trace.gamma_0_2(:,1) = 1/gamrnd(a_gamma,b_gamma,1,1);
    trace.mu_0(:,1) = mvnrnd(zeros(1,n_sub), trace.gamma_0_2(1)*ones(1,n_sub));
    % variances in normals
    trace.sigma_1_2(:,1) = 1./gamrnd(a_sigma_1,b_sigma_1,n_sub,1);
    trace.sigma_0_2(:,1) = 1./gamrnd(a_sigma_0,b_sigma_0,n_sub,1);
    % probability to be in connected component
    trace.delta(:,1) = [start_values(3);start_values(4);zeros(dim_cov_pop-2,1)];
    trace.gamma_d_2(:,1) = 1./gamrnd(a_gamma*ones(dim_cov_ind,1),b_gamma*ones(dim_cov_ind,1));   
    trace.d{1} = mvnrnd(zeros(dim_cov_ind,n_sub),reshape(repelem(trace.gamma_d_2(1), 1,n_sub), 1,n_sub,dim_cov_ind));
    %% posterior estimation
    for l=2:n_it
        l
        % initializing values to use
        % strength
        alpha=trace.alpha(:, l-1);
        gamma_a_2 = trace.gamma_a_2(l-1);
        a = trace.a{l-1};
        a_rep=repelem(a(ind_obs)',n_pairs,1);
        gamma_0_2 = trace.gamma_0_2(l-1);
        mu_0 = trace.mu_0(:,l-1);
        sigma_1_2 = trace.sigma_1_2(:,l-1);
        sigma_1_2_long=repelem(sigma_1_2,2*n_pairs,1);
        sigma_0_2 = trace.sigma_0_2(:,l-1);
        sigma_0_2_long=repelem(sigma_0_2,2*n_pairs,1);
        % probability
        delta = trace.delta(:,l-1);
        gamma_d_2=trace.gamma_d_2(l-1);        
        d = trace.d{l-1};
        d_rep=repelem(d(ind_obs)',n_pairs,1);

        lambda=normcdf(x*delta+sum(x_random.*d_rep,2)); % probability for each of 2*nsub*npairs pairs being connected
        mu_ij = x*alpha+sum(x_random.*a_rep,2); % strength of connections for connected component

        %% update of mixture component indicators
        p_nom=lambda.* lognpdf(y, mu_ij,  sqrt(sigma_1_2_long));
        omega=binornd(1, p_nom./(p_nom+(1-lambda) .* normpdf(y, repelem(mu_0(ind_obs),n_pairs,1),  sqrt(sigma_0_2_long))));

        %%  update hyperparameters 
        for k=1:dim_cov_ind
          gamma_a_2(k) = post_inv_gamma(a_gamma,n_sub,b_gamma,a(k,:)');
          gamma_d_2(k) = post_inv_gamma(a_gamma,n_sub,b_gamma,d(k,:)');
        end         
        gamma_0_2(1)  = post_inv_gamma(a_gamma,n_sub,b_gamma,mu_0);
      

        %% update parameters for means
        N_c = sum(reshape(omega,n_pairs, n_obs))'; % number of connected pairs by person
        ln_y_ij_c=log(y.*omega); % ln_y_ij_c=ln(y) for those pairs that are connected, -Inf otherwise
        ln_y_ij_c(ln_y_ij_c==-Inf)=0; % ln_y_ij_c=ln(y) for those pairs that are connected, 0 otherwise
        
        N_nc=n_pairs-N_c;
        y_nc=reshape(y.*(1-omega),n_pairs,n_obs);
        psi_mu_0=zeros(n_sub,1);
        mu_0=zeros(n_sub,1);
        for i=1:n_sub 
            psi_mu_0(i)=(1/gamma_0_2+sum(N_nc(ind_obs==i))/sigma_0_2(i))^(-1);
            mu_0(i)=normrnd(sum(sum(y_nc(:,ind_obs==i)))/sigma_0_2(i)*psi_mu_0(i), sqrt(psi_mu_0(i)));
        end

        % define z_probit
        mean_prob=x*delta+sum(x_random.*d_rep,2);
        f_limit=normcdf(-mean_prob);
        z_probit=rand(length(omega),1);
        con = omega==1;
        z_probit(con)=f_limit(con)+z_probit(con).*(1-f_limit(con));
        z_probit(~con)=z_probit(~con).*f_limit(~con);
        z_probit=mean_prob+norminv(z_probit);

        % update alpha, a
        sigma_1_2_obs=repelem(sigma_1_2,2,1);
        a1psia1=(x_init./sigma_1_2_obs.*N_c)'*x_init+diag(repelem(1/10^3,dim_cov_pop,1));
        a1psia2=(x_init./sigma_1_2_obs.*N_c)'*a2_short;
        a2psia2=(a2_short./sigma_1_2_obs.*N_c)'*a2_short+diag(repmat(1./gamma_a_2,n_sub,1));
        dthetastar=[a1psia1,a1psia2;a1psia2',a2psia2]^(-1);
        dthetastar=(dthetastar+dthetastar')/2;
        etheta=dthetastar*[(x_init./sigma_1_2_obs)'*sum(reshape(ln_y_ij_c,n_pairs,n_obs),1)';...
                          (a2_short./sigma_1_2_obs)'*sum(reshape(ln_y_ij_c,n_pairs,n_obs),1)'];
        theta=mvnrnd(etheta', dthetastar);

        alpha=theta(1:dim_cov_pop)';
        a=theta(dim_cov_pop+1:length(theta));   
        a_rep=repelem(a(ind_obs)',n_pairs,1);

        % update delta,d
        a1psia1=(x_init*n_pairs)'*x_init+diag(repelem(1/10^3,dim_cov_pop,1));
        a1psia2=(x_init*n_pairs)'*a2_short;
        a2psia2=(a2_short*n_pairs)'*a2_short+diag(repmat(1./gamma_d_2,n_sub,1));
        dthetastar=[a1psia1,a1psia2;a1psia2',a2psia2]^(-1);
        dthetastar=(dthetastar+dthetastar')/2;
        etheta=dthetastar*[(x_init)'*sum(reshape(z_probit,n_pairs,n_obs),1)';...
                           a2_short'*sum(reshape(z_probit,n_pairs,n_obs),1)' ];
        theta=mvnrnd(etheta', dthetastar);

        delta=theta(1:dim_cov_pop)';
        d=theta(dim_cov_pop+1:length(theta));
        %% update variances
        mu_ij_c = (x*alpha+sum(x_random.*a_rep,2)).*omega; % strength of connections for "connected pairs".
        one_minus_omega_prod_z=(1-omega).*y;

         for i=1:n_sub 
            start_mat=2*(i-1)*n_pairs+1;
            end_mat=2*i*n_pairs;
            sigma_1_2(i) = post_inv_gamma(a_sigma_1, N_c(2*i-1)+N_c(2*i), b_sigma_1, ...
                ln_y_ij_c(start_mat:end_mat)-mu_ij_c(start_mat:end_mat)); 
            sigma_0_2(i) = post_inv_gamma(a_sigma_0,(2*n_pairs - N_c(2*i-1)-N_c(2*i)),...
                b_sigma_0, one_minus_omega_prod_z(start_mat:end_mat)-(1-omega(start_mat:end_mat))*mu_0(i));          
        end
        
        %% storage of updates
        trace.alpha(:,l)=alpha;
        trace.a{l}=a;
        trace.gamma_0_2(:,l)=gamma_0_2;
        trace.mu_0(:,l)= mu_0;
        trace.sigma_1_2(:,l)=sigma_1_2;
        trace.sigma_0_2(:,l)=sigma_0_2;
        trace.gamma_a_2(:,l)=gamma_a_2;
        trace.delta(:,l)= delta;
        trace.d{l}=d;        
        trace.gamma_d_2(:,l)=gamma_d_2;
        trace.N_c(:,l)=N_c;  
        if floor(l/100)==l/100
            l
            save(strcat(char(name_trace), sprintf('seed%d', set_rng)),'-struct','trace');
        end
    end
end
