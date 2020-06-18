function z1=read_data_2_t(data_path)
    
    data=importdata(data_path);
    % transform cor from a matrix of pairwise z-transformed correlations to upper triangular (since symmetrical) and then to a vector
    % z1 is a matrix of number of subjects*number of node pairs
    % real data
    if strfind(data_path,'matrixOfCorrelations_ztranform_wo270_272')>0
        [n,n,k]=size(data);
        z1=zeros(k,n*(n-1)/2);
        l=0;
        for i=1:k
            l=l+1;
            v=triu(data(:,:,i))';
            z1(l,:)=v(v~=0);
        end
    else 
        z1=data.z_ij;
    end    
end