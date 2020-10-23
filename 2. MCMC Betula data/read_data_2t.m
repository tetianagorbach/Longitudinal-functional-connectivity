function z = read_data_2t(data_path)
    % READ_DATA_2T reads the matrix of Fisher-transformed correlations 
    %   between node pairs for all subjects from data_path
    %   and transforms it to [number of observations (n_1 + 2*(n-n_1)) ,  number of node pairs]
    %   matrix.
    % Input: 
    %   data_path : path to the file with the matrix of Fisher-transformed  correlations 
    %   of size [number of nodes,  number of nodes, n_1 + 2*(n-n_1)].
    % Output: 
    %   z : matrix of size [n_1 + 2*(n-n_1),  number of node pairs]
    %   of the Fisher-transformed correlations.

    data = importdata(data_path);    
    [n, n, k] = size(data);
    z = zeros(k, n * (n - 1)/2);
    l = 0;
    
    for i = 1:k % for each subject and observation
        l = l + 1;
        v = triu(data(:,:,i))'; % transform data(:,:,i) to upper triangular. Important for the order of node pairs in the final data. 
        z(l,:) = v(v~=0); % transform upper triangular matrix to a vector.
    end  
end