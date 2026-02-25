function [Lambda1, Lambda2, Lambda3] = eig3volume_matlab(Hxx, Hxy, Hxz, Hyy, Hyz, Hzz)
% MATLAB implementation of eig3volume MEX function
% Computes the eigenvalues of 3x3 symmetric matrices for each voxel
%
% Input: Hessian matrix components for each voxel
%   Hxx, Hxy, Hxz, Hyy, Hyz, Hzz - vectors of matrix elements
%
% Output: 
%   Lambda1, Lambda2, Lambda3 - eigenvalues sorted in descending order
%   (Lambda1 >= Lambda2 >= Lambda3)

% Convert to column vectors if needed
Hxx = Hxx(:);
Hxy = Hxy(:);
Hxz = Hxz(:);
Hyy = Hyy(:);
Hyz = Hyz(:);
Hzz = Hzz(:);

n = length(Hxx);

% For better performance with large arrays, process in chunks
chunk_size = 10000;
num_chunks = ceil(n / chunk_size);

Lambda1 = zeros(n, 1);
Lambda2 = zeros(n, 1); 
Lambda3 = zeros(n, 1);

for chunk = 1:num_chunks
    start_idx = (chunk - 1) * chunk_size + 1;
    end_idx = min(chunk * chunk_size, n);
    chunk_indices = start_idx:end_idx;
    chunk_size_actual = length(chunk_indices);
    
    % Extract chunk data
    Hxx_chunk = Hxx(chunk_indices);
    Hxy_chunk = Hxy(chunk_indices);
    Hxz_chunk = Hxz(chunk_indices);
    Hyy_chunk = Hyy(chunk_indices);
    Hyz_chunk = Hyz(chunk_indices);
    Hzz_chunk = Hzz(chunk_indices);
    
    % Process each element in the chunk
    for i = 1:chunk_size_actual
        % Construct 3x3 symmetric matrix
        H = [Hxx_chunk(i), Hxy_chunk(i), Hxz_chunk(i);
             Hxy_chunk(i), Hyy_chunk(i), Hyz_chunk(i);
             Hxz_chunk(i), Hyz_chunk(i), Hzz_chunk(i)];
        
        % Compute eigenvalues
        eigenvals = sort(eig(H), 'descend');
        
        % Assign to output (sorted in descending order)
        actual_idx = start_idx + i - 1;
        Lambda1(actual_idx) = eigenvals(1);
        Lambda2(actual_idx) = eigenvals(2);
        Lambda3(actual_idx) = eigenvals(3);
    end
    
    % Display progress for large datasets
    if n > 50000 && mod(chunk, max(1, floor(num_chunks/10))) == 0
        fprintf('Processing eigenvalues: %d%% complete\n', round(100 * chunk / num_chunks));
    end
end

end 