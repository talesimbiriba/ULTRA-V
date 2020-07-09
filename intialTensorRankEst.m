function [ rank ] = intialTensorRankEst( T , min_diff )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <2 
    min_diff = 0.15;
end

if size(T,3) == 1
    
    error('Error: Input must be a tensor!')
    
elseif size(T,4) == 1
    % abundance tensor:  T is  IxJxK
    disp('Abundance Tensor')
    [I,J,K] = size(T);
    T_mat = zeros(I,J*K);
    count = 1; 
    for j=1:1:J, 
        for k=1:K, 
            T_mat(:,count) = T(:,j,k); 
            count = count+1; 
        end
    end
    s1 = svd(T_mat);
%     aa = [s1 ([1:size(s1)]*0.1)'];
%     [~,idx] = min(diag(aa*aa'));
%     rank = idx;
% rank = knee_pt(s1);
     rank = find(abs(diff(s1,1)) < min_diff);
     rank = rank(1);
     
%     plot(aa(:,2),aa(:,1));
%     hold on;
%     plot(aa(idx,2),aa(idx,1),'ro');
    
else
    % endmember tensor:  T is  IxJxKxL
    disp('Endmember Tensor')
    [I,J,K,L] = size(T);
    T_mat = zeros(I,J*K*L);
    count = 1; 
    for j=1:1:J, 
        for k=1:K, 
            for l=1:L
                T_mat(:,count) = T(:,j,k,l); 
                count = count+1; 
            end
        end
    end
    s1 = svd(T_mat);
    rank = find(abs(diff(s1,1)) < min_diff);
    rank = rank(1);
%     aa = [s1 [1:size(s1)]'];
%     [~,idx] = min(diag(aa*aa'));
%     rank = idx;        
end

end

