function [ A_est ] = tensorRegAbundanceEstimation(Y,M,Q, lambda_A, gamma)
% function [ A_est ] = tensorRegAbundanceEstimation(R,M,Q, lambda_A, gamma)
% 
% Fully constrained abundance estimation using low-rank tensor
% regularization
%
%   Y - N1xN2xL tensor (hypercube)
%   M - N1xN2xLxR tensor with one endmember matrix for pixel
%   Q - N1xN2xR low rank tensor related to the abundance tensor A
%
%
% Author: Tales Imbiriba
% Date: 30/11/2017

if nargin < 5
    gamma = 100;
end
    
% [N1,N2,L] = size(Y);
[N1,N2,L,R] = size(M);
A_est= zeros(N1,N2,R);

ones_R = ones(R,1);

%temp = lambda_A*eye(R) - gamma*(ones_R*ones_R');
temp = lambda_A*eye(R);
L = -eye(R);
k = zeros(R,1);
A = ones(1,R);
b = 1;
%options = optimset('TolX',1e-18,'LargeScale','off');

if (isempty(gcp('nocreate')))
    %disp('seq')
    for n1=1:N1
        for n2=1:N2
            H = squeeze(M(n1,n2,:,:))'*squeeze(M(n1,n2,:,:)) + temp;
            f = -squeeze(M(n1,n2,:,:))'*squeeze(Y(n1,n2,:)) -lambda_A*squeeze(Q(n1,n2,:));% -gamma*ones_R;   
            A_est(n1,n2,:) = qpas(H,f,L,k,A,b);                 
            %A_est(n1,n2,:) = qpas(H,f,L,k);
        end
    end
else
    %disp('parallel')
    parfor n1=1:N1
        for n2=1:N2
            H = squeeze(M(n1,n2,:,:))'*squeeze(M(n1,n2,:,:)) + temp;
            f = -squeeze(M(n1,n2,:,:))'*squeeze(Y(n1,n2,:)) -lambda_A*squeeze(Q(n1,n2,:));% -gamma*ones_R;   
            A_est(n1,n2,:) = qpas(H,f,L,k,A,b);                 
            %A_est(n1,n2,:) = qpas(H,f,L,k);
        end
    end
end
% 
% for n=1:N
%     
%     %f = (-Y(:,n)'*M)';
%     a_est(:,n) = qpas(H,f,L,k,A,b);                 
%     %a_est(:,n) = quadprog(H,f,L,k,[],[],[],[],[],options);
% end

%a_est  = (hyperFcls(Y,M));

end

