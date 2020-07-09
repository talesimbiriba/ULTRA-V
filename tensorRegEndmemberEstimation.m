function [ M_est ] = tensorRegEndmemberEstimation(Y,A,P, lambda_M)
% function [ M_est ] = tensorRegAbundanceEstimation(Y,A,P, lambda_M)
% 
% constrained endmember estimation using low-rank tensor
% regularization
%
%   Y - N1xN2xL tensor (hypercube)
%   A - N1xN2xR tensor with one endmember matrix for pixel
%   P - N1xN2xLxR low rank tensor related to the abundance tensor A
%
%
% Author: Tales Imbiriba
% Date: 1/12/2017



[N1,N2,L,R] = size(P);
M_est = zeros(size(P));
Ir = eye(R);
parfor n1=1:N1
    for n2=1:N2
        M_est(n1,n2,:,:) = (squeeze(Y(n1,n2,:))*squeeze(A(n1,n2,:))' + lambda_M*squeeze(P(n1,n2,:,:))) ...
                /( squeeze(A(n1,n2,:))*squeeze(A(n1,n2,:))' + lambda_M*Ir);
        M_est(n1,n2,:,:) = max(M_est(n1,n2,:,:), 1e-6);
    end
end

end

