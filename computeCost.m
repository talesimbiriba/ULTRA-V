function [ cost ] = computeCost(Y,M,A,P,Q,lambda_M,lambda_A, verbose)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


[N1,N2,~] = size(Y);

c0 = 0;
for n1=1:N1
    for n2=1:N2
       c0 =  c0 + frob(squeeze(Y(n1,n2,:))- squeeze(M(n1,n2,:,:))*squeeze(A(n1,n2,:)),'squared');
    end
end

c1 = lambda_A*frob(A-Q,'squared') ;
c2 = lambda_M*frob(M-P,'squared');


% cost = 0.5*c0 + 0.5* (lambda_A*frob(A-Q)  + lambda_M*frob(M-P));

cost  = 0.5*(c0 + c1 + c2);
if verbose 
    fprintf('Error Norm = %f, low-rank Anorm = %f, low-rank Mnorm = %f\n', c0, c1,c2);
    fprintf('Global cost = %f', cost);
end
end

