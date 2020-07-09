function [ M, A , P, Q] = unmixing_lr_reg(Y, lambda_M, lambda_A, Km, Ka, A0, M0, maxItt, verbose,varepsilon)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 9 
    verbose = true;
end

% Tensorlab opts:
% options.Display = true; % Show progress on the command line.
options.Display = false; % Show progress on the command line.
% options.Initialization = @cpd_rnd; % Select pseudorandom initialization.
options.Initialization = @cpd_gevd; % Select pseudorandom initialization.
% options.Algorithm = @cpd_als;
options.Algorithm = @cpd_nls; % Select ALS as the main algorithm.
options.AlgorithmOptions.LineSearch = @cpd_els; % Add exact line search.

%  options.AlgorithmOptions.TolFun = 1e-12; % Set function tolerance stop criterion
% options.AlgorithmOptions.TolX = 1e-12; % Set step size tolerance stop criterion

% options.AlgorithmOptions.TolFun = 1e-6; % Set function tolerance stop criterion
% options.AlgorithmOptions.TolX = 1e-6; % Set step size tolerance stop criterion

options.AlgorithmOptions.TolFun = 1e-5; % Set function tolerance stop criterion
options.AlgorithmOptions.TolX = 1e-5; % Set step size tolerance stop criterion

options.AlgorithmOptions.TolFun = 1e-4 % Set function tolerance stop criterion
options.AlgorithmOptions.TolX = 1e-4; % Set step size tolerance stop criterion


% varepsilon =.00001;
% varepsilon =1e-10;
if nargin < 10
    varepsilon =1e-3;
end
delta_cost = inf;
i = 1;
A = A0;
% [N1,N2,L] = size(Y);
% R = size(A,3);
% P = zeros([N1,N2,L,R]);
% Q = zeros(size(A0));

% M = zeros([N1,N2,L,R]);
% for n1=1:N1,
%     for n2=1:N2,
%         M(n1,n2,:,:) = M0 ;% + rand(size(M0))*0.0000001;       
%     end
% end
%P = cpdgen(cpd(M,Km,options));
M = M0;
% P = M;
% Q = A0;
% P = ones(size(M))*0.5;
P = cpdgen(cpd(M0 + randn(size(M0))*0.0001,Km,options));
Q = cpdgen(cpd(A0,Ka,options));

cost_old = inf;
cost_init = computeCost(Y,M,A,P,Q,lambda_M,lambda_A, verbose); 
% cost_init = 1;
while delta_cost/cost_init > varepsilon && i <= maxItt 
    if verbose
        fprintf('\nIteration #%d\n', i);
    end
    
    M = tensorRegEndmemberEstimation(Y,A,P, lambda_M);
    A = tensorRegAbundanceEstimation(Y,M,Q, lambda_A);
    
    % polyadic tensor decomposition
%     P = cpdgen((cpd_gevd(M,Km)));
%     Q = cpdgen((cpd_gevd(A,Ka)));
 

    P = cpdgen(cpd(M,Km,options));
    Q = cpdgen(cpd(A,Ka,options));
    

%     U0 = cpd_gevd(P,Km);
%     U0 = {U0{4},U0{3},U0{1},U0{2}};
%     P = cpdgen(cpd_nls(M,U0));
%     
%     U0 = cpd_gevd(Q,Ka);
%     Q = cpdgen(cpd_nls(A,U0));
    
    cost = computeCost(Y,M,A,P,Q,lambda_M,lambda_A, verbose);
%     if i==1, 
%         cost_init = cost;
%     end
    delta_cost = abs(cost_old - cost);
    cost_old = cost;
    if verbose
        fprintf(' delta_cost/cost_init = %f\n', delta_cost/cost_init);
    end
    
    i = i+1;
end

if verbose
    fprintf('Finished!\n')
end

end

