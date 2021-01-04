%---------------------------------------------------------------------%
%  Henry Gas Solubility Optimization (HGSO) source codes demo version %
%---------------------------------------------------------------------%


%---Inputs------------------------------------------------------------
% feat     : features
% label    : labelling
% N        : Number of gases
% max_Iter : Maximum number of iterations
% num_clus : Number of gas types
% K        : Constant
% alpha    : Influence of other gas
% beta     : Constant
% L1       : Initial parameter
% L2       : Initial parameter
% L3       : Initial parameter

%---Outputs-----------------------------------------------------------
% sFeat    : Selected features
% Sf       : Selected feature index
% Nf       : Number of selected features
% curve    : Convergence curve
%---------------------------------------------------------------------


%% Henry Gas Solubility Optimization
clc, clear, close; 
% Benchmark data set 
load ionosphere.mat; 

% Set 20% data as validation set
ho = 0.2; 
% Hold-out method
HO = cvpartition(label,'HoldOut',ho);

% Parameter setting
N        = 10; 
max_Iter = 100; 
num_clus = 2;      % number of gas types / cluster
K        = 1;      % constant
alpha    = 1;      % influence of other gas
beta     = 1;      % constant 
L1       = 5E-3; 
L2       = 100; 
L3       = 1E-2;

% Henry Gas Solubility Optimization
[sFeat,Sf,Nf,curve] = jHGSO(feat,label,N,max_Iter,num_clus,K,alpha,beta,L1,L2,L3,HO);

% Plot convergence curve
plot(1:max_Iter,curve);
xlabel('Number of iterations');
ylabel('Fitness Value'); 
title('HGSO'); grid on;




