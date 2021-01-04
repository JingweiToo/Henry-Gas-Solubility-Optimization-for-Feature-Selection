# Henry Gas Solubility Optimization for Feature Selection

![Wheel](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/cb2ed62c-670f-4619-bd1c-3764b2d36ad2/cd253c4b-2281-4d19-9da3-9a50e2d26d99/images/1601813643.JPG)

## Introduction
* This toolbox offers a Henry Gas Solubility Optimization ( HGSO ) method
* The < Main.m file > illustrates the example of how HGSO can solve the feature selection problem using benchmark data-set.

## Input
* *feat*     : feature vector ( Instances *x* Features )
* *label*    : label vector ( Instances *x* 1 )
* *N*        : number of gas
* *max_Iter* : maximum number of iterations
* *num_clus* : Number of gas types
* *K*        : Constant
* *alpha*    : Influence of other gas
* *beta*     : Constant
* *L1*       : Initial parameter
* *L2*       : Initial parameter
* *L3*       : Initial parameter


## Output
* *sFeat*    : selected features
* *Sf*       : selected feature index
* *Nf*       : number of selected features
* *curve*    : convergence curve


### Example
```code
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
```

## Requirement
* MATLAB 2014 or above
* Statistics and Machine Learning Toolbox

