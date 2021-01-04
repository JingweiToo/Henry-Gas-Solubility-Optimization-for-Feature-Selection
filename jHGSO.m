
function [sFeat,Sf,Nf,curve] = jHGSO(feat,label,N,max_Iter,...
  num_clus,K,alpha,beta,L1,L2,L3,HO)

% Parameters
lb      = 0;
ub      = 1; 
thres   = 0.5; 

Ttheta  = 298.15;
eps     = 0.05; 
c1      = 0.1;
c2      = 0.2; 
 
% Objective function
fun = @jFitnessFunction;
% Number of dimensions
dim = size(feat,2); 
% Number of gas in cluster
Nn  = ceil(N / num_clus); 
% Initial 
X   = zeros(N,dim); 
for i = 1:N
	for d = 1:dim
    X(i,d) = lb + (ub - lb) * rand();
	end
end
% Henry constant & E/R constant 
H = zeros(num_clus,1); 
C = zeros(num_clus,1); 
P = zeros(num_clus,Nn);
for j = 1:num_clus
  H(j) = L1 * rand();
  C(j) = L3 * rand();
  for i = 1:Nn
    % Partial pressure 
    P(j,i) = L2 * rand();
  end
end
% Divide the population into gas/cluster type of gas cluster
Cx = cell(num_clus,1); 
for j = 1:num_clus
  if j ~= num_clus
    Cx{j} = X(((j - 1) * Nn) + 1 : j * Nn, :);
  else
    Cx{j} = X(((num_clus - 1) * Nn + 1 : N), :);
  end
end
% Fitness of each cluster
Cfit  = cell(num_clus,1); 
fitCB = ones(1,num_clus); 
Cxb   = zeros(num_clus,dim); 
fitG  = inf;
for j = 1:num_clus
  for i = 1:size(Cx{j},1)
    Cfit{j}(i,1) = fun(feat,label,(Cx{j}(i,:) > thres),HO);
    % Update best gas
    if Cfit{j}(i) < fitCB(j)
      fitCB(j) = Cfit{j}(i);
      Cxb(j,:) = Cx{j}(i,:);
    end
    % Update global best
    if Cfit{j}(i) < fitG
      fitG = Cfit{j}(i);
      Xgb  = Cx{j}(i,:);
    end
  end
end
% Pre 
S     = zeros(num_clus,Nn); 

curve = zeros(1,max_Iter);
curve(1) = fitG;
t = 2; 
% Iterations
while t <= max_Iter
  % Compute temperature 
  T = exp(-t / max_Iter); 
  for j = 1:num_clus
    % Update henry coefficient 
    H(j) = H(j) * exp(-C(j) * ((1 / T) - (1 / Ttheta)));
    for i = 1:size(Cx{j},1)
      % Update solubility 
      S(j,i) = K * H(j) * P(j,i);
      % Compute gamma 
      gamma  = beta * exp(-((fitG + eps) / (Cfit{j}(i) + eps)));
      % Flag change between - & +
      if rand() > 0.5
        F = -1;
      else
        F = 1; 
      end
      for d = 1:dim     
        % Random constant
        r = rand();
        % Position update 
        Cx{j}(i,d) = Cx{j}(i,d) + F * r * gamma * ...
          (Cxb(j,d) - Cx{j}(i,d)) + F * r * alpha * ...
          (S(j,i) * Xgb(d) - Cx{j}(i,d));
      end
      % Boundary
      XB = Cx{j}(i,:); XB(XB > ub) = ub; XB(XB < lb) = lb;
      Cx{j}(i,:) = XB;
    end
  end
  % Fitness
  for j = 1:num_clus
    for i = 1:size(Cx{j},1)
      % Fitness
      Cfit{j}(i,1) = fun(feat,label,(Cx{j}(i,:) > thres),HO); 
    end
  end
  % Select the worst solution 
  Nw = round(N * (rand() * (c2 - c1) + c1));
  % Convert cell to array
  XX = cell2mat(Cx); 
  FF = cell2mat(Cfit);
  [~, idx] = sort(FF,'descend');
  % Update position of worst solution 
  for i = 1:Nw
    for d = 1:dim
      XX(idx(i),d) = lb + rand() * (ub - lb);
    end
    % Fitness
    FF(idx(i)) = fun(feat,label,(XX(idx(i),:) > thres),HO);
  end
  % Divide the population into gas type of gas cluster back
  for j = 1:num_clus
    if j ~= num_clus
      Cx{j}   = XX(((j - 1) * Nn) + 1 : j * Nn, :); 
      Cfit{j} = FF(((j - 1) * Nn) + 1 : j * Nn);
    else
      Cx{j}   = XX(((num_clus - 1) * Nn + 1 : N), :); 
      Cfit{j} = FF((num_clus - 1) * Nn + 1 : N);
    end
  end  
  % Update best solution
  for j = 1:num_clus
    for i = 1:size(Cx{j},1)
      % Update best gas
      if Cfit{j}(i) < fitCB(j)
        fitCB(j) = Cfit{j}(i);
        Cxb(j,:) = Cx{j}(i,:);
      end
      % Update global best
      if Cfit{j}(i) < fitG
        fitG = Cfit{j}(i);
        Xgb  = Cx{j}(i,:);
      end
    end
  end
  curve(t) = fitG; 
  fprintf('\nIteration %d Best (HGSO)= %f',t,curve(t))
  t = t + 1;
end
% Select features
Pos   = 1:dim; 
Sf    = Pos((Xgb > thres) == 1);
sFeat = feat(:,Sf);
Nf    = length(Sf); 
end


