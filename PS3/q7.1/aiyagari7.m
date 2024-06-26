function K1 = aiyagari7(r0)
% aiyagari.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s N prob b kk kap v dist


%   write wage as a function of interest rate using Cobb-Douglas
%   production function

% r = r0;

wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));

%--------------------------------------------------------------------                               
                                
                              
%   form capital grid
%   
   
maxkap = 130;                    % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.1;                    % size of capital grid increments
nkap = 200;
loggrid = linspace(0,log(maxkap+1-minkap), nkap);
kap = exp(loggrid) - (1-minkap); % state of assets

%  initialize some variables
%

   
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum 

% Seed for the value function:
%v   = zeros(N,nkap);
if isempty(v) % when this is true, we are in the first iteration
   v = ((wage*s' + r0*kap).^(1-mu) -1) / ((1-mu)*(1-beta));
end % when not true, we use the last v of previous iteration as a guess

% Initialising matrices and distance
test = 10;
iter = 0;

% Cash on hand
coh = wage*s' + (1+r0)*kap;


%% GOLDEN SEARCH

% Convert variables in vectors of dimension (N*nkap)
%   *Note: each group of nkap rows has same productivity (i.e. from
%   one row to the next, a always varies, while s only changes once
%   every nkap rows).
    % Cash on hand
        coh_gd = reshape(coh',1,N*nkap);
    % Probability of transition to s'
        prob_gd = reshape(prob,N^2,1); % make prob column vector
        prob_gd = repmat(prob_gd,1,nkap); % replicate nkap columns
        prob_gd = reshape(prob_gd',nkap*N,N); % reshape: s' in columns

% Golden search bounds
    upperK = min(maxkap, coh_gd - 0.001);
    lowerK = minkap*ones(1,N*nkap);

% Solving the value function
while test > 0.01
    % Update the counter
        iter = iter+1;
    % Optimisation: policy functions and expected value
        [kk, tv, indU] = golden7(coh_gd,prob_gd*v,lowerK,upperK);
    % Recover matrix form
        tv=reshape(tv,nkap,N)';
    % Update variables
        test=max(max(abs(tv-v))); % distance
        v=tv; % Stored in a global, so we can use it in the next iteration on r
end

% Recover matrix form
    kk=reshape(kk,nkap,N)';
    indU=reshape(indU,nkap,N)';


display('done with golden search')


%% MESURE

% Initialise auxiliary variables
    indOpt = sparse(nkap, N*nkap); % Indicator: =1 if optimal decision
    Q = []; % Transition matrix
    auxAprime = reshape(kk',1,N*nkap); % now each column represents a combination of s and a
% Create matrix with 1s in the position of optimal decision
    % - row indicates the index of optimal a' on the grid.
    % - column indicates the starting par (s,a) on the grid.
    % Due to interpolation, some a' are not on the grid. Use weights
    % for the nearest elements that belong to the grid.
        % Find nearest elements and their weights
            % Upper neighbour: already found in golden search (indU)
            % Reshape: to convert in row vector
            auxU = reshape(indU',1,N*nkap);
            % Lower neighbour:
            auxL = auxU-1;
            % When a' is the first element in kap, the lower neighbour
            % will not matter, all the weight will be on upper. But we
            % need auxL to remain on the grid, so it cannot be 0. It
            % cannot be 1 either, or we will have auxL=auxU and the
            % weight will not be well defined.
            auxL(auxL==0) = 2;
            % Weight:
            auxW = (kap(auxU) - auxAprime) ./ (kap(auxU) - kap(auxL));
        % Assign weight to the lower neighbour
        indOpt(sub2ind(size(indOpt), auxL, 1:(N*nkap))) = auxW;
        % Assing weight to the upper neighbour
        indOpt(sub2ind(size(indOpt), auxU, 1:(N*nkap))) = 1-auxW;

% Multiplying by the probability of each state, we obtain the
% transition matrix (Q)
    for j=1:N % loop over each productivity state
        % We need to have a vector of probabilities with N*nkap
        % elements (it will be each of the pi(j | s0) repeated nkap
        % times).
            auxProb = reshape(transpose(repmat(prob(:,j),1,nkap)),1,N*nkap);
        % Multiply each column of indOpt by the corresponding
        % probability
            auxQ = indOpt .* auxProb;
        % Add these rows to the bottom of Q
            Q = [Q; auxQ];
    end
    % About the transition matrix Q:
        % - Columns indicate the starting combination of (s,a).
        % - Rows indicate the ending combination of (s',a').
        % Therefore, each cell indicates the probability of reaching
        % (s',a') when starting from (s,a).
% Now that we have Q, use it in a loop to obtain the measure
    % Initialise auxiliary variables
        toler=1e-8; % tolerance to stop the loop
        distance=1; % distance between old and new distribution
        mu0 = ones(N*nkap,1) / (N*nkap); % starting distribution
    % Get stationary distribution
        while distance>=toler
            mu1 = Q*mu0;
            distance=max(max(abs(mu1-mu0)));
            mu0 = mu1;
        end
    % Reshape distribution:
        % - each row represents a state s'.
        % - each column represents the decision on assets a'.
        dist = transpose(reshape(mu1,nkap,N));

display('done with measure')


%% RETURNS

% Aggregate capital
    K1 = sum(kk(:).*dist(:));