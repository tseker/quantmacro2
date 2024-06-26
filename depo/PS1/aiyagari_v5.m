function aggK = aiyagari_q5(r0)
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
   
maxkap = 26;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.1;                    % size of capital grid increments
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

%  initialize some variables
%

   
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum 

% Seed for the value function:
%v   = zeros(N,nkap);
if isempty(v) % when this is true, we are in the first iteration
   v = ((wage*repmat(s',1,nkap) + r0*repmat(kap,N,1)).^(1-mu)-1)/((1-mu)*(1-beta));
end % when not true, we use the last v of previous iteration as a guess

% Initialising matrices and distance
test = 10;
iter = 0;

% Solving the value function
while test > 0.01
    % Update the counter
    iter = iter+1;

    % Optimisation
    for j=1:N %each possible skill s_t
        for i=1:nkap % each asset level a_t
            % Upper limit for golden search
                upperK = min(maxkap, wage*s(j) + (1+r0)*kap(i) - 0.001);
            % Policy function: assets (and expected value)
                auxf = @(x) value_v5(x,i,j,r0,wage); % aux function for golden
                [aprime_VFI(j,i), tv(j,i), auxInterp] = ...
                    golden_v5(auxf,minkap,upperK);
            % We extract the nearest elements on the grid and the weights,
            % so that we do not need to get them again later
                aprime_L(j,i) = auxInterp{1};
                aprime_U(j,i) = auxInterp{2};
                aprime_W(j,i) = auxInterp{3};
        end
    end
    
    % Update variables
        test=max(max(abs(tv-v))) % distance
        v=tv; % Stored in a global, so we can use it in the next iteration on r
end

 
display('done with VF iteration')


% Compute the measure
    % Initialise auxiliary variables
        indOpt = zeros(nkap, N*nkap); % Indicator: =1 if optimal decision
        Q = []; % Transition matrix
        auxAprime = reshape(aprime_VFI',1,N*nkap); % now each column represents a combination of s and a
    % Create matrix with 1s in the position of optimal decision
        % - row indicates the index of optimal a' on the grid.
        % - column indicates the starting par (s,a) on the grid.
        % Due to interpolation, some a' are not on the grid. Use weights
        % for the nearest elements that belong to the grid.
            % Find nearest elements and their weights
                % [auxL, auxU, auxW] = interpWeights(auxAprime,kap);
                % No need to run the previous line, we already have them
                % (we stored them after golden search)
                auxL = reshape(aprime_L',1,N*nkap);
                auxU = reshape(aprime_U',1,N*nkap);
                auxW = reshape(aprime_W',1,N*nkap);
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
            toler=1e-5; % tolerance to stop the loop
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

% Other returns
    % Policy function on assets (as a global)
        kk = aprime_VFI;
    % Aggregate capital
        aggK = sum(kk(:).*dist(:));


display('done with measure')