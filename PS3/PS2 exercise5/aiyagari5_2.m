function [aggK,aggN] = aiyagari5_2(r0)
% aiyagari.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s b N prob b kk kap v dist phi B n_pol c_pol


%   write wage as a function of interest rate using Cobb-Douglas
%   production function

% r = r0;

wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));

%--------------------------------------------------------------------                               
                                
                              
%   form capital grid
%   
   
% Capital Grid
maxkap = 26;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.1;                    % size of capital grid increments
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

% Functions
focn = @(sj,n) B*n.^(phi)-(wage*sj*((1+r0)*kap' + wage*sj*n - kap).^(-mu));

%% Endogenous labour supply
% n(a,a',epsilon)
for j=1:N
    % Initialise variables
        % Bounds
        nlow = (kap-(1+r0)*kap')/(wage*s(j)) + 1e-4;
        nup = ones(nkap,nkap);
        % Loop criterium
        maxtestn = 1;
    % Loop
    while (maxtestn>1e-3)
        % Update guess
            n0 = 0.5*(nup-nlow)+nlow;
        % Evaluate the FOC
            f=focn(s(j),n0);
        % Test distance: how far are we from 0 (satisfied FOC)
            testn=abs(f);
        % Update bounds
            nup(f>0)=n0(f>0);
            nlow(f<0)=n0(f<0);
        % If nlow goes to one, we are looking at a' too high
        % at n = 1, foc could be zero
        % at n=0, foc is always negative, so don't need to check for
        % the boundary condtion.
            nlow(1-nlow<1e-3) = 1;
            nup(1-nlow<1e-3) = 1;
            testn(1-nlow<1e-3) = 0;
        % Update max distance
            maxtestn = max(max(testn));
    end
    nopt(:,:,j)=max(n0,0);
    % Rows: a. Columns: a'. 3rd dimension: current productivity.
end


%% Golden Search

% Seed for the value function:
%v   = zeros(N,nkap);
if isempty(v) % when this is true, we are in the first iteration
   v = ((wage*s' + r0*kap).^(1-mu)) / ((1-mu)*(1-beta));
end % when not true, we use the last v of previous iteration as a guess

% Initialising matrices and distance
test = 10;
iter = 0;

% Solving the value function
while test > 1e-3
    % Update the counter
    iter = iter+1;

    % Optimisation
    for j=1:N %each possible skill s_t
        % Golden search bounds
            upperK = min(maxkap, wage*s(j)+(1+r0)*kap - 0.001);
            lowerK = minkap*ones(1,length(upperK));
        % Policy function: assets and labour supply (and expected value)
            [kk(j,:), n_pol(j,:), tv(j,:), indU(j,:)] = ...
                golden5_2(r0,wage*s(j),lowerK,upperK,prob(j,:),nopt(:,:,j));
    end
    
    % Update variables
        test=max(max(abs(tv-v))); % distance
        v=tv; % Stored in a global, so we can use it in the next iteration on r
end

display('done with golden search')


%% Measure

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

display('done with measure')


%% Other returns
% Aggregate capital
    aggK = sum(kk(:).*dist(:));
    aggN = sum(n_pol(:).*dist(:));

% Policy function: consumption
    c_pol = s' * wage .* n_pol + (1+r0)*kap - kk;