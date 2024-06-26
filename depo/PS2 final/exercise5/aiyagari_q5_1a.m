function [meanK,meanN] = aiyagari_q4_1a(r0,optHoward)
% aiyagari.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s N prob b kk kap v dist phi B n_VFI c_VFI


%% Preliminaries

% write wage as a function of interest rate using Cobb-Douglas
% production function
wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));

% Capital Grid
maxkap = 26;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.1;                    % size of capital grid increments
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

% Functions
focn = @(sj,n) B*n.^(phi)-(wage*sj*((1+r0)*kap' + wage*sj*n - kap).^(-mu));
util = @(cons,n) (cons.^(1-mu))/(1-mu)-B*(n.^(1+phi))/(1+phi);


%% Endogenous labour supply
for j=1:N
    % Initialise variables
        % Bounds
        nlow = (kap-(1+r0)*kap')/(wage*s(j)) + 1e-4;
        nup = ones(nkap,nkap);
        % Loop criterium
        maxtestn = 1;
    % Loop
    while (maxtestn>0.001)
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
            nlow(1-nlow<1e-3) = 1;
            nup(1-nlow<1e-3) = 1;
            testn(1-nlow<1e-3) = 0;
        % Update max distance
            maxtestn = max(max(testn));
    end
    nopt(:,:,j)=max(n0,0);
    % Rows: a. Columns: a'. 3rd dimension: current productivity.
end


%% Value Function Iteration

%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum 

% Seed for the value function:
%v   = zeros(N,nkap);
if isempty(v) % when this is true, we are in the first iteration
   v = ((wage*repmat(s',1,nkap) + r0*repmat(kap,N,1)).^(1-mu)-1)/((1-mu)*(1-beta));
end % when not true, we use the last v of previous iteration as a guess

utilm = zeros(N,nkap,nkap);

for j=1:N%loop over each skill s_t
    for i=1:nkap% loop over each possible a_t            
        cons = s(j)*wage*nopt(i,:,j) + (1+r0)*kap(i) -kap;
        ut = util(cons,nopt(i,:,j));
        ut(cons<= 0) = -10000;      % infeasible consumption choice      
        utilm(i,:,j)=ut;
    end
end

% Initialising matrices and distance
tv = 189*zeros(N,nkap);
tdecis = tv; 
test = 10;
criter_V = 1e-7;
iter=0;

% Solving the value function
while test > 1e-3
    % Update the counter
    iter = iter+1;

    % Optimisation
    for j=1:N%each possible skill s_t
        for i=1:nkap% each asset level a_t             
            ut=utilm(i,:,j); 
            vint = ut + beta*prob(j,:)*v;  % utility given j, i, for all possible k' next period       
            [tv(j,i),tdecis(j,i)] = max(vint);   %what is the optimal asset choice on the grid
        end
    end

    % Policy functions
        % Savings
        kk = kap(tdecis);
        % Labour supply
            % Version with loops
%             for j=1:N
%                 for i=1:nkap
%                     n_VFI(j,i) = nopt(i,tdecis(j,i),j);
%                 end
%             end
            % Version without loops
            n_VFI = nopt(sub2ind(size(nopt), ...
                    reshape(repmat(1:nkap,N,1),1,N*nkap), ...
                    reshape(tdecis,1,N*nkap), ...
                    reshape(repmat(1:N,1,nkap),1,N*nkap)));
            n_VFI = reshape(n_VFI,N,nkap);
        % Consumption
        c_VFI = s' * wage .* n_VFI + (1+r0)*kap - kk;

    % Howard
    if optHoward==1 && iter>3
        % Initialise variable
        dVV=1;
        v1how=tv;
        % Loop
        while dVV>criter_V
            for j = 1:N
                for i = 1:nkap
                    v1how(j,i)=util(c_VFI(j,i),n_VFI(j,i)) + beta*prob(j,:)*v1how(:,tdecis(j,i));
                end
            end
            dVV=max(max(abs(v1how-tv)));
            tv=v1how;
        end
    end

    % Update variables
    test=max(max(abs(tv-v))); % distance
    v=tv; % Stored in a global, so we can use it in the next iteration on r
end

display('done with VF iteration')


%% Measure

% Initialise auxiliary variables
    indOpt = sparse(nkap, N*nkap); % Indicator: =1 if optimal decision
    Q = []; % Transition matrix
    auxDecis = reshape(tdecis',1,N*nkap); % now each column represents a combination of s and a
% Create matrix with one in the position of optimal decision
    % - row indicates the index of optimal a' on the grid.
    % - column indicates the starting par (s,a) on the grid.
    indOpt(sub2ind(size(indOpt), auxDecis, 1:(N*nkap))) = 1;
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
    meanK = sum(kk(:).*dist(:));
    meanN = sum(n_VFI(:).*dist(:));