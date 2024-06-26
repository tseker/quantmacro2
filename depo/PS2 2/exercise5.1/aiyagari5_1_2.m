
function aggregateK = aiyagari5_1_2(r0,optHow)
% aiyagari.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s N prob b kk kap v dist phi B n_VFI c_VFI


%   write wage as a function of interest rate using Cobb-Douglas
%   production function

% r = r0;

wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));

%--------------------------------------------------------------------                               
                                
                              
%   form capital grid
%   
   
maxkap = 50;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.015;                    % size of capital grid increments
kap    = log(minkap+b+0.1):inckap:log(maxkap+b+0.1);   % state of assets 
kap = exp(kap)-b-0.1;
nkap   = length(kap);            % number of grid points
B = 80;
phi    = 2.5;             % Frisch elasticity


%  initialize some variables
criterion = 10;
iter = 0;

%util and foc_n
focn = @(sj,n) B*n.^(phi)-(wage*sj*((1+r0)*kap' + wage*sj*n - kap).^(-mu));
util = @(cons,n) (cons.^(1-mu))/(1-mu)-B*(n.^(1+phi))/(1+phi);


   
%%

% Vectorize to create matrix of states for each (a, epsilon):
    % Create indexes for each:
     %    ind.z= 1;       % Productivity
     %    ind.a = 2;      % Savings
     %    ind.a0 = find (kap>=0);  % Find positions for zero asset
     % 
     % % Create a matrix
     %  mat.ea(:, ind.z) = kron((1:N)', ones(nkap,1));
     %  mat.ea(:, ind.a) = kron(ones(N,1),(1:nkap)');


%% Endogenous Labor Supply

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
    nopt(j,:)=max(n0,0);
end



%% VFI

utilm = 190*zeros(N,nkap,nkap);

for j=1:N%loop over each skill s_t
    for i=1:nkap% loop over each possible a_t            
        cons = s(j)*wage*nopt(j,i,:) + (1+r0)*kap(i) -kap';
        util = util(cons,nopt(j,i,:));
        util(cons<= 0) = -10000;      % infeasible consumption choice      
        utilm(j,i,:)=util;
    end
end

% Initialising matrices and distance
tv = 189*zeros(N,nkap);
tdecis = tv; 
test = 10;
iter=0;

tdecis = zeros(N,nkap);
aprime_VFI = zeros(N,nkap);
%Solve the value function
while test > 0.001
    iter = iter + 1;
    % Interpolation 
    %W(epsilon,a') = beta* \sum_{epsilon'} Prob{epsilon'|epsilon} V(epsilon',a')
    % Only Interpolate W on the asset grid.
    W_matrix = beta * prob * v;
    for j = 1:N
        % use griddedInterpolant: only need to interpolate once.
        W_func = griddedInterpolant(kap,W_matrix(j,:));
        % vectorize the loop over the asset grid
        auxf = @(x) newV(x, W_func, j,wage, r0 );
        upperK = min(max(kap), wage*s(j) + (1+r0)*kap - 0.0001);
        [aprime_VFI(j,:), tv(j,:)] = ...
                        golden(auxf,minkap,upperK);
            % We extract the nearest elements on the grid and the weights,
            for i = 1:nkap
            [~,tdecis(j,i)]= min(abs(aprime_VFI(j,i) - kap));
            end
    end


    c_VFI = s'*wage*nopt + (1+r0)*kap - aprime_VFI;



    test=max(max(abs(tv-v))); % distance
    v=tv;
end



 %%

% Compute the measure
    % Initialise auxiliary variables
        indOpt = zeros(nkap, N*nkap); % Indicator: =1 if optimal decision
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

% Other returns
        % Aggregate capital
        aggregateK = sum(kk(:).*dist(:));

display('done with measure')