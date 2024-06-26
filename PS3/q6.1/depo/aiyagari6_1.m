function [aggK,aggN] = aiyagari6_1(r0,optHoward)
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
focn = @(sj,kap,kapp,n) B*n.^(phi)-(wage*sj*((1+r0)*kap' + wage*sj*n - kapp).^(-mu));
fprima = @(sj,kap,kapp,n) phi*B*n.^(phi-1)+(mu*wage*sj*((1+r0)*kap' + wage*sj*n - kapp).^(-mu-1)*wage*sj);
util = @(cons,n) (cons.^(1-mu))/(1-mu)-B*(n.^(1+phi))/(1+phi);

%% %% Endogenous labour supply
n0=0.9*ones(N,nkap,nkap);
nnew=zeros(N,nkap,nkap);
lambda = 0.1;
f0=zeros(N,nkap,nkap);
fprime=zeros(N,nkap,nkap);
maxtestn = 1;
while (maxtestn>1e-3)
for j=1:N               
    for i=1:nkap        
        for k=1:nkap    
            % Loop
            % while (maxtestn>1e-3)
                % Update guess
                f0(j,i,k) = B*n0(j,i,k)^(phi)-(wage*s(j)*((1+r0)*kap(i) + wage*s(j)*n0(j,i,k) - kap(k))^(-mu));
                fprime(j,i,k) = phi*B*n0(j,i,k).^(phi-1)+(mu*wage*s(j)*((1+r0)*kap(i) + wage*s(j)*n0(j,i,k) - kap(k)).^(-mu-1)*wage*s(j));
                %nnew = max(min(n0-((focn(s(j),kap(i),kap(k),n0))/fprima(s(j),kap(i),kap(k),n0)),1),0);
                % n0 = nnew;
                % % Indicator of convergence
                % maxtestn=focn(s(j),kap(i),kap(k),n0);
                % display(maxtestn)
        end
         nnew = max(min(n0-f0./fprime,1),0);
         maxtestn = max(max(abs(nnew-n0)));
         display(maxtestn)
         n0 = nnew;
        % nopt(i,k,j)=n0;
     end
end
end


%% Value Function Iteration

%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum 

% Seed for the value function:
if isempty(v) % when this is true, we are in the first iteration
   v = ((wage*repmat(s',1,nkap) + r0*repmat(kap,N,1)).^(1-mu))/((1-mu)*(1-beta));
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
end

%Solve the value function
while test > 0.001
    iter = iter + 1;
    % Optimisation
    for j = 1:N
        for i = 1:nkap
        auxf = @(x) newV(x,i, j, wage, r0);
        upperK = min(maxkap, wage*s(j) + (1+r0)*kap(i) - 0.001);
        [aprime_VFI(j,i), tv(j,i)] = ...
                        golden6_1(auxf,minkap,upperK);
            % We extract the nearest elements on the grid and the weights,
            [~,tdecis(j,i)]= min(abs(aprime_VFI(j,i) - kap));
        end
    end
    
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