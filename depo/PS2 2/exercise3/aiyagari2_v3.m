
function aggregateK = aiyagari2_v3(r0)
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
   
maxkap = 50;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.015;                    % size of capital grid increments
kap    = log(minkap+b+0.1):inckap:log(maxkap+b+0.1);   % state of assets 
kap = exp(kap)-b-0.1;
nkap   = length(kap);            % number of grid points

%  initialize some variables
criterion = 10;
iter = 0;
   

% Vectorize to create matrix of states for each (a, epsilon):
    % Create indexes for each:
        ind.z= 1;       % Productivity
        ind.a = 2;      % Savings
        ind.a0 = find (kap>=0);  % Find positions for zero asset

     % Create a matrix
      mat.ea(:, ind.z) = kron((1:N)', ones(nkap,1));
      mat.ea(:, ind.a) = kron(ones(N,1),(1:nkap)');


     
%%

% Iterate on Bellman's equation and get the decision 
% %  rules and the value function at the optimum 
% 
% for j=1:N
%     for i=1:nkap
%         cons(j,i) = s(j)*wage + r0*kap(i);
%         utilm(j,i) = (cons(j,i).^(1-mu)-1)/(1-mu);
%         v(j,i) = utilm(j,i)/(1-beta);
%     end
% end
% 
% utilm = 190*zeros(N,nkap,nkap);
% 
% for j=1:N%loop over each skill s_t
%     for i=1:nkap% loop over each possible a_t            
%         cons = s(j)*wage + (1+r0)*kap(i) -kap';
%         util = (cons.^(1-mu)-1)/(1-mu);
%         ii = find( cons<= 0);              % infeasible consumption choice
%         util(ii) = -10000;            
%         utilm(j,i,:)=util;
%     end
% end
% 
% % Initialize matrices and distance
% tv = 189*zeros(N,nkap); 
% tdecis = tv; 
% test    = 10;
% 
% 
% tdecis = zeros(N,nkap);
% aprime_VFI = zeros(N,nkap);
% %Solve the value function
% while test > 0.001
%     iter = iter + 1;
%     % Optimisation
%     for j = 1:N
%         for i = 1:nkap
%         auxf = @(x) newV(x,i, j, wage, r0);
%         upperK = min(max(kap), wage*s(j) + (1+r0)*kap(i) - 0.0001);
%         [aprime_VFI(j,i), tv(j,i)] = ...
%                         golden(auxf,minkap,upperK);
%             % We extract the nearest elements on the grid and the weights,
%             [~,tdecis(j,i)]= min(abs(aprime_VFI(j,i) - kap));
%         end
%     end
% 
%     c_VFI = s'*wage + (1+r0)*kap - aprime_VFI;
% 
% 
% 
%     test=max(max(abs(tv-v))); % distance
%     v=tv;
% end
  
 
% display('done with VF iteration')

%%
% Endogenous Grid Method

% Initial guess for the policy function
c0 = wage*s'+r0*kap; 

% Initialising matrices and distance
test = 10;
iter = 0;

while test>1e-8
    iter= iter+1;
    %Implied consumption:
    c_imp = (beta*(1+r0)*prob*(c0).^(-mu)).^(1/-mu);

    %Implied savings:
    % Given c_imp, a', and productivity, obtain a by budget constraint
    a_imp = (kap+c_imp-wage*s')/(1+r0);

    
    % Make sure a (savings) on the grid
    for j=1:N
        kk(j,:) = interp1(a_imp(j,:),kap,kap, 'linear','extrap');
    end

    % Borrowing constraint bounds to prevent error in the measure while
    % calculating the measure
    kk(kk<minkap) = minkap;
    kk(kk>maxkap) = maxkap;

    % Consumption policy function
    c_egm = wage*s'+(1+r0)*kap-kk;

    % Distance measure
    test = max(max(abs(c_egm-c0)));
    c0 = c_egm;
    
end

display('done with EGM')

tdecis = zeros(N,nkap);
    for j = 1:N
        for i = 1:nkap
            [~,tdecis(j,i)]= min(abs(kk(j,i) - kap));
        end
    end
   






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