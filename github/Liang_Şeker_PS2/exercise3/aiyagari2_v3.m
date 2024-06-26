
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
maxkap = max(kap);
minkap = min(kap);

%  initialize some variables
criterion = 10;
iter = 0;

     
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

disp('done with EGM')

indU = zeros(N,nkap);
    for j = 1:N
        for i = 1:nkap
            if isempty(find(kap>=kk(j,i),1,'first'))
                disp('wrong i,j');
                [i,j]
            end
            indU(j,i)= find(kap>=kk(j,i),1,'first');
        end
    end
   






% Compute the measure
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

aggregateK = sum(kk(:).*dist(:));