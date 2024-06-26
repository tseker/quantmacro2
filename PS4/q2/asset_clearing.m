function [asset_mkt] = asset_clearing( real_i, labor_t )

% bewley_ss.m compute the steady state distribution, policyfunctions and
% the asset market clearing condition, taken a guess for interest
% rate as given.

% a guess for interest rate imply a gueess for the tax rate
% rB = tau*w*L

tau_t = par.B./(par.A.*labor_t*par.tech_labor);
c_t = zeros(par.nkap,par.nkap,par.T);
kk_t =  zeros(par.nkap,par.nkap,par.T);
kk = zeros(par.nkap,par.nkap);
wage = par.A;


%% EGM to perform the value function iteration

% Endogenous Grid Method

for t = T-(1:(T-1))
    c_imp = (par.beta*(1+r_guess)*par.prob*c_t(:,:,t).^(-par.gamma)).^(-1/par.gamma);
    a_imp = (par.kap+c_imp-(1-tau_t(:,t))*wage*par.L*par.s')/(1+r_guess);
    for j=1:par.N
        kk(j,:) = interp1(a_imp(j,:),par.kap,par.kap, 'linear','extrap');
    end
    kk(kk<par.minkap) = par.minkap;
    kk(kk>par.maxkap) = par.maxkap;
    kk_t(:,:,t) = kk;
    c_t(:,:,t) = wage*(1-tau_t(t))*par.L*par.s'+(1+real_i(t))*par.kap-kk;
end



disp('done with EGM, transition')


indU = zeros(par.N,par.nkap);
    for j = 1:par.N
        for i = 1:par.nkap
            if isempty(find(par.kap>=policy.kk(j,i),1,'first'))
                disp('wrong i,j');
                [i,j]
                indU(j,i) = par.nkap;
            else
                indU(j,i)= find(par.kap>=policy.kk(j,i),1,'first');
            end
        end
    end





%% Measure

% Initialise auxiliary variables
    indOpt = sparse(par.nkap, par.N*par.nkap); % Indicator: =1 if optimal decision
    Q = []; % Transition matrix
    auxAprime = reshape(policy.kk',1,par.N*par.nkap); % now each column represents a combination of s and a
% Create matrix with 1s in the position of optimal decision
    % - row indicates the index of optimal a' on the grid.
    % - column indicates the starting par (s,a) on the grid.
    % Due to interpolation, some a' are not on the grid. Use weights
    % for the nearest elements that belong to the grid.
        % Find nearest elements and their weights
            % Upper neighbour: already found in golden search (indU)
            % Reshape: to convert in row vector
            auxU = reshape(indU',1,par.N*par.nkap);
            % Lower neighbour:
            auxL = auxU-1;
            % When a' is the first element in kap, the lower neighbour
            % will not matter, all the weight will be on upper. But we
            % need auxL to remain on the grid, so it cannot be 0. It
            % cannot be 1 either, or we will have auxL=auxU and the
            % weight will not be well defined.
            auxL(auxL==0) = 2;
            % Weight:
            auxW = (par.kap(auxU) - auxAprime) ./ (par.kap(auxU) - par.kap(auxL));
        % Assign weight to the lower neighbour
        indOpt(sub2ind(size(indOpt), auxL, 1:(par.N*par.nkap))) = auxW;
        % Assing weight to the upper neighbour
        indOpt(sub2ind(size(indOpt), auxU, 1:(par.N*par.nkap))) = 1-auxW;

% Multiplying by the probability of each state, we obtain the
% transition matrix (Q)
    for j=1:par.N % loop over each productivity state
        % We need to have a vector of probabilities with N*nkap
        % elements (it will be each of the pi(j | s0) repeated nkap
        % times).
            auxProb = reshape(transpose(repmat(par.prob(:,j),1,par.nkap)),1,par.N*par.nkap);
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
        mu0 = ones(par.N*par.nkap,1) / (par.N*par.nkap); % starting distribution
    % Get stationary distribution
        while distance>=toler
            mu1 = Q*mu0;
            distance=max(max(abs(mu1-mu0)));
            mu0 = mu1;
        end
    % Reshape distribution:
        % - each row represents a state s'.
        % - each column represents the decision on assets a'.
        dist = transpose(reshape(mu1,par.nkap,par.N));

display('done with measure')


%% Other returns
% Aggregate capital
    k1 = sum(policy.kk(:) .* dist(:));
    
% asset market
    asset_mkt = k1-par.B;
    

end