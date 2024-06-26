function [K_imp, v_t] = deterministicseq7(r_t, A_t, T, optV)

global beta mu delta alpha s N prob kk kap nkap v dist
util = @(c) (c.^(1-mu)) / (1-mu);                                     
w_t     = (1-alpha)*(A_t.*(alpha./(r_t+delta)).^alpha).^(1/(1-alpha));  % wage sequence


%% HOUSEHOLD'S PROBLEM

% We will operate with 2 dimensions: rows for (s,a) and columns for (t).

% Matrix of states
states(:,1) = reshape(repmat([1:N],nkap,1),N*nkap,1);
states(:,2) = repmat([1:nkap]',N,1);
s_vec = s(states(:,1))';
kap_vec = kap(states(:,2))';

% Last period: convergence to steady state
aprime_t(:,T) = reshape(kk',N*nkap,1);                                    % convergence to steady state savings
c_t(:,T) = (1+r_t(T))*kap_vec + w_t(T)*s_vec - aprime_t(:,T);

% Computing backwards until initial period
for t = T-[1:(T-1)]
    
    % Implied consumption and savings
    auxC = reshape(c_t(:,t+1),nkap,N)';
    c_imp = ( beta*(1+r_t(t)) * prob * auxC.^(-mu) ).^(1/-mu);  % Euler Equation
    a_imp = ( kap + c_imp - w_t(t)*s' )/(1+r_t(t));             % Budget Constraint

    % Policy function: savings (by interpolation)
    for j=1:N
        auxAprime(j,:) = linInterp(kap,a_imp(j,:),kap);
    end
    % Bounds
    auxAprime(auxAprime<kap(1))=kap(1);
    auxAprime(auxAprime>kap(end))=kap(end);

    % Policy functions 
        % Savings
        aprime_t(:,t) = reshape(auxAprime',N*nkap,1);
        clear auxAprime
        % Consumption
        c_t(:,t) = w_t(t)*s_vec + (1+r_t(t))*kap_vec - aprime_t(:,t);
end


%% MEASURE (as in the Aiyagari function)

lambda = reshape(dist',N*nkap,1);                   % Steady state lambda
for t = 1:T-1
    Q = [];                                         % Transition Matrix
    indOpt = sparse(nkap, N*nkap);                  % Indicator: =1 if optimal decision
    % Find nearest neighbours
    for k=1:(N*nkap)
        indU(k,t) = find(kap >= aprime_t(k,t),1,'first');
    end
    auxL = indU(:,t)-1;
    auxL(auxL==0) = 2;
    indL(:,t) = auxL;
    indW(:,t) = (kap(indU(:,t))' - aprime_t(:,t)) ./ (kap(indU(:,t))' - kap(auxL)');
    % Create matrix with 1s in the position of optimal decision
    indOpt(sub2ind(size(indOpt), auxL, [1:(N*nkap)]')) = indW(:,t);
    indOpt(sub2ind(size(indOpt), indU(:,t), [1:(N*nkap)]')) = 1-indW(:,t);
    % Multiplying by the probability of each state, we obtain the transition matrix (Q)
    for j=1:N
        auxProb = reshape(transpose(repmat(prob(:,j),1,nkap)),1,N*nkap);
        auxQ = indOpt .* auxProb;
        Q = [Q; auxQ];
    end
    lambda(:,t+1) = Q*lambda(:,t);
    % Reshape distribution into matrix form
    %dist_t(:,:,t) = transpose(reshape(lambda(:,t),nkap,N));
end


%% VALUE FUNCTION

if optV
    % Probability of transition to s'
    prob_vec = reshape(prob,N^2,1); % make prob column vector
    prob_vec = repmat(prob_vec,1,nkap); % replicate nkap columns
    prob_vec = reshape(prob_vec',nkap*N,N); % reshape: s' in columns
    % reshape technique: the order of elements is according to columns, the
    % first colum comes first
    % Convergence to steady state valu  e
    v_t(:,T) = reshape(v',N*nkap,1);
    % Get all the sequence
    for t=T-[1:T-1]
        expV = prob_vec * reshape(v_t(:,t+1),nkap,N)';
        % Interpolate continuation value
        cont = (1-indW(:,t)).*expV(sub2ind(size(expV), [1:N*nkap]', indU(:,t))) +...
               indW(:,t).*expV(sub2ind(size(expV), [1:N*nkap]', indL(:,t)));
        v_t(:,t) = util(c_t(:,t)) + beta*cont;
    end
end


%% RETURNS

% Implied capital
K_imp = sum(repmat(kap',5,1) .* lambda);

% Value
if ~optV
    v_t=[];
end