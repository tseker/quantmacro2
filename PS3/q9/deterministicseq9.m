function [k_imp,prodN_imp, lambda_imp,v_t,lambda_dist]   = deterministicseq9(r_t,lambda_t,T, v, v_post, dist,dist_post, kk,optV);

global beta mu delta alpha s N prob kap nkap A phi tau2 B G
util = @(c) (c.^(1-mu)) / (1-mu);                                     
w_t     = (1-alpha)*(A.*(alpha./(r_t+delta)).^alpha).^(1/(1-alpha));  % wage sequence


%% HOUSEHOLD'S PROBLEM

% We will operate with 2 dimensions: rows for (s,a) and columns for (t).

% Matrix of states
states(:,1) = reshape(repmat([1:N],nkap,1),N*nkap,1);
states(:,2) = repmat([1:nkap]',N,1);
s_vec = s(states(:,1))'; 
kap_vec = kap(states(:,2))';

% Last period: convergence to steady state
aprime_t(:,T) = reshape(kk',N*nkap,1);   % kk needs to be ajusted  

% Solve labor at period T.
% lower bound: imposed so that c>0
nlow = ( max(  (aprime_t(:,T)-(1+r_t(T)).*kap_vec)./lambda_t(T),0  )).^(1/(1-tau2))./(w_t(T).*s_vec) ;
        nlow = max(nlow,0);
        nup = ones(N*nkap,1);
        nlow(nlow>nup) = 1-1e-5;
        % Loop criterium
        maxtestn = 1;
        n0 = 0.3*ones(N*nkap,1);
        prodwage = w_t(T).* s_vec;
    % Loop: we combine newton method and bisection method.
    % bisection method is robust in the sense that it will not create a negative consumption
    % monotonicity of the focn function ensures that bisection works.
    % Newton is faster but not robust.
    while (maxtestn>0.0001)
        foc = focn(n0,kap_vec,aprime_t(:,T), lambda_t(T), prodwage, r_t(T));
        first_dev = first_dev_focn(n0,kap_vec,aprime_t(:,T), lambda_t(T), prodwage, r_t(T));
        n1 = n0 - foc./first_dev; 
        % if n1 out of the bound, 
        n1(n1<nlow | n1>nup ) = 0.5*nlow(n1<nlow | n1>nup) + 0.5*nup(n1<nlow|n1>nup);
        foc1 = focn(n1,kap_vec,aprime_t(:,T), lambda_t(T), prodwage, r_t(T));
        nlow(foc1<0 & n1<=nup & n1>=nlow) = n1(foc1<0 & n1<=nup & n1>=nlow);
        nup(foc1>=0 & n1<=nup & n1>=nlow) = n1(foc1>=0 & n1<=nup & n1>=nlow);
        maxtestn = max(abs(foc1));
        n0 = n1;
    end
l_t = zeros(N*nkap,T);
l_t(:,T) = n0;
c_t = zeros(N*nkap,T);
c_t(:,T) = (1+r_t(T))*kap_vec + lambda_t(T).*( w_t(T).*s_vec.*l_t(:,T)).^(1-tau2) - aprime_t(:,T); %c_t needs to be adjusted





% Computing backwards until initial period
for t = T-[1:(T-1)]
    
    % Implied consumption and savings
    auxC = reshape(c_t(:,t+1),nkap,N)';
    c_imp = ( beta*(1+r_t(t)) * prob * (auxC.^(-mu)) ).^(-1/mu);  % Euler Equation doesn't change
    c_imp = reshape(c_imp',N*nkap,1);
    l_imp = (1/B.*lambda_t(t).*(1-tau2) .*(w_t(t).*s_vec).^(1-tau2) .* c_imp.^(-mu) ).^(1/(phi+tau2));
    a_imp = ( kap_vec + c_imp - lambda_t(t).* (w_t(t).*s_vec .*l_imp).^(1-tau2))./(1+r_t(t));             % Budget Constraint change
    
    
    % Policy function: savings (by interpolation)
    a_imp = reshape(a_imp,nkap,N)';
    c_imp = reshape(c_imp,nkap,N)';
    for j=1:N
        auxAprime(j,:) = interp1(a_imp(j,:),kap,kap,'linear','extrap');
    end
    % Bounds
    auxAprime(auxAprime<kap(1))=kap(1);
    auxAprime(auxAprime>kap(end))=kap(end);
    

    % Policy functions 
        % Savings
        aprime_t(:,t) = reshape(auxAprime',N*nkap,1);
        clear auxAprime
        % Consumption & labor
nlow = ( max(  (aprime_t(:,t)-(1+r_t(t)).*kap_vec)./lambda_t(t),0  )).^(1/(1-tau2))./(w_t(t).*s_vec) ;
        nlow = max(nlow,0);
        nup = ones(N*nkap,1);
        nlow(nlow>nup) = 1-1e-5;
        % Loop criterium
        maxtestn = 1;
        n0 = 0.3*ones(N*nkap,1);
        prodwage = w_t(t).* s_vec;
    % Loop: we combine newton method and bisection method.
    % bisection method is robust in the sense that it will not create a negative consumption
    % monotonicity of the focn function ensures that bisection works.
    % Newton is faster but not robust.
    while (maxtestn>0.0001)
        foc = focn(n0,kap_vec,aprime_t(:,t), lambda_t(t), prodwage, r_t(t));
        first_dev = first_dev_focn(n0,kap_vec,aprime_t(:,t), lambda_t(t), prodwage, r_t(t));
        n1 = n0 - foc./first_dev; 
        % if n1 out of the bound, 
        n1(n1<nlow | n1>nup ) = 0.5*nlow(n1<nlow | n1>nup) + 0.5*nup(n1<nlow|n1>nup);
        foc1 = focn(n1,kap_vec,aprime_t(:,t), lambda_t(t), prodwage, r_t(t));
        nlow(foc1<0 & n1<=nup & n1>=nlow) = n1(foc1<0 & n1<=nup & n1>=nlow);
        nup(foc1>=0 & n1<=nup & n1>=nlow) = n1(foc1>=0 & n1<=nup & n1>=nlow);
        maxtestn = max(abs(foc1));
        n0 = n1;
    end
l_t(:,t) = n0;
c_t(:,t) = (1+r_t(t))*kap_vec + lambda_t(t).*( w_t(t).*s_vec.*l_t(:,t)).^(1-tau2) - aprime_t(:,t); %c_t needs to be adjusted

end



%% MEASURE (as in the Aiyagari function)

lambda_dist = zeros(N*nkap,T);
lambda_dist(:,1) = reshape(dist',N*nkap,1);                   % Steady state lambda; sets the initial distribution at time period 1
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
    lambda_dist(:,t+1) = Q*lambda_dist(:,t);
    % Reshape distribution into matrix form
    %dist_t(:,:,t) = transpose(reshape(lambda(:,t),nkap,N));
end

disp('Done with backward and forward iteration')


%% VALUE FUNCTION
% code iteratively calculates the value function backward in time by combining the utility from consumption, 
% the disutility from labor, and the discounted expected continuation value. 
% It considers the transition probabilities and the interpolation of continuation values between nearest neighbors in the capital grid.

if optV
    % Probability of transition to s'
    prob_vec = reshape(prob,N^2,1); % make prob column vector
    prob_vec = repmat(prob_vec,1,nkap); % replicate nkap columns
    prob_vec = reshape(prob_vec',nkap*N,N); % reshape: s' in columns
    % reshape technique: the order of elements is according to columns, the
    % first colum comes first
    % Convergence to steady state value
    v_t(:,T) = reshape(v_post',N*nkap,1); % v needs to be adjusted
    % Get all the sequence
    for t=T-[1:T-1]
        expV = prob_vec * reshape(v_t(:,t+1),nkap,N)';
        % Interpolate continuation value
        cont = (1-indW(:,t)).*expV(sub2ind(size(expV), [1:N*nkap]', indU(:,t))) +...
               indW(:,t).*expV(sub2ind(size(expV), [1:N*nkap]', indL(:,t)));
        v_t(:,t) = util(c_t(:,t)) - B*l_t(:,t).^(1+phi)./(1+phi) + beta*cont;
    end
end

%% RETURNS

% Implied capital
% something wrong
k_imp = sum(repmat(kap',5,T).*lambda_dist);
prodN_imp = sum((l_t.*repmat(s_vec,1,T)).*lambda_dist);

lambda_imp = (w_t.*prodN_imp-G)./sum(  (w_t.*l_t.*repmat(s_vec,1,T)).^(1-tau2).*lambda_dist);



% Value
if ~optV
    v_t=[];
end