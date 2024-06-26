%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kap  phi B c_pol tau1 tau2 maxkap nkap minkap G

tic

%% PRELIMINARIES

%  Parameter values
    mu     = 2;               % risk aversion. CRRA              
    beta   = 0.95;            % subjective discount factor 
    delta  = 0.1;             % depreciation
    A      = 1;               % production technology
    alpha  = 0.36;            % capital's share of income
    b      = 0.2;             % Borrowing constraint
    B      = 150;             % disutility of labour
    phi    = 2.5;             % Frisch elasticity
    tau1   = 0.18;            % Tax progressivity, initial
    tau2   = 0.3;             % Tax progressivity, post reform
    G      = 0.06;            % Government spending
    
    % Capital Grid
    maxkap = 50;                     % maximum value of capital grid  
    minkap = -b;                     % borrowing constraint                    
    nkap = 150;                      % number of asset grids
    kap = logspace(0,log10(maxkap+1-minkap),nkap)-1+minkap;

% approximate labor endowment shocks with 5 states Markov chain
% log(s_t) = rho*log(s_t-1)+sigma*error_t
    N        = 5;             % number of discretized states
    rho      = 0.9;           % first-order autoregressive coefficient
    sigma    = 0.2;           % standard deviation of error_t
    
    % prob is transition matrix of the Markov chain
    % logs is the discretized states of log labor earnings
    % invdist is the invariant distribution of Markov chain

    % s productivity
    
    [logs,prob] = tauchen(N,0,rho,sigma,3); logs = logs';
        % prob contains probability of transition from s0 (row) to s1 (column)
    s = exp(logs);

% Compute invariant distribution
    invdist = ones(1,N)/N; test = 1; 
    
    while (test>0.0000001)
        invdist2 = invdist*prob;
        test = max(abs(invdist2-invdist));
        invdist = invdist2;
    end
    invdist = invdist';
    prodN0 = s*invdist;

%% GENERAL EQUILIBRIUM (steady state)

% Steady State: before reform
% Loop on r
    % Initial bounds
        minrate = 0.02;
        maxrate = 0.03;
    % Initial guess for interest rate
        r0 = 0.5*(maxrate-minrate)+minrate;
    % Initialise distance test
        testr=1;
        w0=(1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));
        v = 1/(1-beta) * ((s'*w0 + r0*kap).^(1-mu))/(1-mu);
while (testr>0.001)
% Wages
    w0=(1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha)); % wage from Cobb-Douglas production function
% Loop on lambda
    % Initial guess
    if ~exist("lambda0",'var')
        lambda0 = 0.7;
    end
    % Initialise distance test
    testlambda=1;
    while (testlambda>0.001)
    % Update lambda
%         lambda0 = 1/2*maxlambda + 1/2*minlambda;
    % Aiyagari
        % Technique: use previous guess for v
        [K1,prodN1,dist,v,n_pol,kk]=aiyagari9(r0,w0,lambda0, tau1, v); % calculating asset and labour supply
    % Quantitites
       % K0 = ((r0+delta)/(alpha*A*prodN1^(1-alpha)))^(1/(alpha-1)); % capital demand for that guess
       % Y0 = A * K0^alpha * prodN0^(1-alpha); % production
       % G = Gshare*Y0; % government spending
    % Implied lambda
        den = sum(sum( (w0 * s' .* n_pol).^(1-tau1) .* dist )); % denominator
        lambda1 = (w0*prodN1 - G)/den;
        disp('lambda guess, lambda implied  = ')
        [lambda0 lambda1]
    % Update test
        testlambda=abs(lambda1-lambda0);
%     % Update bounds
%         if lambda1 > lambda0
%             minlambda = lambda0;
%         else
%             maxlambda = lambda0;
%         end
    % Update lambda
        lambda0 = 0.8*lambda0 + 0.2*lambda1;
    end
% Implied interest rate
    r1 = alpha*A * max(0.001,K1)^(alpha-1) * prodN1^(1-alpha) - delta;
    disp('r guess, r implied  = ')
    [r0 r1]
% Update test on r
    testr=abs(r1-r0); 
% Updating the interest rate
    r0 = 0.9*r0 + 0.1*r1; % guess for interest rate
end

% Relevant quantities
    % Production
    Y1 = A * K1^alpha * prodN1^(1-alpha);
    % Gov spending share
    Gshare = G/Y1;
    % Hours worked
    meanH = sum(n_pol(:) .* dist(:));
    % Aggregate welfare
    aggW = sum(v(:) .* dist(:));
    % Show results
    fprintf("Lambda = %4.2f.\n" + ...
            "The ratio of government spending to GDP is: %6.4f.\n" + ...
            "Average hours worked are: %6.4f.\n" + ...
            "Aggregate welfare is: %6.4f.\n" + ...  
            "Aggregate capital is: %6.4f.\n", ...
            lambda0,Gshare,meanH,aggW,K1);
        
        
 
%% Steady State: after reform
% Loop on r
    % Initial bounds
        minrate = 0.02;
        maxrate = 0.03;
    % Initial guess for interest rate
        r0_post = 0.5*(maxrate-minrate)+minrate;
    % Initialise distance test
        testr=1;
        v_post = 1/(1-beta) * ((s'*w0 + r0*kap).^(1-mu))/(1-mu);
while (testr>0.001)
% Wages
    w0_post=(1-alpha)*(A*(alpha/(r0_post+delta))^alpha)^(1/(1-alpha)); % wage from Cobb-Douglas production function
% Loop on lambda
    % Initial guess
    if ~exist("lambda0_post",'var')
        lambda0_post = 0.7;
    end
    % Initialise distance test
    testlambda=1;
    while (testlambda>0.001)
    % Update lambda
%         lambda0 = 1/2*maxlambda + 1/2*minlambda;
    % Aiyagari
        [K1_post,prodN1_post,dist_post,v_post,n_pol_post, kk_post]=aiyagari9(r0_post,w0_post,lambda0_post, tau2, v_post); % calculating asset and labour supply
    % Quantitites
       % K0 = ((r0+delta)/(alpha*A*prodN1^(1-alpha)))^(1/(alpha-1)); % capital demand for that guess
       % Y0 = A * K0^alpha * prodN0^(1-alpha); % production
       % G = Gshare*Y0; % government spending
    % Implied lambda
        den_post = sum(sum( (w0_post * s' .* n_pol_post).^(1-tau2) .* dist_post )); % denominator
        lambda1_post = (w0_post*prodN1_post - G)/den_post;
        disp('lambda guess, lambda implied  = ')
        [lambda0_post lambda1_post]
    % Update test
        testlambda=abs(lambda1_post-lambda0_post);
%     % Update bounds
%         if lambda1 > lambda0
%             minlambda = lambda0;
%         else
%             maxlambda = lambda0;
%         end
    % Update lambda
        lambda0_post = 0.8*lambda0_post + 0.2*lambda1_post;
    end
% Implied interest rate
    r1_post = alpha*A * max(0.001,K1_post)^(alpha-1) * prodN1_post^(1-alpha) - delta;
    disp('r guess, r implied  = ')
    [r0_post r1_post]
% Update test on r
    testr=abs(r1_post-r0_post); 
% Updating the interest rate
    r0_post = 0.9*r0_post + 0.1*r1_post; % guess for interest rate
end

% Relevant quantities
    % Production
    Y1_post = A * K1_post^alpha * prodN1_post^(1-alpha);
    % Gov spending share
    Gshare_post = G/Y1_post;
    % Hours worked
    meanH_post = sum(n_pol_post(:) .* dist_post(:));
    % Aggregate welfare
    aggW_post = sum(v_post(:) .* dist_post(:));
    % Show results
    fprintf("Lambda = %4.2f.\n" + ...
            "The ratio of government spending to GDP is: %6.4f.\n" + ...
            "Average hours worked are: %6.4f.\n" + ...
            "Aggregate welfare is: %6.4f.\n" + ...  
            "Aggregate capital is: %6.4f.\n", ...
            lambda0_post,Gshare_post,meanH_post,aggW_post,K1_post);       

%% Tax reform: progressivity

% Parameters
T       = 100;                                                          % number of periods

% Initialise variables
testr           = 1;
testlambda      = 1;
tolr            = 1e-5;
tollambda       = 1e-4;
r_t             = r0_post*ones(1,T) ;      % initial guess for r: sudden jump to new SS
lambda_t        =  lambda0_post.*ones(1,T); % initial guess for lambda


% Finding the sequence
while testr > tolr & testlambda> tollambda
    % Implied variables
    [k_imp,prodN_imp, lambda_imp]  = deterministicseq9(r_t,lambda_t,T, v, v_post, dist, dist_post,kk_post, false);
    r_imp   = alpha * A .* max(0.001,k_imp).^(alpha-1) .* prodN_imp.^(1-alpha) - delta;
    % Update test and guess
    testlambda = max(abs(lambda_t-lambda_imp));
    lambda_t = 0.7* lambda_t + 0.3 * lambda_imp;
    testr   = max(abs(r_t-r_imp));
    r_t = 0.95 * r_t + 0.05 * r_imp;
end

% Value right after the shock
[~,~,~,v_t] = deterministicseq9(r_t,lambda_t,T, v, v_post, dist, dist_post,kk_post, true);
v0 = reshape(v_t(:,1),nkap,N)';

% Consumption equivalence gain
c_eq = (((1-mu)*v0+1).^(1/(1-mu)) -  ((1-mu)*v+1).^(1/(1-mu)) )./((1-mu)*v+1).^(1/(1-mu));
c_eq_mean = sum(c_eq(:).*dist(:));
fprintf("The average consumption equivalent is %3.2f%%.\n",100*c_eq_mean);

% Wealth gain



toc



%% GRAPHS (deterministic path)

w_t     = (1-alpha)*(A.*(alpha./(r_t+delta)).^alpha).^(1/(1-alpha));  % wage sequence

fig=figure(2); % Transitions

% General options
sim=1; % take the first simulation
fsize=12;
lwidth=1.5;

% Capital
subplot(3,1,1)
    hold on
    plot(k_imp, 'LineWidth', lwidth) % IRF
    yline(K1_post,'--','LineWidth', 2/3*lwidth) % SS value
    % Options
    title("Path of capital")
    grid on
    xlabel('Time')
    ylabel('Stock of capital')
    set(gca,'FontSize',fsize)
    hold off
% Interest rate
subplot(3,1,2)
    hold on
    plot(r_t, 'LineWidth', lwidth) % IRF
    yline(r0_post,'--','LineWidth', 2/3*lwidth); % SS value
    % Options
    title("Path of interest rate")
    grid on
    xlabel('Time')
    ylabel('Interest rate')
    set(gca,'FontSize',fsize)
    hold off
% Wage
subplot(3,1,3)
    hold on
    plot(w_t, 'LineWidth', lwidth) % IRF
    yline(w0_post,'--','LineWidth', 2/3*lwidth); % SS value
    % Options
    title("Path of wages")
    grid on
    xlabel('Time')
    ylabel('Wage')
    set(gca,'FontSize',fsize)
    hold off

saveas(fig,'main_9.png');


fig2= figure(3); % Consumption equivalents
    plot(kap,[c_eq(1,:)' c_eq(3,:)' c_eq(5,:)'],'LineWidth',lwidth)
    % Options
    grid on
    ylabel('consumption equivalent')
    xlabel('asset holdings')
    set(gca,'FontSize',fsize)
    hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northeast');

saveas(fig2,'ceq_9.png');