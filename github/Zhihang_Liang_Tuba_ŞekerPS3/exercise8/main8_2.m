%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kk kap v dist phi B n_pol c_pol tau

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
    tau    = 0.3;            % Tax progressivity
    G      = 0.06;            % Government spending

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

% Loop on r
    % Initial bounds
        minrate = 0.02;
        maxrate = 0.03;
    % Initial guess for interest rate
        r0 = 0.5*(maxrate-minrate)+minrate;
    % Initialise distance test
        testr=1;
while (testr>0.0001)
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
        [K1,prodN1]=aiyagari8_1(r0,w0,lambda0); % calculating asset and labour supply
    % Quantitites
       % K0 = ((r0+delta)/(alpha*A*prodN1^(1-alpha)))^(1/(alpha-1)); % capital demand for that guess
       % Y0 = A * K0^alpha * prodN0^(1-alpha); % production
       % G = Gshare*Y0; % government spending
    % Implied lambda
        den = sum(sum( (w0 * s' .* n_pol).^(1-tau) .* dist )); % denominator
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


%% FIGURES

nkap=length(kap);
kopt=reshape(kk,N,nkap);

figure;

subplot(3,2,1);
plot(kap,[kap' kopt(1,:)' kopt(3,:)' kopt(5,:)'])
ylabel('next period asset holdings')
xlabel('current asset holdings')
hleg1 = legend('Constant Asset Line','Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northwest');

subplot(3,2,2);
plot(kap,[dist(1,:)' dist(5,:)'])
ylabel('distribution')
xlabel('asset holdings')
hleg1 = legend('Asset Distribution - Lowest Skill Shock', 'Asset Distribution - Highest Skill Shock');
 
subplot(3,2,3);
plot(s',[dist(:,1) dist(:,nkap)])
ylabel('distribution')
xlabel('Skill Shock')
hleg1 = legend('Asset Distribution - Lowest Asset Level', 'Asset Distribution - Highest Asset Level');

subplot(3,2,4);
plot(kap,[c_pol(1,:)' c_pol(3,:)' c_pol(5,:)'])
ylabel('consumption')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northwest');

subplot(3,2,5);
plot(kap,[n_pol(1,:)' n_pol(3,:)' n_pol(5,:)'])
ylabel('labour supply')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northeast');

subplot(3,2,6);
plot(kap,[v(1,:)' v(3,:)' v(5,:)'])
ylabel('value function')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','southeast');

toc