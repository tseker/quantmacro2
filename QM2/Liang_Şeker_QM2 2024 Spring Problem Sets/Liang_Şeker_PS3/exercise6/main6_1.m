%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kk kap v dist phi B n_pol c_pol Ntarget

tic

%  set parameter values

mu     = 2;               % risk aversion. CRRA              
beta   = 0.95;            % subjective discount factor 
delta  = 0.1;             % depreciation
A      = 1;               % production technology
alpha  = 0.36;            % capital's share of income
B      = 150;              % disutility of labour
phi    = 2.5;             % Frisch elasticity
Ntarget = 0.3;

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

% Compute invariant distribution
invdist = ones(1,N)/N; test = 1; 

while (test>0.0000001)
    invdist2 = invdist*prob;
    test = max(abs(invdist2-invdist));
    invdist = invdist2;
end
invdist = invdist';

s = exp(logs);
labor = s*invdist;


b = 0.2;  % BORROWING CONSTRAINT

% setting limits for interest rate search

minrate = 0.02;
maxrate = 0.03;

% setting test values
testr=1;
testn=1;


% we find the optimal interest rate by bisection
r0 = 0.015;


while (testr>0.0001 )
    k0=((r0+delta)/(alpha*A*Ntarget^(1-alpha)))^(1/(alpha-1));%capital demand for that guess
    [k1,n1]= aiyagari6_1(r0); %calculating asset supply
    r1= alpha*A*max(0.001,k1)^(alpha-1)*n1^(1-alpha)-delta;%which interest rate implies that 

    % Implied interest rate:
    r0 = 0.9*r0 + 0.1*r1;  %Guess the interest rate

    %capital demand equals this value of capital supply 
    display('k0, k1, r guess, r implied, target n, computed n  = ')
    [k0 k1 r0 r1 Ntarget n1]

    testr=abs(r1-r0); 
end



% Mean hours worked
meanN = sum(n_pol(:).*dist(:));


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