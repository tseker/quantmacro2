%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kk kap v dist

tic

%  set parameter values

mu     = 2;               % risk aversion. CRRA              
beta   = 0.95;            % subjective discount factor 
delta  = 0.1;             % depreciation
A      = 1;               % production technology
alpha  = 0.36;            % capital's share of income

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

%minrate = -delta;
%maxrate = (1-beta)/beta;
minrate = 0.02;
maxrate = 0.03;

% setting test values
testr=1;

% INITIALIZE V (globals)

% we find the optimal interest rate by bisection
r0 = 0.5*(maxrate-minrate)+minrate;%guess for interest rate
while (testr>0.001)
    k0=((r0+delta)/(alpha*A*labor^(1-alpha)))^(1/(alpha-1)); %capital demand for that guess
    k1= aiyagari_v6(r0); %calculating asset supply
    r1=alpha*A*max(0.001,k1)^(alpha-1)*labor^(1-alpha)-delta;%which interest rate implies that 

    %capital demand equals this value of capital supply 
    display('r guess, r implied  = ')
    [r0 r1]

    %pause
     
    testr=abs(r1-r0); 
    %updating the interest rate
    r0 = 0.9*r0 + 0.1*r1;%guess for interest rate
end

nkap=length(kap);
kopt=reshape(kk,N,nkap);

%
wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));
for j=1:N
    for i=1:nkap              
        cons(j,i) = s(j)*wage + (1+r0)*kap(i) -kopt(j,i);
    end
end

figure(1);

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
plot(kap,[cons(1,:)' cons(3,:)' cons(5,:)'])
ylabel('consumption')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northwest');

subplot(3,2,5);
plot(kap,[v(1,:)' v(3,:)' v(5,:)'])
ylabel('value function')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','southeast');

toc

kap0 = -b:0.1:26;
figure(2)
hold on
    plot(kap,1/3*ones(1,length(kap)),'or')
    plot(kap0,2/3*ones(1,length(kap0)), 'ob')
    ylim([0 1])
    yticks([])
hold off