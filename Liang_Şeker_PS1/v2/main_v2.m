
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kk kap nkap v dist maxkap minkap inckap

tic

%  set parameter values

mu     = 2;               % risk aversion. CRRA              
beta   = 0.95;            % subjective discount factor 
delta  = 0.1;             % depreciation
A      = 1;               % production technology
alpha  = 0.36;            % capital's share of income
b = 0.2;                  % BORROWING CONSTRAINT

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




% setting limits for interest rate search

%minrate = -delta;
%maxrate = (1-beta)/beta;
minrate = 0.02;
maxrate = 0.03;

% setting test values
testr=1;
r=testr;




% Choose Howard improvement for VFI
optH=1;

% we find the optimal interest rate by bisection

while (testr>0.001)
    
    r0 = 0.5*(maxrate-minrate)+minrate;%guess for interest rate
    k0=((r0+delta)/(alpha*A*labor^(1-alpha)))^(1/(alpha-1));%capital demand for that guess
    k1= aiyagari_v2(r0,optH); %calculating asset supply
    r1=alpha*A*max(0.001,k1)^(alpha-1)*labor^(1-alpha)-delta;%which interest rate implies that 

    %capital demand equals this value of capital supply 
    display('k0, k1, r guess, r implied  = ')
    [k0 k1 r0 r1]
     
    testr=abs(r1-r0);    
    %updating the interest rate
    if k1>k0
        maxrate=r0;
    else
        minrate=r0;
    end
     
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

k1
r1
toc

%%
fig = figure;

set(fig, 'Position', [100, 100, 1050, 1050]);  % [left, bottom, width, height]


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

saveas(fig,'aiyagari_v2.png');
