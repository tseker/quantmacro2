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

N        = 6;             % number of discretized states
rho      = 0.9;           % first-order autoregressive coefficient
sigma    = 0.2;           % standard deviation of error_t

% prob is transition matrix of the Markov chain
% logs is the discretized states of log labor earnings
% invdist is the invariant distribution of Markov chain


% s productivity

[logs,prob] = rouwenhorst(N,0,rho,sigma); logs = logs';

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





% Q2: Choose Howard improvement for VFI
optH=0;

% we find the optimal interest rate by bisection
r0 = 0.5*(maxrate-minrate)+minrate;%guess for interest rate

while (testr>0.001)
    k0=((r0+delta)/(alpha*A*labor^(1-alpha)))^(1/(alpha-1));%capital demand for that guess
    k1= aiyagari2_v3(r0); %calculating asset supply
    r1=alpha*A*max(0.001,k1)^(alpha-1)*labor^(1-alpha)-delta;%which interest rate implies that 

    
    %Q4: Implied interest rate:
    r0 = 0.8*r0 + 0.2*r1;  %Guess the interest rate

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


%% sort by wealth

col = sum(dist);
quintind = zeros(5,1);
rip = cumsum(col);

for i = 1:5
    quintind(i) = find(rip <= 0.20 * i, 1, 'last');
end

wealth = kap.*col;
wealthquint = zeros(5,1);
wealthcum = cumsum(wealth);
wealthquint(1) = wealthcum(quintind(1));
for i=2:5
    wealthquint(i) = wealthcum(quintind(i))- wealthcum(quintind(i-1));
end
share = wealthquint./(sum(wealthquint));
%mean income by wealth quantile 
inc = wage*s'; % + (1+r1)*repmat(kap,N,1);
incdistrib = inc.*dist;
ripdistrib = cumsum(col);
probquint = zeros(5,1);
probquint(1) = ripdistrib(quintind(1));
for i=2:5
   probquint(i) = ripdistrib(quintind(i)) - ripdistrib(quintind(i-1));
end
colinc = sum(incdistrib,1);
ripincdistrib = cumsum(colinc);
sumincquint = zeros(5,1);
sumincquint(1) = ripincdistrib(quintind(1));
for i=2:5
   sumincquint(i) = ripincdistrib(quintind(i)) - ripincdistrib(quintind(i-1));
end
meaninc = sumincquint./probquint;
correlation = corr(meaninc,wealthquint);

%%
% Plotting the share of wealth for each quintile without numbers on top
f6 = figure('visible','on');
bar(1:5, share, 'b');
title('Share of Wealth by Wealth Quintile');
xlabel('Wealth Quintile');
ylabel('Share of Total Wealth');
xticks(1:5);
yticks(0:0.1:1);
grid on;




%%
fig=figure;

set(fig, 'Position', [100, 100, 900, 900]);  % [left, bottom, width, height]


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

% subplot(3,2,5);
% plot(kap,[v(1,:)' v(3,:)' v(5,:)'])
% ylabel('value function')
% xlabel('asset holdings')
% hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','southeast');

saveas(fig,'aiyagari_v6.png');


