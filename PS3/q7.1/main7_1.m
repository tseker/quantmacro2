%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%				Aiyagari's model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

global beta mu A delta alpha s b N prob kk c_pol kap nkap v dist

tic

%% PRELIMINARIES

%  Parameter values
    mu     = 2;               % risk aversion. CRRA              
    beta   = 0.95;            % subjective discount factor 
    delta  = 0.1;             % depreciation
    A      = 1;               % production technology
    alpha  = 0.36;            % capital's share of income
    b      = 0.2;             % Borrowing constraint

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
    labor = s*invdist;
%% GE

% Loop on r
    % Initial bounds
        minrate = 0.02;
        maxrate = 0.03;
    % Initial guess for interest rate
        r0 = 0.5*(maxrate-minrate)+minrate;
    % Initialise distance test
        testr=1;
while (testr>1e-8)
    k0=((r0+delta)/(alpha*A*labor^(1-alpha)))^(1/(alpha-1)); %capital demand for that guess
    k1= aiyagari7(r0); %calculating asset supply
    r1=alpha*A*max(0.001,k1)^(alpha-1)*labor^(1-alpha)-delta;%which interest rate implies that 

    %capital demand equals this value of capital supply 
    display('r guess, r implied  = ')
    [r0 r1]

    testr=abs(r1-r0); 
    %updating the interest rate
    r0 = 0.9*r0 + 0.1*r1;%guess for interest rate
end

nkap=length(kap);

% Other variables
    % Wage
    w0 = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));
    % Policy function: consumption
    c_pol = s' *w0 + (1+r0)*kap -kk;

disp("Done with steady state computations")

%% TFP SHOCK

% Parameters
T       = 100;                                                          % number of periods
shock   = 1.01;                                                         % initial shock
rho     = 0.9;                                                          % shock persistence

% TFP sequence
aux_t   = [0:(T-1)];                                                    % aux variable for t-1
A_t     = 1 - rho.^aux_t + rho.^aux_t * shock;                          % overall TFP shock; rho.^aux_t : decay of the shock over time; 

% Initialise variables
testr           = 1;
tolr            = 1e-5;
r_t             = r0*ones(1,T);                                         % initial guess for r
r_t(1)          = alpha*A_t(1)*(k1/labor)^(alpha-1)-delta;

% Finding the sequence
while testr > tolr
    % Implied variables
    k_imp   = deterministicseq7(r_t,A_t,T,false);
    r_imp   = alpha * A_t .* max(0.001,k_imp).^(alpha-1) * labor^(1-alpha) - delta;
    % Update test and guess
    testr   = max(abs(r_t-r_imp));
    r_t(2:end) = 0.7 * r_t(2:end) + 0.3 * r_imp(2:end);
end

% Value right after the shock
[~,v_t] = deterministicseq7(r_t,A_t,T,true);
v0 = reshape(v_t(:,1),nkap,N)';

% Consumption equivalent
c_eq = (v0 ./ v).^(1/(1-mu)) - 1;
c_eq_mean = sum(c_eq(:).*dist(:));
fprintf("The average consumption equivalent is %3.2f%%.\n",100*c_eq_mean);
toc

%%

fig=figure(1);

subplot(3,2,1);
plot(kap,[kap' kk(1,:)' kk(3,:)' kk(5,:)'])
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
plot(kap,[v(1,:)' v(3,:)' v(5,:)'])
ylabel('value function')
xlabel('asset holdings')
hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','southeast');
saveas(fig,'main_71.png');

%%
%% GRAPHS (deterministic path)

w_t     = (1-alpha)*(A_t.*(alpha./(r_t+delta)).^alpha).^(1/(1-alpha));  % wage sequence

fig= figure(2); % Transitions

% General options
sim=1; % take the first simulation
fsize=12;
lwidth=1.5;

% Capital
subplot(3,1,1)
    hold on
    plot(k_imp, 'LineWidth', lwidth) % IRF
    yline(k1,'--','LineWidth', 2/3*lwidth) % SS value
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
    yline(r0,'--','LineWidth', 2/3*lwidth); % SS value
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
    yline(w0,'--','LineWidth', 2/3*lwidth); % SS value
    % Options
    title("Path of wages")
    grid on
    xlabel('Time')
    ylabel('Wage')
    set(gca,'FontSize',fsize)
    hold off

saveas(fig,'detpath_7.png');


fig2= figure(3); % Consumption equivalents
    plot(kap,[c_eq(1,:)' c_eq(3,:)' c_eq(5,:)'],'LineWidth',lwidth)
    % Options
    grid on
    ylabel('consumption equivalent')
    xlabel('asset holdings')
    set(gca,'FontSize',fsize)
    hleg1 = legend('Lowest Skill Shock', 'Intermediate Skill Shock', 'Highest Skill Shock','Location','northeast');

    saveas(fig2,'cev7.png');
