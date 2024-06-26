%% PS 
clear
% close all figures
close all
% change location to the folder where this m file is saved
mfile_name          = mfilename('fullpath');  % This line of code retrieves the full path of the current MATLAB file and assigns it to the variable mfile_name. The mfilename function with the argument 'fullpath' returns the full path of the currently executing file, including the file name and extension. The full path is then assigned to the variable mfile_name.
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

global par T

par.beta = 0.99;    % discount
par.gamma = 2;      % risk aversion
par.psi = 1;        % frisch elasticity
par.kappa = 0.1;    % NKPC kappa
par.sigma_epsilon = 0.1;    % idiosyncratic shock s.d.
par.rho = 0.96;     % idiosyncratic producitivyt shock persistence
par.b = 0;          % borrowing limit
par.mu = 2;         % labor union markup
par.psi_pi = 1.5;    % Taylor rule
par.rhonu = 0.5;    % monetary shock persistence
par.sigma_nu = sqrt(0.0001);  % monetary policy shock s.d.
par.A = 1;           % aggregate technology
par.L = 1;          % aggregate labor supply



par.N = 5; % states of idiosyncratic labor income shocks

par.B = 2;          % government bond supply

% asset grid
par.maxkap =70;                     % maximum value of capital grid  
par.minkap = -par.b;                     % borrowing constraint                    
par.nkap = 200;                      % number of asset grids
par.kap = logspace(0,log10(par.maxkap+1-par.minkap),par.nkap)-1+par.minkap;

par.T = 100;
par.monetary_shock = 0.0025* par.rhonu.^(0:(T-1));

% idiosyncratic shock
[logs,par.prob] = tauchen(par.N,0,par.rho,par.sigma_epsilon,3); logs = logs';

        % prob contains probability of transition from s0 (row) to s1 (column)
 % Compute invariant distribution
invdist = ones(1,par.N)/par.N; test = 1; 

while (test>0.0000001)
    invdist2 = invdist*par.prob;
    test = max(abs(invdist2-invdist));
    invdist = invdist2;
end
invdist = invdist';

 par.s = exp(logs);
par.tech_labor = par.s*invdist; % technology augmented labor supply
 
 
 
%% General Equilibrium: Steady-State
r_max = 0.1;
r_min = 0.00001;
tol = 1e-6;

% bisection to impose asset market clearing condition
while abs(r_max-r_min)>tol
   
   r_guess = 0.5*r_max + 0.5*r_min;
   [ss_asset_mkt, ss_policy, ss_dist] = bewley_ss(r_guess, par );
   
   if ss_asset_mkt>0
       r_max = r_guess;
   else
       r_min = r_guess;
   end
   
end
% 
% %% figures
% 
% fig = figure;
% 
% 
% plot(par.kap,[ss_dist(1,:)',ss_dist(3,:)', ss_dist(5,:)'])
% ylabel('distribution')
% xlabel('asset holdings')
% hleg1 = legend('Asset Distribution - Lowest Skill Shock',  'Asset Distribution - Intermediate Skill Shock','Asset Distribution - Highest Skill Shock');
% 
% 

ss_R = r_guess;
tau = r_guess*par.B/(par.A*par.L*par.tech_labor);
par.phi = 1/par.mu*(1-tau)*(par.A*par.tech_labor)^(-par.gamma);




%% Question 2










