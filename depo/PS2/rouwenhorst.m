function [Z,P] = rouwenhorst(~,~,rho,sigma,m)
% Calculates the transition probabilities matrix for a Rouwenhorst approximation
% of an AR(1) process with persistence parameter rho and standard deviation of
% the error term sigma, using m points in the discrete approximation.
% m is the number of states, exp is the expectation, rho is the auto-correlation coefficient, sigma is the standard deviation of the error term epsilon.

m = round(m); % ensure that m is an integer

% Calculate the standard deviation of the underlying process
s = sqrt(sigma^2 / (1 - rho^2));

% Calculate the first and last points
y1 = -s*sqrt(m-1) + exp;
yM = s*sqrt(m-1) + exp;

% Calculate the intermediate points
Z = y1:(yM-y1)/(m-1):yM;

% Calculate the transition probabilities

p = (1+rho) /2;
q=p;

P_Rouw=[ p  (1-p);
        (1-q) q];
    
    for i_R=2:(m-1)
    a1R=[P_Rouw zeros(i_R, 1); zeros(1, i_R+1)];
    a2R=[zeros(i_R, 1) P_Rouw; zeros(1, i_R+1)];
    a3R=[zeros(1,i_R+1); P_Rouw zeros(i_R,1)];
    a4R=[zeros(1,i_R+1); zeros(i_R,1) P_Rouw];
    P_Rouw=p*a1R+(1-p)*a2R+(1-q)*a3R+q*a4R;
    P_Rouw(2:i_R, :) = P_Rouw(2:i_R, :)/2;
    end

    
P = P_Rouw;

