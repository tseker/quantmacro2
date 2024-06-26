function [Z,TP] = rouwenhorst(N,mu,rho,sigma)
% Calculates the transition probabilities matrix for a Rouwenhorst approximation
% of an AR(1) process with persistence parameter rho and standard deviation of
% the error term sigma.
% Finds a Markov chain whose sample paths approximate those of
% the AR(1) process
%                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)

% m : the number of states, 
% exp : the expectation, 
% rho : the auto-correlation coefficient, 
% sigma : the standard deviation of the error term epsilon.

%Output:     Z       N*1 vector, nodes for Z
%            PI      N*N matrix, transition probabilities

% Advantage of Rouwenhorst: better match the un-&-conditional variances and
% autocorrelation of the AR(1) process.
% Disadvantage: Errors further away from the normal distribution.



% Calculate the standard deviation of the underlying process
s = sqrt(sigma^2 / (1 - rho^2));


% Calculate the transition probabilities
p = (1+rho) /2;
TP=[ p  (1-p);
        (1-p) p];
    
    for n=3:N
    a1R=[TP zeros(n-1, 1); zeros(1, n)];
    a2R=[zeros(n-1, 1) TP; zeros(1,n)];
    a3R=[zeros(1,n); TP zeros(n-1,1)];
    a4R=[zeros(1,n); zeros(n-1,1) TP];
    TP=p*a1R+(1-p)*a2R+(1-p)*a3R+p*a4R;
    TP(2:end-1, :) = TP(2:end-1, :)/2;
    end

f= sqrt(N-1)*s;
Z = linspace(-f,f,N)';
Z = Z+mu;

