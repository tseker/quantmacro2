function [value_q5,auxInterp] = value_q5(aprime,inda,inds,r,wage)

% The function takes:
% - a scalar with initial assets (a),
% - a scalar for the productivity state (s),
% - a scalar for the savings decision (aprime),
% - prices (r, wage),
% - and some globals.

% The function returns:
% - expected value.
% - interpolation information (neighbouring elements and weights).

global beta mu s N kap prob v

% Consumption: (from budget constraint)
    cons = s(inds)*wage + (1+r)*kap(inda) - aprime;

% Utility:
    util = (cons^(1-mu)-1)/(1-mu);
    util(cons<=0) = -10000; % infeasible consumption choice

% Value function
    % Interpolate for first state and save weights
    [auxV(1,1),auxInterp] = linInterp(aprime,kap,v(1,:));
    % We use the weights previously computed to save time:
    for sprime=2:N
        auxV(sprime,1) = linInterp(aprime,kap,v(sprime,:),auxInterp);
        %auxV(sprime,1) = interp1(kap,v(sprime,:),aprime);
    end

    %[vprime,auxInterp] = linInterp(aprime,agrid,v);
    value_q5 = util + beta * prob(inds,:)*auxV;