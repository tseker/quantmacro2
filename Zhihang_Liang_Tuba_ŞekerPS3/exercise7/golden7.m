function [x1,f1,indU1] = golden7(coh,expV,amin,amax)

global beta mu kap

% Updated:
%   - to return the index of the nearest upper neighbour too (ind1).
%   - and to deal with vectors.
% Warning: coh, amin and amax should be vectors of the same size.


%% PRELIMINARIES

% Parameters
    %tol = optget('golden','tol',sqrt(eps));
    % Vector size
    siz = length(coh);
    % Tolerance
    tol=0.001;
    % Weights
    alpha1 = (3-sqrt(5))/2;
    alpha2 = (sqrt(5)-1)/2;
    % Distance
    d  = amax-amin;
    % Initial guesses
    x1 = amin+alpha1*d;
    x2 = amin+alpha2*d;

% Functions
    util = @(c) (c.^(1-mu)) / (1-mu);


%% FIRST GUESS

% Current utility
u1 = util(max(coh-x1,0.01));
indU1 = zeros(size(coh));
% Find closest upper element for each element in x1
for k=1:siz
    indU1(k) = find(kap >= x1(k),1,'first');
end
% Index for closest lower element
indL1 = indU1-1;
indL1(indL1==0) = 2; % to avoid error (searching position 0 of vector)
% Interpolate expected continuation value
v1 = (x1-kap(indL1)) ./ (kap(indU1)-kap(indL1)) ... % weight of upper neighbour
     .* expV(sub2ind(size(expV), [1:siz], indU1)) ... % value of upper neighbour
     + (kap(indU1)-x1) ./ (kap(indU1)-kap(indL1)) ... % weight of lower neighbour
     .* expV(sub2ind(size(expV), [1:siz], indL1)); % value of upper neighbour
% Value of the first guess
f1 = u1 + beta*v1;


%% SECOND GUESS (analogous)

u2 = util(max(coh-x2,0.01));
indU2 = zeros(size(coh));
for k=1:siz
    indU2(k) = find(kap >= x2(k), 1,'first');
end
% Index for closest lower element
indL2 = indU2-1;
indL2(indL2==0) = 2;
% Interpolate expected continuation value
v2 = (x2-kap(indL2)) ./ (kap(indU2)-kap(indL2)) .* expV(sub2ind(size(expV), [1:siz], indU2)) ...
     + (kap(indU2)-x2) ./ (kap(indU2)-kap(indL2)) .* expV(sub2ind(size(expV), [1:siz], indL2));
% Value of the first guess
f2 = u2 + beta*v2;


%% GOLDEN SEARCH

% Update distance
    d = alpha1*alpha2*d;

% Optimisation loop
while max(abs(d))>tol
% Update distance
    d = d*alpha2;
% Indicator: 1  in positions where f1>f2
    auxInd = f1 > f2;
% if f2<f1:
    % x2 is new upper bound
        x2(auxInd) = x1(auxInd);
        f2(auxInd) = f1(auxInd);
        indU2(auxInd) = indU1(auxInd);
    % Update lower bound
        x1(auxInd) = x1(auxInd)-d(auxInd);
        u1 = util(max(coh-x1,0.01));
        % Create vector that only includes the index of elements to
        % change
            indChange = [1:siz];
            indChange = indChange(auxInd);
        % Find them
        for k=indChange
            indU1(k) = find(kap>=x1(k),1,'first');
        end
        indL1 = indU1-1;
        indL1(indL1==0) = 2;
        % Interpolate expected continuation value
        v1 = (x1-kap(indL1)) ./ (kap(indU1)-kap(indL1)) .* expV(sub2ind(size(expV), [1:siz], indU1)) ...
             + (kap(indU1)-x1) ./ (kap(indU1)-kap(indL1)) .* expV(sub2ind(size(expV), [1:siz], indL1));
        % Compute value
        f1 = u1 + beta*v1;
% else, f2>f1:
    % Update indicator: 1  in positions where f2>f1
        auxInd = ~auxInd;
    % x1 is new lower bound
        x1(auxInd) = x2(auxInd);
        f1(auxInd) = f2(auxInd);
        indU1(auxInd) = indU2(auxInd);
    % Update upper bound
        x2(auxInd) = x2(auxInd) + d(auxInd);
        u2 = util(max(coh-x2,0.01));
        % Create vector that only includes the index of elements to
        % change
            indChange = [1:siz];
            indChange = indChange(auxInd);
        % Find them
        for k=indChange
            indU2(k) = find(kap>=x2(k),1,'first');
        end
        indL2 = indU2-1;
        indL2(indL2==0) = 2;
        % Interpolate expected continuation value
        v2 = (x2-kap(indL2)) ./ (kap(indU2)-kap(indL2)) .* expV(sub2ind(size(expV), [1:siz], indU2)) ...
             + (kap(indU2)-x2) ./ (kap(indU2)-kap(indL2)) .* expV(sub2ind(size(expV), [1:siz], indL2));
        % Compute value
        f2 = u2 + beta*v2;
end


%% RETURNS

% Return the larger of the two
    auxInd = f2>f1;
    x1(auxInd) = x2(auxInd);
    f1(auxInd) = f2(auxInd);
    indU1(auxInd) = indU2(auxInd);