function [x1,n1,a1,f1,indU1] = golden6_2(prodwage,rkap,cmin,cmax,expV)

global beta mu kap phi B

% Updated:
%   - to return the index of the nearest upper neighbour too (ind1).
%   - and to deal with vectors.

% Parameters
    %tol = optget('golden','tol',sqrt(eps));
    % Vector nkap
    siz = length(expV);
    % Tolerance
    tol=0.001;
    % Weights
    alpha1 = (3-sqrt(5))/2;
    alpha2 = (sqrt(5)-1)/2;
    % Distance
    d  = cmax-cmin;

% Functions
    % Labour supply
    labsup = @(c) (prodwage .* c.^(-mu) / B) .^ (1/phi);
    % Utility
    util = @(c,n) (c.^(1-mu))/(1-mu)-B*(n.^(1+phi))/(1+phi);
    % Savings
    aprime = @(c,n) prodwage.*n + rkap - c;

%% First guess

% Initial guess (consumption)
x1 = cmin+alpha1*d;

% Labour supply
n1 = min(labsup(x1),1);

% Savings
a1 = aprime(x1,n1);

% Current utility
u1 = util(x1,n1);

% Find closest upper element in grid for each element in a1
indU1 = zeros(1,siz); % initialise variable
for k=1:siz
    indU1(k) = find(kap>=a1(k),1,'first');
end

% In case of ind1=1... error with kap(ind1-1). Correct it with extra index for lower element:
indL = indU1-1;
indL(indL==0) = 2; 

% Interpolate expected continuation value
v1 = (a1-kap(indL)) ./ (kap(indU1)-kap(indL)).* expV(sub2ind(size(expV), [1:siz], indU1))+ (kap(indU1)-a1) ./ (kap(indU1)-kap(indL)).* expV(sub2ind(size(expV), [1:siz], indL)); % value of lower neighbour

% Value of the first guess
f1 = u1 + beta*v1;

% Bounds on a' 
if max(a1) > kap(end)
    f1(a1>kap(end)) = -1000 + a1(a1>kap(end));
    a1(a1>kap(end)) = kap(end);
elseif min(a1) < kap(1)
    f1(a1<kap(1)) = -1000 + a1(a1<kap(1));
    a1(a1<kap(1)) = kap(1);
end

%% Second guess
    
% Initial guess
x2 = cmin+alpha2*d; % consumption
n2 = min(labsup(x2),1); % labour supply
a2 = aprime(x2,n2); % savings
u2 = util(x2,n2); % current utility

% Find closest upper element in grid for each element in a2
indU2 = zeros(1,siz); % initialise variable
for k=1:siz
    indU2(k) = find(kap>=a2(k),1,'first');
end

indL = indU2-1;
indL(indL==0) = 2;

% Interpolate expected continuation value
v2 = (a2-kap(indL)) ./ (kap(indU2)-kap(indL)).* expV(sub2ind(size(expV), [1:siz], indU2))+ (kap(indU2)-a2) ./ (kap(indU2)-kap(indL)).* expV(sub2ind(size(expV), [1:siz], indL));

% Value of the first guess
f2 = u2 + beta*v2;

% Bounds on a' (need the if command to avoid error due to empty vector)
if max(a2) > kap(end)
    f2(a2>kap(end)) = -1000 + a2(a2>kap(end));
    a2(a2>kap(end)) = kap(end);
elseif min(a2) < kap(end)
    f2(a2<kap(1)) = -1000 + a2(a2<kap(1));
    a2(a2<kap(1)) = kap(1);
end


%% GOLDEN SEARCH

% Update distance
    d = alpha1*alpha2*d;

% Optimisation loop
while max(abs(d))>tol
    % Update distance
        d = d*alpha2;

    % if f2<f1:
        % Indicator: 1  in positions where f1>f2
            auxInd = f1 > f2;
        % x2 is new upper bound
            x2(auxInd) = x1(auxInd);
            f2(auxInd) = f1(auxInd);
            indU2(auxInd) = indU1(auxInd);
        % Update lower bound
            x1(auxInd) = x1(auxInd)-d(auxInd);
            n1 = min(labsup(x1),1); % labour supply
            a1 = aprime(x1,n1); % savings
            u1 = util(x1,n1); % current utility
            % Create vector that only includes the index of elements to
            % change
                indChange = [1:siz];
                indChange = indChange(auxInd);
            % Find them
            for k=indChange
                try
                    indU1(k) = find(kap>=a1(k),1,'first');
                catch % if a1(k) > maxkap
                    indU1(k) = length(kap);
                end
            end
            % Index for lower neighbour
            indL = indU1-1;
            indL(indL==0) = 2; % it does not matter, all weight goes to upper
            % Interpolate expected continuation value
            v1 = (a1-kap(indL)) ./ (kap(indU1)-kap(indL)) ...
                 .* expV(sub2ind(size(expV), [1:siz], indU1)) ...
                 + (kap(indU1)-a1) ./ (kap(indU1)-kap(indL)) ...
                 .* expV(sub2ind(size(expV), [1:siz], indL));
            % Value of the first guess
            f1 = u1 + beta*v1;
            % Bounds on a' (need the if command to avoid error due to empty vector)
            if max(a1) > kap(end)
                f1(a1>kap(end)) = -1000 + a1(a1>kap(end));
                a1(a1>kap(end)) = kap(end);
            elseif min(a1) < kap(end)
                f1(a1<kap(1)) = -1000 + a1(a1<kap(1));
                a1(a1<kap(1)) = kap(1);
            end

    % else, f2>f1:
        % Update indicator: 1  in positions where f2>f1
            auxInd = ~auxInd;
        % x1 is new lower bound
            x1(auxInd) = x2(auxInd);
            n1(auxInd) = n2(auxInd);
            a1(auxInd) = a2(auxInd);
            f1(auxInd) = f2(auxInd);
            indU1(auxInd) = indU2(auxInd);
        % Update upper bound
            x2(auxInd) = x2(auxInd) + d(auxInd);
            n2 = min(labsup(x2),1); % labour supply
            a2 = aprime(x2,n2); % savings
            u2 = util(x2,n2); % current utility
            % Create vector that only includes the index of elements to
            % change
                indChange = [1:siz];
                indChange = indChange(auxInd);
            % Find them
            for k=indChange
                try
                    indU2(k) = find(kap>=a2(k),1,'first');
                catch
                    indU2(k) = length(kap);
                end
            end
            % Index for lower neighbour
            indL = indU2-1;
            indL(indL==0) = 2; % it does not matter, all weight goes to upper
            % Interpolate expected continuation value
            v2 = (a2-kap(indL)) ./ (kap(indU2)-kap(indL)) ...
                 .* expV(sub2ind(size(expV), [1:siz], indU2)) ...
                 + (kap(indU2)-a2) ./ (kap(indU2)-kap(indL)) ...
                 .* expV(sub2ind(size(expV), [1:siz], indL));
            % Value of the first guess
            f2 = u2 + beta*v2;
            % Bounds on a' (need the if command to avoid error due to empty vector)
            if max(a2) > kap(end)
                f2(a2>kap(end)) = -1000 + a2(a2>kap(end));
                a2(a2>kap(end)) = kap(end);
            elseif min(a2) < kap(end)
                f2(a2<kap(1)) = -1000 + a2(a2<kap(1));
                a2(a2<kap(1)) = kap(1);
            end
end


%% RETURNS

% Return the larger of the two
    auxInd = f2>f1;
    x1(auxInd) = x2(auxInd);
    n1(auxInd) = n2(auxInd);
    a1(auxInd) = a2(auxInd);
    f1(auxInd) = f2(auxInd);
    indU1(auxInd) = indU2(auxInd);