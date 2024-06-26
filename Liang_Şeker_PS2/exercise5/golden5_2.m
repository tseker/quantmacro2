function [x1,n1_int,f1,ind1] = golden5_2(r0,prodwage,amin,amax,prob,nopt)

global beta mu kap v phi B

util = @(cons,n) (cons.^(1-mu))/(1-mu)-B*(n.^(1+phi))/(1+phi);


    % Vector size
    nkap = length(kap);
    % Tolerance
    tol=1e-7;
    % Weights
    alpha1 = (3-sqrt(5))/2;
    alpha2 = (sqrt(5)-1)/2;
    % Distance
    d  = amax-amin;
    % Initial guesses
    x1 = amin+alpha1*d;
    x2 = amin+alpha2*d;
    
%%

% Evaluating first guess
    %  find the index of the first element in the vector kap that is greater than or equal to the value x1(k)
    for k=1:nkap
        ind1(k) = find(kap>=x1(k),1,'first');
    end

    % Interpolate value function
    vf_int = (x1-kap(ind1-1)) ./ (kap(ind1)-kap(ind1-1)).* v(:,ind1)+ (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* v(:,ind1-1); 
             % weight of upper neighbour*value of upper neighbour + weight of lower neighbour*value of lower neighbour

    % Interpolate labor supply
    n1_int = (x1-kap(ind1-1))./ (kap(ind1)-kap(ind1-1)).* nopt(sub2ind(size(nopt),1:nkap,ind1)) + (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* nopt(sub2ind(size(nopt),1:nkap,ind1-1));

    % Cash in hand
    cih = prodwage*n1_int + (1+r0)*kap;

    % Utility
    u1 = util(max(cih-x1,0.01),n1_int);
    
    % Expectation
    v1 =  sum(prob' .* vf_int);

    % Value of the first guess
    f1 = u1 + beta*v1;


% Evaluating second guess (analogous)
% Evaluate the second guess x2 (analogous)
    for k=1:nkap
        ind2(k) = find(kap>=x2(k),1,'first');
    end
    vf_int=(x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)) .* v(:,ind2)+ (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)) .* v(:,ind2-1);
    n2_int = (x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)).* nopt(sub2ind(size(nopt),1:nkap,ind2)) + (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)).* nopt(sub2ind(size(nopt),1:nkap,ind2-1));
    u2 = util(max(cih-x2,0.01),n2_int);
    v2 =  sum(prob' .* vf_int);
    f2 = u2 + beta*v2;

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
            ind2(auxInd) = ind1(auxInd);
        % Update lower bound
            x1(auxInd) = x1(auxInd)-d(auxInd);
            % Create vector that only includes the index of elements to change
                indChange = [1:nkap];
                indChange = indChange(auxInd);
            % Find them
            for k=indChange
                ind1(k) = find(kap>=x1(k),1,'first');
            end

            vf_int = (x1-kap(ind1-1)) ./ (kap(ind1)-kap(ind1-1)).* v(:,ind1)+ (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* v(:,ind1-1); 
            n1_int = (x1-kap(ind1-1))./ (kap(ind1)-kap(ind1-1)).* nopt(sub2ind(size(nopt),1:nkap,ind1)) + (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* nopt(sub2ind(size(nopt),1:nkap,ind1-1));
            cih = prodwage*n1_int + (1+r0)*kap;
            u1 = util(max(cih-x1,0.01),n1_int);
            v1 =  sum(prob' .* vf_int);
            f1 = u1 + beta*v1;

    % else, f2>f1:
        % Update indicator: 1  in positions where f2>f1
            auxInd = ~auxInd;
        % x1 is new lower bound
            x1(auxInd) = x2(auxInd);
            f1(auxInd) = f2(auxInd);
            ind1(auxInd) = ind2(auxInd);
        % Update upper bound
            x2(auxInd) = x2(auxInd) + d(auxInd);
            % Create vector that only includes the index of elements to change
                indChange = [1:nkap];
                indChange = indChange(auxInd);
            % Find them
            for k=indChange
                ind2(k) = find(kap>=x2(k),1,'first');
            end

        vf_int=(x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)) .* v(:,ind2)+ (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)) .* v(:,ind2-1);
        n2_int = (x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)).* nopt(sub2ind(size(nopt),1:nkap,ind2)) + (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)).* nopt(sub2ind(size(nopt),1:nkap,ind2-1));
        cih = prodwage*n2_int + (1+r0)*kap;
        u2 = util(max(cih-x2,0.01),n2_int);
        v2 =  sum(prob' .* vf_int);
        f2 = u2 + beta*v2;
end

% Return the larger of the two
    auxInd = f2>f1;
    x1(auxInd) = x2(auxInd);
    n1_int(auxInd) = n2_int(auxInd);
    f1(auxInd) = f2(auxInd);
    ind1(auxInd) = ind2(auxInd);