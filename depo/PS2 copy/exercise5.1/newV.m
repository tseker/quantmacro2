function [newV] = newV(aprime, W_func, j,wage, r0,nopt )

global beta mu delta A alpha s N prob b kk kap v dist
    % vectorization of loop for asset
    % aprime is a vector of the same size as kap
    c = s(j)*wage*nopt(j,i,:) + (1+r0).*kap -aprime;
    util = (c.^(1-mu)-1)./(1-mu)-B*(n^(1+phi)/(1+phi));
    if c <=0
        util = -inf;
    end
        
    newV = util + W_func(aprime);
end