function [newV] = newV(aprime, W_matrix, j,wage, r0 )

global beta mu delta A alpha s N prob b kk kap v dist
    % vectorization of loop for asset
    % aprime is a vector of the same size as kap
    c = s(j)*wage + (1+r0).*kap -aprime;
    util = (c.^(1-mu)-1)./(1-mu);
    if c <=0
        util = -inf;
    end
        
    newV = util + W_func(aprime);
end