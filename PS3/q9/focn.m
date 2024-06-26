function [foc] = focn(n,kap_vec,aprime, lambda, prodwage, r)


global mu  phi B tau2

    foc = B*n.^ (phi+tau2) - ...
        lambda*(1-tau2) .*prodwage.^(1-tau2) .*( (1+r)*kap_vec +...
        lambda*(prodwage.*n).^(1-tau2) - aprime).^(-mu);
    
end