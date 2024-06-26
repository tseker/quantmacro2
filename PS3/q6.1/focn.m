function [foc] = focn(sj,n, wage, r0)


global beta mu A delta alpha s b N prob kk kap v dist phi B n_pol c_pol Ntarget

    foc = B*n.^(phi)-(wage*sj*((1+r0)*kap' + wage*sj*n - kap).^(-mu));

end