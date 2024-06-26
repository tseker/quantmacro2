function [result] = first_dev_focn(sj,n, wage, r0)

global beta mu A delta alpha s b N prob kk kap v dist phi B n_pol c_pol Ntarget

   result = phi*B*n.^(phi-1)+mu* wage.^2 * sj.^2 *((1+r0)*kap' + wage*sj*n - kap).^(-mu-1);
   
end