function [error ] = NKPC_error(labor_t,par,wage_inflation,real_i)

    % tau = r*par.B/(par.A*par.L*par.tech_labor);
    wage_inflation_1 = [wage_inflation(2:end),0];
    tau_t = real_i.* par.B ./(par.A*par.L*par.tech_labor);
    
    error = wage_inflation - par.kappa.*...
            (  par.phi .* (labor_t).^(par.psi) - ...
        1/par.mu.*(1-tau_t).*par.A.*(labor_t * par.tech_labor).^(-par.gamma)    )-...
            par.beta .*wage_inflation_1;

end