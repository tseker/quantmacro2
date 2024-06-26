function [labor_t] = NKPC(wage_inflation, par,real_i)
    
    fun = @(x) NKPC_error(x,par,wage_inflation,real_i);
    labor_t = fsolve(fun, [1,1]);
    
end