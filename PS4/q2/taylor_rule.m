function [nominal_i] = taylor_rule(inflation, par, ss_R)
    
    nominal_i = (1+ss_R) .* (1+inflation).^par.psi_pi.*par.monetary_shock -1;
    
end