function [real_i] = fisher_eq( nominal_i,inflation )
   
    real_i = (1+nominal_i)./(1+inflation)-1;
    
end