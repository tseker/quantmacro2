function [newV] = newV(aprime, i, j, wage, r0 )

global beta mu delta A alpha s N prob b kk kap v dist

    c = s(j)*wage + (1+r0)*kap(i) -aprime;
    util = (c^(1-mu)-1)/(1-mu);
    if c <=0
        util = -inf;
    end
        
    temp = kap - aprime;
    temp_pos = min(temp(temp>=0));
    temp_neg = max(temp(temp<0));
    
    index2 = find(temp == temp_pos);
    index1 = find(temp==temp_neg);
    weight = (aprime-kap(index1))/(kap(index2)-kap(index1));
    tempV =  weight * v(:,index2) + (1-weight)* v(:,index1);
    
    newV = util + beta * prob(j,:)* tempV;
end