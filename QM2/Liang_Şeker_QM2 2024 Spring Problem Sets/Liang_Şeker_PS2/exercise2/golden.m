function [x1,f1,indU] = golden(W_matrix, j,wage,r0,a,b)

global beta mu delta A alpha s N prob  kk kap v dist nkap

    W_func = griddedInterpolant(kap,W_matrix(j,:));


% lower the tolerance compared to previous case to increase the speed
tol = 0.00001;

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
d  = b-a;
x1 = a+alpha1*d;
x2 = a+alpha2*d;


c1 = s(j)*wage + (1+r0).*kap -x1;
util1 = (c1.^(1-mu)-1)./(1-mu);
util1(c1<0)= -inf;
f1 = util1 + W_func(x1);


c2 = s(j)*wage + (1+r0).*kap -x2;
util2 = (c2.^(1-mu)-1)./(1-mu);
util2(c2<0)= -inf;
f2 = util2 + W_func(x2);

d = alpha1*alpha2.*d;
while max(d)>tol
    d = d.*alpha2;
    x1_temp = x1;
    x2_temp = x2;
    x1 = (x1_temp-d).*(f1-f2>0) + x2_temp.* (f1-f2<=0);
    x2 =  x1_temp.*(f1-f2>0) + (x2_temp+d).*(f1-f2<=0);

    c1 = s(j)*wage + (1+r0).*kap -x1;
    util1 = (c1.^(1-mu)-1)./(1-mu);
    util1(c1<0)= -inf;
    f1 = util1 + W_func(x1);


    c2 = s(j)*wage + (1+r0).*kap -x2;
    util2 = (c2.^(1-mu)-1)./(1-mu);
    util2(c2<0)= -inf;
    f2 = util2 + W_func(x2);
  
end

% Return the larger of the two
 x1(f2>f1) = x2(f2>f1); 
 f1(f2>f1) = f2(f2>f1); 
 indU = zeros(1,nkap);
     for k = 1:nkap
         indU(k) = find(kap>=x1(k),1,'first');
     end
end

