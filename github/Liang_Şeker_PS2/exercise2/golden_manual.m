function [x1,f1] = golden_manual(W_matrix, j,wage,r0,a,b)

global beta mu delta A alpha s N prob  kk kap v dist nkap



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
ind1 = zeros(1,nkap);
for k=1:nkap
    ind1(k) = find(kap>=x1(k),1,'first');
end
f1 = util1 + (x1-kap(ind1-1)) ./ (kap(ind1)-kap(ind1-1)).* W_matrix(j,ind1)+ (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* W_matrix(j,ind1-1);


c2 = s(j)*wage + (1+r0).*kap -x2;
util2 = (c2.^(1-mu)-1)./(1-mu);
util2(c2<0)= -inf;
ind2 = zeros(1,nkap);
for k=1:nkap
    ind2(k) = find(kap>=x2(k),1,'first');
end
f2 = util2 + (x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)).* W_matrix(j,ind2)+ (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)).* W_matrix(j,ind2-1);

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
    for k=1:nkap
        ind1(k) = find(kap>=x1(k),1,'first');
    end
    f1 = util1 + (x1-kap(ind1-1)) ./ (kap(ind1)-kap(ind1-1)).* W_matrix(j,ind1)+ (kap(ind1)-x1) ./ (kap(ind1)-kap(ind1-1)).* W_matrix(j,ind1-1);


    c2 = s(j)*wage + (1+r0).*kap -x2;
    util2 = (c2.^(1-mu)-1)./(1-mu);
    util2(c2<0)= -inf;
    for k=1:nkap
        ind2(k) = find(kap>=x2(k),1,'first');
    end
    f2 = util2 + (x2-kap(ind2-1)) ./ (kap(ind2)-kap(ind2-1)).* W_matrix(j,ind2)+ (kap(ind2)-x2) ./ (kap(ind2)-kap(ind2-1)).* W_matrix(j,ind2-1);

  
end

% Return the larger of the two
 x1(f2>f1) = x2(f2>f1); 
 f1(f2>f1) = f2(f2>f1); 
 
end

