
function meank = aiyagari_v2(r0, optHoward)
% aiyagari.m is a function file which computes aggregate savings 
% given aggregate interest rate in Aiyagari's QJE 1994 paper

% r is interest rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global beta mu delta A alpha s N prob b kk kap v dist


%   write wage as a function of interest rate using Cobb-Douglas
%   production function

% r = r0;

wage = (1-alpha)*(A*(alpha/(r0+delta))^alpha)^(1/(1-alpha));

%--------------------------------------------------------------------                               
                                
                              
%   form capital grid
%   
   
maxkap = 26;                     % maximum value of capital grid  
minkap = -b;                     % borrowing constraint
inckap = 0.1;                    % size of capital grid increments
kap    = minkap:inckap:maxkap;   % state of assets 
nkap   = length(kap);            % number of grid points

%  initialize some variables
criterion = 10;
iter = 0;
   
%  iterate on Bellman's equation and get the decision 
%  rules and the value function at the optimum 


v   = zeros(N,nkap);

utilm = 190*zeros(N,nkap,nkap);

for j=1:N%loop over each skill s_t
    for i=1:nkap% loop over each possible a_t            
        cons = s(j)*wage + (1+r0)*kap(i) -kap';
        util = (cons.^(1-mu)-1)/(1-mu);
        ii = find( cons<= 0);              % infeasible consumption choice
        util(ii) = -10000;            
        utilm(j,i,:)=util;
    end
end

% Initialize matrices and distance
tv = 189*zeros(N,nkap); 
tdecis = tv; 
test    = 10;
criter_V = 1e-7;

%Solve the value function
while test > 0.001
    iter = iter + 1;

    % Optimisation
    for j=1:N                                    %each possible skill s_t
        for i=1:nkap                             % each asset level a_t             
            util=reshape(utilm(j,i,:),nkap,1); 
            vint = util' + beta*prob(j,:)*v;     % utility given j, i, for all possible k' next period       
           [tv(j,i),tdecis(j,i)] = max(vint);    %what is the optimal asset choice on the grid
        end
    end

    % Policy functions
    aprime_VFI = kap(tdecis);
    c_VFI = s'*wage + (1+r0)*kap - aprime_VFI;
    % Double-check: these two should be equal
        % s(4)*wage + (1+r0)*kap(15) - aprime(4,15)
        % c_VFI(4,15)

    % Howard
    if optHoward==1 && iter>3
        % Initialise variable
        dVV=1;
        v1how=tv;
        % Loop
        while dVV>criter_V
            % V(k) = u(c(k)) + beta*V(k')
            % v1how = (c_VFI.^(1-mu) - 1) / (1-mu) + beta*prob(j,:)*v1how(tdecis);

            for j = 1:N
                for i = 1:nkap
                    v1how(j,i)=(c_VFI(j,i)^(1-mu)-1)/(1-mu) + beta*prob(j,:)*v1how(:,tdecis(j,i));
                end
            end

            dVV=max(max(abs(v1how-tv)));
            %disp(dVV)
            tv=v1how;
        end
    end

    test=max(max(abs(tv-v))); % distance
    v=tv;
end
  
 
display('done with VF iteration')


% Compute the measure

% Slow method: with simulations

rng(12345) % fix seed
nsim = 7500;  % number of unique HH's we simulate

asim=(floor(1/3*maxkap)+3)*ones(nsim,1);  % integer initial HH asset grid
zsim=(floor(N/2)+1)*ones(nsim,1); % integer initial z-shocks

err_sim=1;
while err_sim>1e-4
    
    unidraw = rand(nsim,1);
    
    zsim_prime = 0*zsim; asim_prime = 0*asim;
    
    for i=1:nsim
        for zp=1:N        
            if ( unidraw(i)<= sum( prob(zsim(i),1:zp)) )
                zsim_prime(i) = zp;
                break
            end
        end
    end
    
    zsim=zsim_prime;  % value of z-shock at beginning of time t
    
    for i=1:nsim
        asim_prime(i) = tdecis(zsim(i),asim(i)) ;
    end
    
    meank = mean(kap(asim_prime));
    
    err_sim = max(abs(mean(kap(asim_prime))-mean(kap(asim))));
    asim = asim_prime; % update savings for use at t+1

end


dist = ones(N,nkap); kk=[];
for p = 1:N
    for i = 1:nkap
        dist(p,i) = sum((zsim(:)==(p)).*asim(:)==i);
        kk(p,i) = kap(tdecis(p,i));
    end
end
dist = dist/nsim;
meank2 = sum(kk(:).*dist(:));



display('done with measure')

