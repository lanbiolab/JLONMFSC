function [Z,S,U,E,obj_error] = LRR(X,Z_ini,H_ini,lambda1,lambda2,max_iter,Ctg)
% subcode of LRRNMF_main
[m,n] = size(X);
% ---------- Initilization -------- %
miu = 1e-2;
rho = 1.3;
max_miu = 1e8;
tol  = 1e-5;
tol2 = 1e-2;
C1 = zeros(m,n);
C2 = zeros(n,n);
C3 = zeros(n,n);
E  = zeros(m,n);
dist = zeros(n,n);
for iter = 1:max_iter
    if iter == 1
        Z = Z_ini;
        S = Z_ini;
        U = Z_ini;
        dist=L2_distance_1(H_ini,H_ini);
        %dist = zeros(n,n);
    end
    S_old = S;
    U_old = U;
    Z_old = Z;
    E_old = E;
   
    % -------- Update Z --------- %
    Z = Ctg*(X'*(X-E+C1/miu)+S+U-(C2+C3)/miu);
    Z = Z- diag(diag(Z));
    % -------- Update S --------- %
    S     = Z+(C2-dist)/miu;
    S     = S - diag(diag(S));
    for ic = 1:n
        idx    = 1:n;
        idx(ic) = [];
        S(ic,idx) = EProjSimplex_new(S(ic,idx));% 
    end
    
    % -------- Update U --------- %
    [AU,SU,VU] = svd(Z+C3/miu,'econ');
    AU(isnan(AU)) = 0;
    VU(isnan(VU)) = 0;
    SU(isnan(SU)) = 0;
    SU = diag(SU);    
    SVP = length(find(SU>lambda1/miu));
    if SVP >= 1
        SU = SU(1:SVP)-lambda1/miu;
    else
        SVP = 1;
        SU = 0;
    end
    U = AU(:,1:SVP)*diag(SU)*VU(:,1:SVP)';    
    % ------- Update E ---------- %
    temp1 = X-X*Z+C1/miu;
    temp2 = lambda2/miu;
    E = max(0,temp1-temp2) + min(0,temp1+temp2);   
%     temp1 = X-X*Z+C1/miu;
%     E = solve_l1l2(temp1,lambda2/miu);
    
    % -------- Update C1 C2 C3 miu -------- %
    L1 = X-X*Z-E;
    L2 = Z-S;
    L3 = Z-U;
    C1 = C1+miu*L1;
    C2 = C2+miu*L2;
    C3 = C3+miu*L3;
    LL1 = norm(Z-Z_old,'fro')/norm(Z,'fro');
    LL2 = norm(S-S_old,'fro')/norm(S,'fro');
    LL3 = norm(U-U_old,'fro')/norm(U,'fro');
    LL4 = norm(E-E_old,'fro')/norm(E,'fro');
    SLSL = max(max(max(LL1,LL2),LL3),LL4);
    if miu*SLSL < tol2
        miu = min(rho*miu,max_miu);
    end
    % --------- obj ---------- %
    leq1 = max(max(abs(L1(:))),max(abs(L2(:))));
    stopC = max(leq1,max(abs(L3(:))));
    obj_error(iter)=LL1+LL2+LL3;
    if stopC < tol
        %iter
        break;
    end   
end
%     t=1:iter;
%     figure
%     plot(t,obj_error,'r-'),xlabel('µü´ú´ÎÊý');ylabel('Îó²î');
end



function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end

