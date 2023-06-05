function [Bk,Hk,xk,alpha] = BFGS(fun,grad,x0,B,H)
if nargin == 5
    dk = -H*grad(x0);
else
    dk = -B\grad(x0);
end

para = [1,0.1,0.9];
c2 = 0.95;
assert(para(3)<c2);
%armijo backtracking
alpha = para(1); rho = para(2); c=para(3);
while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk && ...
        -dk'*grad(x0+alpha*dk) <= -c2*dk'*grad(x0))
    alpha = rho*alpha;
    if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk  && ...
        -dk'*grad(x0+alpha*dk) <= -c2*dk'*grad(x0)
        break;
    end
end
epsilon = 1e-10; % to remove 0/0 case when B is updated.
sk = alpha*dk ;
xk = x0 + sk;
yk = grad(xk) - grad(x0);
rho = (1+epsilon)/(sk'*yk + epsilon);
Bk = B + (yk*yk' + epsilon)/(yk'*sk + epsilon) - (B*sk)*(sk'*B')/(sk'*B*sk);
if nargin ==5
Hk = (eye(size(H)) - sk*yk'*rho)*H*(eye(size(H)) - sk*yk'*rho) + (rho + epsilon)*sk*sk'; 
end
Hk= Bk;
end