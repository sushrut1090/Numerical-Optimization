function [Bk,xk,alpha] = BFGS(fun,grad,x0,B)
dk = -B\grad(x0);

para = [1,0.1,0.99];
%armijo backtracking
alpha = para(1); rho = para(2); c=para(3);
while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk)
    alpha = rho*alpha;
    if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk
        break;
    end
end
epsilon = 1e-19; % to remove 0/0 case when B is updated.
sk = alpha*dk ;
xk = x0 + sk;
yk = grad(xk) - grad(x0);
Bk = B + (yk*yk' + epsilon)/(yk'*sk + epsilon) - (B*sk)*(sk'*B')/(sk'*B*sk);

end