% Steepest decent with armijo backtracking. 
% para = [alpha, rho, c] 
function [xk,alpha] = SteepestDecentDirection(fun,grad,x0,para)
assert(para(1)>0); assert(para(2)>0); assert(para(2)<1); assert(para(3)<1);
assert(para(3)>0);
dk  = -grad(x0)/norm(grad(x0));

%armijo backtracking
alpha = para(1); rho = para(2); c=para(3);
while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk)
    alpha = rho*alpha;
    if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk
        break;
    end
end
xk = x0 + alpha*dk;


end