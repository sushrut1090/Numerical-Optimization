% Steepest decent with armijo backtracking. 
% para = [alpha, rho, c] 
function [xk,alpha,funcounter,gradcounter] = SteepestDecentDirection(fun,grad,x0,para)
assert(para(1)>0); assert(para(2)>0); assert(para(2)<1); assert(para(3)<1);
assert(para(3)>0);
dk  = -grad(x0)/norm(grad(x0));
funcounter = 0;
gradcounter = 2;
%armijo backtracking
alpha = para(1); rho = para(2); c=para(3);
while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk)
    alpha = rho*alpha;
    if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk
        funcounter = funcounter+2;
        gradcounter = gradcounter+1;
        break;
    end
    funcounter = funcounter+4;
    gradcounter = gradcounter+2;
end
xk = x0 + alpha*dk;
end