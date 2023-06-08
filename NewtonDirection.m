% Steepest decent with armijo backtracking. 
% para = [alpha, rho, c] 
function [xk,alpha,funcounter,gradcounter] = NewtonDirection(fun,grad,hess,x0,para)
assert(para(1)>0); assert(para(2)>0); assert(para(2)<1); assert(para(3)<1);
assert(para(3)>0);
funcounter = 0;
gradcounter = 1;
dk  = -hess(x0)\grad(x0);

%wolfe 
% alpha = para(1); rho = para(2); c=para(3);
% c2 = 0.95;
% assert(para(3)<c2);
% while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk && ...
%         -dk'*grad(x0+alpha*dk) <= -c2*dk'*grad(x0))
%     alpha = rho*alpha;
%     if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk  && ...
%         -dk'*grad(x0+alpha*dk) <= -c2*dk'*grad(x0)
%         break;
%     end
%     funcounter = funcounter+4;
%     gradcounter = gradcounter+6;
% end

%armijo backtracking
alpha = para(1); rho = para(2); c=para(3);
while (fun(x0 + alpha*dk) > fun(x0) + c*alpha*grad(x0)'*dk)
    alpha = rho*alpha;
    if fun(x0 + alpha*dk) <= fun(x0) + c*alpha*grad(x0)'*dk
        break;
    end
    funcounter = funcounter +4;
    gradcounter = gradcounter+2;
end
xk = x0 + alpha*dk;


end