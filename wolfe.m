function bool = wolfe(fun,grad,x0,dk,alpha,C,delta,sigma)
%WOLFE check if Non monotone wolfe conditions hold
% The output is boolean
if (fun(x0+alpha*dk) <= C + delta*alpha*grad(x0)'*dk && grad(x0 + alpha*dk)'*dk >= sigma*grad(x0)'*dk)
    bool = true;
else
    bool = false;
end
end

