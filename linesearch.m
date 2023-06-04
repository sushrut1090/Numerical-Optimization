function xk = linesearch(fun,x0,para,method) 
assert(para(1)>0); assert(para(2)>0); assert(para(2)<1); assert(para(3)<1);
assert(para(3)>0);
a = fun;
a.x = x0;
dk = method(fun,x0);


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