function [xk,alpha] = NonMonotoneLineSearch(x0,fun,grad,para,C,hess)
%NONMONOTONELINESEARCH implementation of NLSE with wolfe's condition 
% para =[eta_min,eta_max,delta,sigma,rho,mu]
% Description in Hager & Zhang 
%dk = direction(grad,x0);
if nargin == 6
    dk = -hess(x0)\grad(x0);
else
    dk  = -grad(x0)/norm(grad(x0));
end
alpha =1;
n_min = para(1); n_max = para(2); delta = para(3); sigma = para(4); rho=para(5);
mu = para(6);
kk = 1;
%wolfe's condition
while (wolfe(fun,grad,x0,dk,alpha,C,delta,sigma) == false)
    alpha = alpha/2;
    kk = kk+1;
    if (kk == 100000)
       error("Wolfe condition not met in 100000 iterations") 
    end
end

xk = x0 + alpha*dk;
end

