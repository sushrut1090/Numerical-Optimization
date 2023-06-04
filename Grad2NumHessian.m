function NumHess = Grad2NumHessian(fun,x0,h)
%UNTITLED Gets analytic function and outputs numerical hessian for a given
%point x0. h is the step size of FD.

sz = size(x0);
NumHess = eye(sz(1));
xlen = sz(1);
dx = zeros(xlen,1);

end