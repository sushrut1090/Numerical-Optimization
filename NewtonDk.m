function dk = NewtonDk(grad,hess,x0)
dk  = -hess(x0)\grad(x0);
end