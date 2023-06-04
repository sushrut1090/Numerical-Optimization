function dk = SteepestDescentDk(grad,x0)
dk  = -grad(x0)/norm(grad(x0));
end