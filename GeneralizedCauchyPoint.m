function [xc] =  GeneralizedCauchyPoint(x,l,u,grad,theta,W,M)
n = size(x,1);
t = zeros(n,1);
d = zeros(n,1);
g = grad(x); % grad is provided analytically hence we compute g here
for i = 1:n
    % initilize t
    if(g(i)<0)
        t(i) = (x(i)-u(i))/g(i);
    elseif (g(i)>0)
        t(i) = (x(i)-l(i))/g(i);
    else
        t(i) = 1e16;
    end
    
    % initilize d
    if (abs(t(i)) < 1e-15)
        d(i) = 0;
    else
        d(i) = -g(i);
    end
end
 % set F 
 F = [];
 for i = 1:n
     if (t(i)>0)
         F = [F;i];
     end
 end
 p = W'*d;
 c = 0;
df = -d'*d; % f' computation
d2f = -theta*df' - p'*M*p; % f"
dt_min = -df/d2f;
t_old = 0;
T  = min(t); % assimuing all i are in F -  code it to ensure
% flag = flase;
% while flag = false
%     a =  find(min(t));
b = find(T,t);
F = F(F~=b);
dt = t-0;

while (dt_min>dt)

end


end