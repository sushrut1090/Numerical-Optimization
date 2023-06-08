clc;clear all;close all;
% new
% f = @(x) 100*( x(2)^2 - 2*x(2)*x(1) + x(1)^2) + (1 - 2*x(1) + x(1)^2);
% grad = @(x) [100*(-2*x(2)+2*x(1))-2+2*(x(1));
%              200*x(2)-200*x(1)];
% hess = @(x) [202 -200;
%             -200 200];
% x0 = [-1.8;-1.8];

% rosenbruck
f = @(x) 100*(x(2).^2-x(1).^2).^2 + (1-x(1)).^2;

grad = @(x)[-400*x(1)*(x(2).^2 -x(1).^2) - 2*(1 - x(1));
            400*x(2)*(x(2).^2 - x(1).^2)];
hess = @(x) -[400*(x(2).^2-3*x(1)^2)+2 800*x(1)*x(2);
            800*x(1)*x(2) -400*(3*x(2)^2 -x(1)^2)];
x0 = [-1.8;-1.8];
sol = [1;1];

% ex2
% f = @(x) x(1).^2 + 4*x(2).^2 + 2*x(1)*x(2);
% grad = @(x) [2*x(1)+2*x(2);
%               8*x(2) + 2*x(1)];
% hess = @(x) [2 2;
%             2 8];
% x0 = [-3;-3];
% sol = [0;0];

% ex3
% f = @(x) (x(1)+2*x(2)-7).^2 + (2*x(1)+x(2)-5).^2;
% grad = @(x) [10*x(1)+8*x(2)-34;
%             8*x(1)+10*x(2)-38];
% hess = @(x) [10 8;
%             8 10];
% x0 = [-9.5;9.5];

% ex4
% f = @(x) 5*x(1)^4 + 6*x(2)^4 - 6*x(1)^2 + 2*x(1)*x(2) + 5*x(2)^2 + 15*x(1) - 7*x(2) + 13;
% grad = @(x) [20*x(1)^3-12*x(1)+2*x(2)+15;
%              24*x(2)^3+2*x(1)+10*x(2)-7];
% hess = @(x) [60*x(1)^2-12 2;
%              2 72*x(2)^2+10];
% x0 = [1.9;-1.9];

% cube
%  f = @(x) 100*(x(2)-x(1).^3).^2 + (1-x(1)).^2;
%  grad = @(x) [100*(6*x(1)^5-6*x(2)*x(1)^2)+2*x(1)-2;
%               200*(x(2)-x(1)^3)];
%  hess = @(x)[100*(30*x(1)^4-12*x(1)*x(2))+2 -600*x(1)^2;
%              -600*x(1)^2 200];
%  x0 = [-1.2;-1];
% sol = [1;1];

% powell
% f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4 + 10*(x(1)-x(4))^4;
% grad = @(x) [2*(x(1)+10*x(2))+40*(x(1)-x(4))^3;
%              20*(x(1)+10*x(2))+4*(x(2)-2*x(3))^3;
%               10*(x(3)-x(4))-8*(x(2)-2*x(3))^3;
%               -5*(x(3)-x(4))-40*(x(1)-x(4))^3];
% hess = @(x) [2+120*(x(1)-x(4))^2, 20, 0, -120*(x(1)-x(4))^2;
%             20, 200+12*(x(2)-2*x(3))^2, -24*(x(2)-2*x(3))^2, 0;
%             0, -24*(x(2)-2*x(3))^2, 10+48*(x(2)-2*x(3))^2, -10;
%             -120*(x(1)-x(4))^2, 0, -5, 5+120*(x(1)-x(4))^2];
% x0 = [3;-1;0;1];
% sol = [0;0;0;0];

% wood
% f = @(x) 100*(x(1)^2-x(2))^2 + (x(1)-1)^2 + (x(3)-1)^2 + 90*(x(3)^2-x(4))^2;
% grad = @(x)[200*(x(1)^2-x(2))*(2*x(1))+2*(x(1)-1);
%             -200*(x(1)^2-x(2));
%             2*(x(3)-1)+90*(x(3)^2-x(4))*(2*x(3));
%             -180*(x(3)^2-x(4))];
% x0 = [1.5;1.5;1.5;1.5];
% sol = [1;1;1;1];

% trig n=2
% f = @(x) (1 + (1-cos(x)) - sin(x) - cos(x))^2;
% grad = @(x) (sin(x)-cos(x)+sin(x))*2*(1 + (1-cos(x)) - sin(x) - cos(x));
% x0 = 1/5;


para = [0.1,0.9,0.1,0.9,1.3,0.9];
para2 = [1,0.1,0.9];
nk = (para(1)+para(2))/2;
C = f(x0); Q = 1;
x = x0; 
B = eye(size(x0,1));
H=B;
for i = 1:500
%[xx,alpha(i),fun_count(i),grad_count(i)] = SteepestDescentDirection(f,grad,x,para2);
%[xx,alpha(i),fun_count(i),grad_count(i)] = NonMonotoneLineSearch(x,f,grad,para,C,hess); 
[B,H,xx,alpha(i)] = BFGS(f,grad,x,B);
%[xx,alpha(i),fun_count(i),grad_count(i)] = NewtonDirection(f,grad,hess,x,para2);
x = xx;

QQ = nk*Q +1;
CC = (nk*Q*C + f(xx))/QQ;
%fun_count(i) = fun_count(i)+1;
error(i) = norm(xx - sol);
if (norm(grad(xx)) < 1e-8 || error(i)<1e-4)
    break;
end
Q = QQ; C = CC;
end
figure(1)
plot(alpha)
title('Step size alpha over iterations')
xlabel('iteration')
ylabel('Step size')
figure(2)
plot(error)
title('Error ||analytic - final||')
xlabel('iteration')
ylabel('Error')