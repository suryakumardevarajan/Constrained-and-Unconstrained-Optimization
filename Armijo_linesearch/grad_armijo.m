function [grad] = grad_armijo(x)
syms x1 x2
%f = 2*x1^4+3*x2^4+2*x1^2+4*x2^2+x1*x2-3*x1-2*x2;
%f = x1^2 + x2^2;
v =[x1,x2];
df_dx = diff(fun_armijo(v), x1);
df_dy = diff(fun_armijo(v), x2);
grad = [subs(df_dx,[x1,x2], [x(1),x(1)]) subs(df_dy, [x1,x2], [x(1),x(2)])]; 
grad = grad';
end
