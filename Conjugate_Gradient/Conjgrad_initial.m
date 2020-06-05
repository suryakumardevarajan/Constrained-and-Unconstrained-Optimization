clc 
clear all 
close all

% Function definition
syms x1 x2;
%f=-(sqrt((x1^2 +1)*(2*x2^2 + 1)))/(x1^2+x2^2+0.5);
f = (x1^2 + x2^2);
% f = abs(x1-2) + abs(x2-2);
%f = x1 - x2 + 2*x1^2 + 2*x1*x2 + x2^2;



% Initial values
x(1) = 4;
y(1) = 4;
err = 1e-2; 
i = 1; 
alpha = 0.01;


% Gradient Computation
df_dx = diff(f, x1);
df_dy = diff(f, x2);
grad = [subs(df_dx,[x1,x2], [x(1),y(1)]) subs(df_dy, [x1,x2], [x(1),y(1)])]; 
grad = grad';
s_j = -grad;



while(norm(s_j)>err)
    %s_j = -grad;
    % finding the next value
    x(i+1) = x(i)+alpha*s_j(1); 
    y(i+1) = y(i)+alpha*s_j(2); 
    grad_new=[subs(df_dx,[x1,x2], [x(i+1),y(i+1)]) subs(df_dy, [x1,x2], [x(i+1),y(i+1)])];
    grad_new= grad_new';
    beta = ((grad_new- grad)' * grad_new)./ norm(grad)^2;
    s_j = -grad_new + beta * s_j;
    grad = grad_new;
    i = i+1;
end

% Plots
fcontour(f, 'Fill', 'On');
hold on;
plot(x,y,'*-r');

% Displaying results
fprintf('Number of Iterations for the cost function convergence: %d\n\n', i);
fprintf('Mininimum point: [%d,%d]\n\n', x(i), y(i));

