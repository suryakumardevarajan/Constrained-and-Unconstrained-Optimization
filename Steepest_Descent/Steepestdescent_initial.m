% This is the code of steepest descent algorithm
clc
clear
format long

% Function definition
syms x1 x2;
f=(sqrt((x1^2 +1)*(2*x2^2 + 1)))/(x1^2+x2^2+0.5); 
%f = x1 - x2 + 2*x1^2 + 2*x1*x2 + x2^2;

% Initial values
x(1) = 4;
y(1) = 4;
err = 1e-6; 
i = 1; 

% Gradient Computation
df_dx = diff(f, x1);
df_dy = diff(f, x2);
grad = [subs(df_dx,[x1,x2], [x(1),y(1)]) subs(df_dy, [x1,x2], [x(1),y(1)])]; 

% Assigning the steepest descent search direction
s_j= -(grad);

% Algorithm of steepest descent
while (norm(s_j)>err)    
    del_V_j = [x(i),y(i)]';  
    
    % Assigning alpha
    syms alpha; % Step size
    g = subs(f, [x1,x2], [x(i)+s_j(1)*alpha,y(i)+alpha*s_j(2)]);
    dg_dalpha = diff(g,alpha);
    alpha = solve(dg_dalpha, alpha);    
    % finding the next value
    x(i+1) = del_V_j(1)+alpha*s_j(1); 
    y(i+1) = del_V_j(2)+alpha*s_j(2); 
    i = i+1;
    % Updating Gradient
    grad = [subs(df_dx,[x1,x2], [x(i),y(i)]) subs(df_dy, [x1,x2], [x(i),y(i)])];
    % New Search Direction
    s_j = -(grad); 
end

% Plots
fcontour(f, 'Fill', 'On');
hold on;
plot(x,y,'*-r');

% Displaying results
if (norm(s_j)<err)
    % when x_j satisfies the necessary condition for local minimum   
    disp("Minimum point is obtained");
end
fprintf('Number of Iterations for the cost function convergence: %d\n\n', i);
fprintf('Mininimum point: [%d,%d]\n\n', x(i), y(i));


