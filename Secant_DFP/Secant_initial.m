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
x_old = 4;
y_old = 4;
err = 1e-4; 
i = 1; 
alpha = 0.1;


% Gradient Computation
df_dx = diff(f, x1);
df_dy = diff(f, x2);
grad = [subs(df_dx,[x1,x2], [x_old,y_old]) subs(df_dy, [x1,x2], [x_old,y_old])]; 
grad = grad';
H=eye(2);
s_j = -H* grad;




%loop
while(norm(s_j)>err)
    s_j = -H* grad;
    g = subs(f, [x1,x2], [x_old+s_j(1)*alpha,y_old+alpha*s_j(2)]);
    h = subs(f, [x1,x2], [x_old,y_old]);
    %if (g<h)
    % finding the next value
    x_new = x_old+alpha*s_j(1);
    y_new = y_old+alpha*s_j(2);
   % end
    del_x = x_new - x_old;
    del_y = y_new - y_old;
    del_xy = [del_x, del_x]';
    %del_xy =  alpha* s_j;
    grad_new=[subs(df_dx,[x1,x2], [x_new,y_new]) subs(df_dy, [x1,x2], [x_new,y_new])];
    grad_new = grad_new';
    del_g = grad_new - grad;
    H = H + ((del_xy* del_xy')./(del_xy' * del_g)) - (((H* del_g) * (H * del_g)')./(del_g' * H * del_g));
    x_old = x_new;
    y_old= y_new;
    grad = grad_new;
    i=i+1;
end



% Plots
fcontour(f, 'Fill', 'On');
hold on;
plot(x_old,y_old,'*-r');

% Displaying results
%if (norm(s_j)<err)
    % when x_j satisfies the necessary condition for local minimum   
 %   disp("Minimum point is obtained");
%end
fprintf('Number of Iterations for the cost function convergence: %d\n\n', i);
fprintf('Mininimum point: [%d,%d]\n\n', x_old, y_old);
