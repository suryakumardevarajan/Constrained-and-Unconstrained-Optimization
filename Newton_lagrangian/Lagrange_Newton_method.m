% This function uses Lagrange Newton method to calculate the optimal point
% of a function.
% Only Equality constraints are taken into account.

% The method to give the input is as follows. 
% x = sym('x',[2,1]); (example 1)
% a = 1;
% b = [1;2];
% c = [12,3;3,10];
% fx = a + (b'*x)+(0.5*x'*c*x);
% cx = -1+x(1)^2+x(2)^2;
% (example 2)
% a = [9, 1, 7, 5, 4, 7; 1, 11, 4, 2, 7, 5; 7, 4, 13, 5, 0, 7; 5, 2, 5, 17, 1, 9; 4, 7, 0, 1, 21, 15; 7, 5, 7, 9, 15, 27];
% b = [1, 4, 5, 4, 2, 1];
% c = 5;
% fx = x'*a*x + b*x+c;
% cx = -1+x(1)^2+x(2)^2;

% [xhat1,i1]= Lagrange_Newton_Eq(fx,cx,[1;0],1e-4)
% where fx is the function we want to optimize, cx is constraints, initial
% point and tolerance

function [xold,Iterations] = Lagrange_Newton_Eq(func,constraint,x_initial,eps)

n = length(x_initial);
m = length(constraint);
Iterations = 0; % Initialization 
x = sym('x',[n,1]); 
A = sym('A',[n,m]);
lambda = ones(m,1); 
xold = x_initial;
y = zeros(n,n);
for i = 1:m 
    % Calculation of A and lambda
    A(:,i) = gradient(constraint(i),x);
    y = y + lambda(i)*hessian(constraint(i),x);
end
% Calculating Wx
Wx = hessian(func,x) - y;
% It can be viewed to deliver a K-T point
% Expression for L
L = [Wx,-A;-A',zeros(m,m)]; 
% Y
Y = [-gradient(func,x);constraint]; 
% Stopping conditions
LxGrad = gradient(func,x) - sum(lambda.*gradient(constraint,x)); 
LyGrad = constraint;
while (((norm(subs(LxGrad,x,xold)) > eps) || (norm(subs(LyGrad,x,xold)) > eps)) && (Iterations<50)) % Start iteration 
    % Calculate deltax and lambda
    dxdl = inv(double(subs(L,x,xold)))*double(subs(Y,x,xold));
    % Update the variables
    delta = dxdl(1:n,1);
    lambdanew = dxdl(n+1:length(dxdl),1);
    xnew = xold + delta;
    xold = xnew;
    lambda = lambdanew;
    y = zeros(n,n);
    for i = 1:m
        A(:,i) = gradient(constraint(i),x);
        y = y + lambda(i)*hessian(constraint(i),x);
    end
    Wx = hessian(func,x) - y;
    L = [Wx,-A;-A',zeros(m,m)];
    Y = [-gradient(func,x);constraint];
    LxGrad = gradient(func,x) - sum(lambda.*gradient(constraint,x));
    Iterations = Iterations+1;
end
end
