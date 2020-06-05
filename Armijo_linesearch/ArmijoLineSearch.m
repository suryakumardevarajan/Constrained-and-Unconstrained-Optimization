% This code is for Armijo line search which is implemeted with Steepest descent algorithm
% We need {fun_armijo and grad_armijo}.m files to run this program

function [argmin, iterations] = ArmijoLineSearch(fun_armijo,x0,sigma,beta, eps)
if (nargin == 2)
    sigma = 1e-4;
    beta = 0.5;
    eps = 1e-4;
elseif (beta>=1 || beta<=0 || eps <=0)
    disp('Wrong paramaters used beta should be Beta and sigma in (0,1), eps >0');
end

x = x0;
k = 0;
maxIterations = 1e8;
g = grad_armijo(x);

while norm(g) > eps
    g = grad_armijo(x);
    d = -g;
    l = 0; 
    t = 0;
    while(fun_armijo(x+(beta^l)*d)>fun_armijo(x)+(beta^l)*sigma*g'*d)
        l = l+1;
    end
    t = beta^l;
    k = k+1;
    x = x+t*d;    
    if k == maxIterations
        break;
    end
end

disp("Minimum point is obtained");
fprintf('Number of Iterations for the cost function convergence: %d\n\n', k);
fprintf('Mininimum point: \n');
disp(x(1))
disp(x(2))

