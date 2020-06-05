% Box’s evolutionary optimization is a simple optimization technique developed by G. E. P. Box1
% in 1957. The algorithm requires (2N +1) points, of which 2N
% are corner points of an N-dimensional hypercube2
% centred on the other point.
% All (2N + 1) function values are compared and the best point is identified.
% In the next iteration, another hypercube is formed around this best point. If
% at any iteration, an improved point is not found, the size of the hypercube is
% reduced. This process continues until the hypercube becomes very small.

% Referred from optimization for engineering design by kalyanmoy deb

function [xnew] = runevolution(x,R)

del = [2,2];          
x_initial =x;
tol = 10^-3;
y = x_initial;
del1 = (sqrt(sum(del.^2)));
j=1;

while del1 > tol
    
    z1 = x_initial - del./2;
    z2 = x_initial + ([-(del(1)),del(2)])./2;
    z3 = x_initial + ([del(1),-(del(2))])./2;
    z4 = x_initial + del./2;    
    x = [x_initial;z1;z2;z3;z4]';
    f = zeros (1,5);
    for i=  1:5
        f(i) = func_pb(x(:,i),R);
    end
    [~ , idx] = min(f(:,:));
    y = x(:,idx);
    i=1;                        
    if (y == x_initial')
        del = del./2;
    else
        x_initial = y;
        x_initial =x_initial';
    end
    del1 = (sqrt(sum(del.^2)));
    j = j+1;
end
xnew = y;
