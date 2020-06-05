function [fun] = fun_armijo(x)
%fun = 2*x(1)^4+3*x(2)^4+2*x(1)^2+4*x(2)^2+x(1)*x(2)-3*x(1)-2*x(2);
fun = x(1)^2 + x(2)^2;
%fun = abs(x(1)-2)+abs(x(2)-2);
    % @return the gradient 

end
