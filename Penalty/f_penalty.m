function sol = f_penalty(x)
% objective function
sol = x(1)^4 - 2*x(1)*x(1)*x(2) +x(1)*x(1) + x(1)*x(2)*x(2) -2*x(1) + 4;
%sol = log(x(1))-x(2); %(Part 2, 3rd problem)
%sol = abs(x(1)-2)+abs(x(2)-2); %(Part 2, 1st problem)