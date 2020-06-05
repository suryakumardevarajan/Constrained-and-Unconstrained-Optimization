function sol = hfun_penalty(x)

% the equality constraint vector

sol = [(x(1)*x(1) + x(2)*x(2) -2)]; 

%sol = [(x(1)*x(1) + x(2)*x(2) -4)]; %(Part 2, 3rd problem)

%sol = x(1)*x(1) + x(2)*x(2) -1; %(Part 2, 1st problem)