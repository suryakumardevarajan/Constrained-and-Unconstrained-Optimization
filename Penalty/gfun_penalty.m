function sol = gfun_penalty(x)

% Inequality constraint

sol = [(0.25*x(1)*x(1) + 0.75*x(2)*x(2) -1)];

%sol = x(1)-1; %(Part 2, 3rd problem)

%sol = x(1)-x(2)*x(2); %(Part 2, 1st problem)