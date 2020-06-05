function sol = gfun_alm(x)

% Inequality constraint
% Example 1
%sol = (0.25*x(1)*x(1) + 0.75*x(2)*x(2) -1);

% Example 2
%sol = x(1)-1; %(Part 2, 3rd problem)

% Example 3
sol = x(1)-x(2)*x(2); %(Part 2, 1st problem)