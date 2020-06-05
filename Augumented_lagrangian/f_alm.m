function sol = f_alm(x)
% objective function

% Example 1
%sol = x(1)^4 - 2*x(1)*x(1)*x(2) +x(1)*x(1) + x(1)*x(2)*x(2) -2*x(1) + 4;

% Example 2
%sol = log(x(1))-x(2); %(Part 2, 3rd problem)

% Example 3
sol = abs(x(1)-2)+abs(x(2)-2); %(Part 2, 1st problem)