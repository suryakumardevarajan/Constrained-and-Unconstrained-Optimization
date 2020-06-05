
% This function calculates the unconstrained fuction for Penalty method

function sol = func_penalty(x)

global m leq
global rh rg
sol = f_penalty(x);

if leq > 0
   hval = hfun_penalty(x);
   sol = sol +  rh*(hval*hval);
end

if m > 0
   gval = gfun_penalty(x);
   sol = sol + rg*(max(0,gval))^2;
end
