function sol = func_alm(x)
% 
% calculates the unconstrained fuction for ALM method
%
global lamda beta
global m leq
global rh rg
sol = f_alm(x);

if leq > 0
   hval = hfun_alm(x);
   sol = sol + lamda*hval'+ rh*(hval*hval');
end

if m > 0
   gval = gfun_alm(x);
   for j = 1:m
      g(j) = max(gval(j),-beta(j)/(2*rg));
   end
   sol = sol + beta*g'+ rg*(g*g');
end
