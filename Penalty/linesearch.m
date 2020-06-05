% This code gives the line search technique using Golden Section Method
function ReturnValue = linesearch(functname,tol,x,s,lower_bound,intr,trials)
format compact;

% finding the upper bound
upper_val = minfuncheck(functname,x,s,lower_bound,intr,trials);
u_b=upper_val(1);	

% if upper bound is close to lowbound reverse the direction of the search vector 
if (u_b <= 1.0e-06)
   l_b = lower_bound;  
   xL = x + l_b*s;  
   f_lb =feval(functname,xL);
   ReturnValue =[l_b f_lb x];
   return
end

if (tol == 0)
    tol = 0.0001;  
end

eps = tol/(u_b - lower_bound);
tau = 0.38197;
nmax = round(-2.078*log(eps)); 

l_b = lower_bound; 
a1 = (1-tau)*l_b + tau*u_b; 
x1 = x + a1*s; 
fa1 = feval(functname,x1);
a2 = tau*l_b + (1 - tau)*u_b;
x2 = x + a2*s; 
fa2 = feval(functname,x2);

for i = 1:nmax

	if fa1 >= fa2
   	l_b = a1;	
   	a1 = a2;
    fa1 = fa2;
    a2 = tau*l_b + (1 - tau)*u_b;
    x2 = x + a2*s; 
    fa2 = feval(functname,x2);
	else
  	   u_b = a2;
       a2 = a1;
       fa2 = fa1;
       a1 = (1-tau)*l_b + tau*u_b;	
       x1 = x + a1*s; 
       fa1 = feval(functname,x1);
	end
end
% returns the value
ReturnValue =[a1 fa1 x1];