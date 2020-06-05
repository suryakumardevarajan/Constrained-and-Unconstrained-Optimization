% This code is for Penalty method invovling one inequality
% and one equality constraint. We use secant method for minimization
% We need {f_penalty, gfun_penalty, hfun_penalty, func_penalty, secant_DFP, linesearch,
% minfuncheck}.m files to run this program

clear all
clc

% Global values to be used in other functions
global n m leq
global rh rg

% Initialising the values
format compact
format short e
warning off  
x_old = [3 2]; % Initial value for the problem
Xmin = [-5 -5]; % less than the initial value (~ 3 times)
Xmax = [5 5]; % more than the  initial value (~ 3 times)
n = length(x_old);
leq = 1; % number of equality constraints
m = 1; % number of inequality constraints
iter_penalty = 20; % limit of iter used in this program
iter_secant = 20; % limit of iter used in secant method

%Following are the variables used in secant method
tol = 1.0e-08; 
lb = 0; 
int = 1.0; 
ntrials = 20;
epsilon = 1.0e-016; % error value

rh = 1;
rg = 1;
ch = 10; % scaling factor for lamda
cg = 10; % scaling factor for gamma

% calculating and storing the values
f_old = f_penalty(x_old);
g_old = gfun_penalty(x_old);
h_old = hfun_penalty(x_old);

% store values   
iter = 1;  
x_s(iter,:) = x_old;
rh_s(iter) = rh;
rg_s(iter) = rg;
f_s(iter) = f_old;
uncon_f_s(iter)=func_penalty(x_old);

% method starts here
x_cal = x_old;
count = 0;

for loop = 1:iter_penalty			
 
   count = count + 1;
   x_initial = x_cal ;			
   functname = char('func_penalty');
   fresult = secant_DFP(functname,x_initial,iter_secant,tol,lb,int,ntrials);
   for ix = 1:n
      xnew(ix) = fresult(ix);
   end
   Fnew = fresult(n+1);
   fend = f_penalty(xnew);
   iter = iter + 1;
   x_s(iter,:) = xnew;
   f_s(iter) = fend;
   uncon_f_s(iter) = Fnew; 
   gVend = gfun_penalty(xnew);
   hVend = hfun_penalty(xnew);  
   x_initial = xnew;
      
   % convergence
   
   fdiff = Fnew - uncon_f_s(iter - 1);   
   % stop if iterations are exceeded
   if loop == iter_penalty
   fprintf('maximumum number of iterations reached : %6i \n',iter_penalty);
   fprintf('\n The values for x  : \n');
   disp(x_s);
  break;
   end
   
   % convergence in changes in x      
   xdiff = (xnew-x_cal)*(xnew-x_cal)';
   if xdiff < epsilon
       fprintf('\n The values for x are : \n');
       disp(x_s);
       fprintf('\n Optimal value x :\n');
       disp(xnew);  
  break;
   end
   
   if ((fdiff)^2 < epsilon) 
       fprintf('Convergence in f : % 14.3E  reached in %6i iterations \n', ...
           abs(fdiff), loop);
       fprintf('\n The values for x :\n');
       disp(x_s);
       fprintf('\n Optimal value x :\n');
       disp(xnew);             
      break;
   else
       fprintf('\n Penalty iteration: '),disp(count)
       fprintf('Value of (X) :  '),disp(x_initial);
   end
   
   x_cal = x_initial; 
   rh = rh*ch;
   rg = rg*cg;
   rh_s(iter) = rh;
   rg_s(iter) = rg;
end


