% This code is for Augumented Lagrangian method invovling one inequality
% and one equality constraint. We use secant method for minimization
% We need {f_alm, gfun_alm, hfun_alm, func_alm, secant_DFP, linesearch,
% minfuncheck}.m files to run this program

clear all
clc
% Global values to be used in other functions
global lamda beta
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
lamda = 1; % default value used in equality constraint
m = 1; % number of inequality constraints
beta = 1; % default value used in inequality constraint
iter_alm = 20; % limit of iter used in this program
iter_secant = 20; % limit of iter used in secant method

%Following are the variables used in secant method
tol = 1.0e-08; 
lb = 0; 
int = 1.0; 
ntrials = 20;
epsilon = 1.0e-016; % error value

ch = 10; % scaling factor for lamda
cg = 10; % scaling factor for gamma

% calculating and storing the values
f_old = f_alm(x_old);

g_old = gfun_alm(x_old);
g_err = g_old*g_old';
gErr(1) = g_err;
rg = f_old/g_err; %Penalty factor

h_old = hfun_alm(x_old);
h_err = h_old*h_old';
hErr(1) = h_err;
rh = f_old/h_err; %Penalty factor

% store values   
iter = 1;  
x_s(iter,:) = x_old;
lamda_s(iter,:) =lamda; 
beta_s(iter,:)= beta; 
rh_s(iter) = rh;
rg_s(iter) = rg;
f_s(iter) = f_old;
uncon_f_s(iter)=func_alm(x_old);

% method starts here
x_cal = x_old;
count = 0;

for loop = 1:iter_alm			
 
   count = count + 1;
   x_initial = x_cal ;			
   
   functname = char('func_alm');
   fresult = secant_DFP(functname,x_initial,iter_secant,tol,lb,int,ntrials);
   
   for ix = 1:n
      xnew(ix) = fresult(ix);
   end
   Fnew = fresult(n+1);
   
   fend = f_alm(xnew);
   
   iter = iter + 1;
   x_s(iter,:) = xnew;
   f_s(iter) = fend;
   uncon_f_s(iter) = Fnew;
   
   gVend = gfun_alm(xnew);
   g_err = gVend*gVend';
   gErr(iter) = g_err;
   
   hVend = hfun_alm(xnew);
   h_err = hVend*hVend';
   hErr(iter) = h_err;
    
   x_initial = xnew;
      
   % convergence
   
   fdiff = Fnew - uncon_f_s(iter - 1);
   gdiff = gErr(iter) - gErr(iter -1);
   hdiff = hErr(iter) - hErr(iter -1);
   
   % stop if iterations are exceeded
   if loop == iter_alm
   fprintf('maximumum number of iterations reached : %6i \n',iter_alm);
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
   
   if abs(fdiff) < epsilon && abs(gdiff) <epsilon && abs(hdiff) < epsilon
       fprintf('Convergence in f : % 14.3E  reached in %6i iterations \n', ...
           abs(fdiff), loop);
       fprintf('\n The values for x :\n');
       disp(x_s);
       fprintf('\n Optimal value x :\n');
       disp(xnew);             
      break;
   else
       fprintf('\n ALM iteration: '),disp(count)
       fprintf('Value of (X) :  '),disp(x_initial);
   end
   
   x_cal = x_initial; 
   lamda = lamda + 2*rh*hVend;
   beta = beta + 2*rg*(gVend-beta/(2*rg));
   lamda_s(iter,:) =lamda;
   beta_s(iter,:)= beta;

   rh = rh*ch;
   rg = rg*cg;
   rh_s(iter) = rh;
   rg_s(iter) = rg;

end


