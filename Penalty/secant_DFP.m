% This is the code for secant algorithm with Davidon Fletcher Powell(DFP) Method
% Requires:  minfuncheck.m, linsearch.m, gradfunction.m and function file
% Parameters:  initial value of x, number of iterations, tolerance, the initial value of stepsize,
% incremental value and number of trials
% eg: Conjugategradient('fun',[0 0],30, 0.0001, 0,1 ,15)

function sol = surya_secant_ALM(functname,x_initial,iter,tol,lb,int,ntrials)

clf 
% convergence
err = 1.0e-04;  
n = length(x_initial);% number of variables


% Algorithm
x_i(1,:) = x_initial;
x = x_initial;
line_col = 'r';
opt_f(1) = feval(functname,x); 
opt_x(1)=0;
grad = (gradfunction(functname,x)); 
H = eye(n);  % Initial H is identity matrix function
conv(1)=grad*grad';

for i = 1:iter-1
   % determine search direction
   s = (-H*grad')'; 
   x_o = linesearch(functname,tol,x,s,lb,int,ntrials);
   opt_x(i+1) = x_o(1);
   opt_f(i+1) = x_o(2);
   for k = 1:n
      x_i(i+1,k)=x_o(2+k);
      x(k)=x_o(2+k);
   end
   grad= (gradfunction(functname,x));
   conv(i+1)=grad*grad';
   
  
if(conv(i+1)<= err) 
    break; 
end 
del_x = (x - x_i(i,:))';
del_g = (grad -gradfunction(functname,x_i(i,:)))'; 
H = H + ((del_x*del_x')/(del_x'*del_g)) - (((H*del_g)*(H*del_g)')/(del_g'*H*del_g));%DFP
end

len=length(opt_x);
opt_value=x_i(length(opt_x),:);

%Displaying results
%fprintf('\nNo. of iterations:  '),disp(len) % Comment for ALM
%fprintf('\nOptimal point, Value of func at minimal pt \nduring the
%iterations\n')  % Comment for ALM
%disp([x_i opt_f'])  % Comment for ALM
sol = [opt_value opt_f(len)];