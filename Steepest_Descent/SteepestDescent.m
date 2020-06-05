% Steepest Descent Method for n variables
% Requires:  minfuncheck.m, linesearch.m, gradfunction.m and function file
% Parameters:  initial value of x, number of iterations, tolerance, the initial value of stepsize,
% incremental value and number of trials
% eg: SteepestDescent('func',[0 0],30, 0.0001, 0,1 ,15)

function sol = SteepestDescent(functname,x_initial,iter,tol,lb,int,ntrials) 
clf 
err = 1.0e-08;   
n = length(x_initial); % number of variables 
if (n == 2)
%  plotting contours
	x1 = -5:0.1:5;
	x2 = -5:0.1:5;
	len_x1 = length(x1);
	len_x2 = length(x2);
	for i = 1:len_x1
   	for j = 1:len_x2
      	x1_x2 =[x1(i) x2(j)];
      	fun(j,i) = feval(functname,x1_x2);
   	end
    end
	c1 = contour(x1,x2,fun,[3.1 3.25 3.5 4 6 10 15 20 25],'k');
	grid
	xlabel('x1')
	ylabel('x2')
    title('Steepest Descent method');
end
%  Algorithm starts
x_i(1,:) = x_initial;
x = x_initial;
line_col = 'r';
opt_f(1) = feval(functname,x); % value of function at start
opt_x(1)=0;
s_j = -(gradfunction(functname,x)); % steepest descent

conv(1)=s_j*s_j';
for i = 1:iter-1
   
   x_o = linesearch(functname,tol,x, s_j,lb,int,ntrials);
   opt_x(i+1) = x_o(1);
   opt_f(i+1) = x_o(2);
   for k = 1:n
      x_i(i+1,k)=x_o(2+k);
      x(k)=x_o(2+k);
   end
   s_j = -(gradfunction(functname,x)); % steepest descent
   conv(i+1)=s_j*s_j';
   
   % lines in contour
   if (n == 2)
      line([x_i(i,1) x_i(i+1,1)],[x_i(i,2) x_i(i+1,2)],'LineWidth',2, 'Color',line_col)
   if strcmp(line_col,'r') 
         line_col = 'k';
   else
         line_col = 'r';
   end    
    pause(1)  
   end
   
   if(conv(i+1)<= err) 
       break; 
   end 
   
end
len=length(opt_x);
opt_pt=x_i(length(opt_x),:);

% Displaying Results
fprintf('\nNo. of iterations:  '),disp(len)
fprintf('\n - Optimal point , Value of function at Optimal point \nduring the iterations\n')
disp([x_i opt_f'])
sol = [opt_pt opt_f(len)];