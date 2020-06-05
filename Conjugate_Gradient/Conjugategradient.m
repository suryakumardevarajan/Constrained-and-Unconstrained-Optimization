% Conjugate gradient for n variables
% Requires:  minfuncheck.m, linsearch.m, gradfunction.m and function file
% Parameters:  initial value of x, number of iterations, tolerance, the initial value of stepsize,
% incremental value and number of trials
% eg: Conjugategradient('fun',[0 0],30, 0.0001, 0,1 ,15)

function sol = Conjugategradient(functname,x_initial,iter,tol,lb,int,ntrials)
clf 
format short e   
% convergence
err = 1.0e-04;  
n = length(x_initial); %  number of variables
if (n == 2)

% Contour plot
	x1 = -5:0.1:5;
	x2 = -5:0.1:5;
	len_x1 = length(x1);
	len_x2 = length(x2);
	for i = 1:len_x1
        for j = 1:len_x2
            x1x2 =[x1(i) x2(j)];
            fun(j,i) = feval(functname,x1x2);
        end
    end
    c1 = contour(x1,x2,fun,[5 10 50 100 200 500 1000],'k');
     % [3.1 3.25 3.5 4 6 10 15 20 25],'k');
	grid
	xlabel('x1')
    ylabel('x2')
	title('Conjugate Gradient Method')
end

%Algorithm Starts
x_i(1,:) = x_initial;
x = x_initial;
line_col = 'r';
opt_f(1) = feval(functname,x); % value of function at start
opt_x(1)=0;
s_j = -(gradfunction(functname,x)); 
conv(1)=s_j(1)*s_j(1)';

for i = 1:iter-1
   % determine search direction
   if (i > 1)
      beta = grad*grad'/(s_j*s_j');
      s_j = -grad + beta*s_j;
   end
   x_o = linesearch(functname,tol,x, ...
      s_j,lb,int,ntrials);
   opt_x(i+1) = x_o(1);
   opt_f(i+1) = x_o(2);
   for k = 1:n
      x_i(i+1,k)=x_o(2+k);
      x(k)=x_o(2+k);
   end
   grad= (gradfunction(functname,x));
   conv(i+1)=grad*grad';
   
   % Line draw
   if (n == 2)
      line([x_i(i,1) x_i(i+1,1)],[x_i(i,2) x_i(i+1,2)],'LineWidth',2,'Color',line_col)
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
opt_value=x_i(length(opt_x),:);

% Displaying Results
fprintf('\nNo. of iterations:  '),disp(len)
fprintf('\nOptimal point, Value of func at minimal pt \nduring the iterations\n')
disp([x_i opt_f'])
sol = [opt_value opt_f(len)];