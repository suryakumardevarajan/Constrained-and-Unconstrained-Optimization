% This code gives the minimum of a function 
function sol = minfuncheck(functname,x,s,initial_step,inc_step,ns)

format compact
if (ns ~= 0) 
    ntrials = ns;
else
    ntrials = 10;  
end

if (inc_step ~= 0) 
    step = inc_step;
else
    step = 1;  
end

% finds a value of function greater than or equal to the previous value
for i = 1:ntrials
   j = 0;	
   del_a = j*step;	
   a_old = initial_step + del_a;  
   dx0 = a_old*s;	
   x0 = x + dx0;  
   f0 = feval(functname,x0);
   j = j+1;	
   del_a = j*step;	
   a_new = initial_step + del_a;
   dx1 = a_new*s;	
   x1 = x + dx1;	
   f1 = feval(functname,x1);
   f1s = f1;
   if f1 < f0 
         for j = 2:ntrials
         	a_new = initial_step + j*step;	
            dx1 = a_new*s;	
            x1 = x + dx1;	
            f1 = feval(functname,x1);
            f1s = min(f1s,f1);
            if f1 > f1s 
      			sol = [a_new f1 x1];
               return;
            end
         end
         fprintf('\nFunction cannot be increased in ntrials')
      	sol = [a_new f1 x1];
         return;	         
   else	f1 >= f0;
      step = 0.5*step;
   end
end
sol =[initial_step f0 x0];
   
   