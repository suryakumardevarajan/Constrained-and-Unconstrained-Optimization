% This code gives the method for penalty and barrier functions
% Reference text book "OPTIMIZATION FOR ENGINEERING DESIGN" by Kalyanboy
% This code needs func_pb.m and runevolution.m files

eps1 = 10e-5;       
eps2 = 10e-5;
rg_initial = 0.1;           % initial Penalty parameter
x_initial =[rand*2 rand*2];
% x0 = [0 0];
c = input('Enter 0 for Barrier method, 1 for penalty Method = ');
 if c==0
     c = rand;              
 else
     c = rand*1.9;
 end
 t = 1;
 x(t,:)= x_initial;
 rg(t)= rg_initial;
 x1 = x(1);
 x2 = x(2);
 constraint = ((x1-5).^2 +x2.^2 -26); % change the constraint here for every example
 if (constraint < 0)
     for t = 1
         x1 = x(t,1);
         x2 = x(t,2);
         constraint = (x1-x2*x2);
         k = rg.*constraint^2;
         F = func_pb(x(t,:),0);
         P(1,t) = func_pb(x(t,:),rg(1,t));
         G(t,:) = runevolution(x(t,:),rg(1,t));     % to perform unconstrained search for given R
         x(t+1,:) = G(t,:);
         t = t+1;
         P(1,t) = func_pb(G(t-1,:),rg(1,t-1));
         err = abs(P(1,t));
         rg(1,t) = c*rg(1,t-1);
         x1 = x(t,1);
         x2 = x(t,2);
         constraint = ((x1-5).^2 +x2.^2 -26);
         
     end
 else
     x1 = x(t,1);
     x2 = x(t,2);
     k = 0;
     F = func_pb(x(t,:),0);
     P(1,t) = func_pb(x(t,:),0);
     G(t,:) = runevolution(x(t,:),0);
     t = t+1;
     P(1,t) = func_pb(G(t-1,:),0);
     err = abs(P(1,t));
     rg(1,t) = c*rg(1,t-1);
     x(t,:) = G(t-1,:);
 end
 
 a(:,1) = G(t-1,:)';

 while err>eps2
     G(t,:) = runevolution(x(t,:),rg(1,t));
     x(t+1,:) = G(t,:);
     t = t+1;
     rg(1,t) = c*rg(1,t-1);
     P(1,t) = func_pb(G(t-1,:),rg(1,t-1));
     err = abs(P(1,t)-P(1,t-1));
       x1 = x(t,1);
       x2 = x(t,2);
       constraint = ((x1-5).^2 +x2.^2 -26);
       
     a(:,1) = G(t-1,:)';
     if abs(constraint) <eps1
         break;
     elseif (abs(x(t,1)-x(t-1,1))<10e-5 && abs(x(t,2)-x(t-1,2))<10e-5)
         break;
     end
     
 end
fval = func_pb(x(t,:),0);
fprintf('The minimum point is (%4f,%4f) with a function value of %.4f  \n',a,fval)