% function used in penalty and barrier function method
function P = func_pb(x,R)
     x1 = x(1);
     x2 = x(2);
     
     %f = R.*(((x1-5).^2 +x2.^2 -26).^2);       % R * constraint function
     %P = (x1.^2 +x2 -11).^2 + (x1 +x2.^2-7).^2 +f;  % Objective function

     %f = R.*(-x1*x2);       % R * constraint function
     %P = (-x1+x2^2+1).^2 + (x1 +x2).^2 +f;  % Objective function
     
     f =  R.*(x1-x2*x2);
     P =  abs(x1-2)+abs(x2-2)+f; 
end