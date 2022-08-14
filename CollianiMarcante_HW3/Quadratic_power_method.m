function [x4, iter, err]= Quadratic_power_method(A , tol , nmax , x0, p) 
 
[n,m] = size(A); 
if n ~= m 
   error ('Only for square matrices '); 
end 
 
if nargin == 1 
    tol = 1e-06; 
    x0 = ones(n,1); 
    nmax = 100; 
end 
%Initialization 
iter = 1; 
x1 = zeros(n, 1); 
x2 = zeros(n, 1); 
x3 = x0/norm(x0); 
err(1) = 1; 
 
while (err(iter) > tol) && (iter <= nmax)  % Check on norm and max iteration
    x4 = A*x3;  
    if mod(iter, p)==0 && iter>4 % Calling Quadratic extrapolation every p iteration
        x4 = QuadraticExtrapolation(x1, x2, x3, x4); 
    end
    x1 = x2; 
    x2 = x3; 
    x4 = x4/norm(x4);
    iter = iter +1;
    err(iter) = norm(x4-x3, 1);
    x3 = x4;
end 
x4 = x4/norm(x4);
end