function [x1 , iter, err]= power_method_markov(A , tol , nmax , x0 )
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
err(1) = 1;
iter = 1;

while (err(iter) > tol) && (iter <= nmax) % Check on norm and max iteration
    x0 = x0/norm(x0);
    x1 = A*x0;
    iter = iter + 1;
    err(iter) = norm(x0 -x1, 1);
    x0 = x1;
end
end