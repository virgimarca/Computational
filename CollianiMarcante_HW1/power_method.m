function [ lambda ,x , iter ]= power_method(A , tol , nmax , x0 )

%power_method(A,tol,nmax,x0) computes the eigenvalue with maximum modulus
%of a real matrix .

%lambda = power_method(A) computes with the power method the eigenvalue of
%A of maximum modulus from an initial guess which by defaul is an all one
%vector using a default absolute error tolerance TOL (1e-6) and a default
%number of iterations nmax (100).

%lambda = power_method(A, tol, nmax, x0) uses an absolute error tolerance
%tol, a maximum number of iterations nmax and an initial vector x0 given as
%input by the user.

%[lambda, v, iter] = power_method(A,tol, nmax, x0) also returns the
%eigenvector V such that A*v = lambda*v and the iteration number at which v
%was compouted.

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
x0 = x0/norm(x0);
pro = A*x0;
lambda = x0'*pro ;
err = tol*abs(lambda)+1;
iter = 0;

while (err > tol*abs ( lambda )) && (iter <= nmax)
    x = pro;
    x = x/norm(x);
    pro = A*x;
    lambdanew = x'*pro ;
    err = abs(lambdanew-lambda);
    lambda = lambdanew ;
    iter = iter + 1;
end


end

