function [ lambda, x , iter ]= inverse_power_shift(A , mu , tol , nmax , x0 )

%inverse_power_shift(A,mu,tol,nmax,x0) performs the inverse power method
%with shift

%lambda = inverse_power_shift(A) computes the eigenvalue of A of minimum
%modulus with the inverse power method using the default error tolerance 
%(1e-6) the default maximum number of iterations (100) and the default initial
%vector.

%lambda = inverse_power_shift(A,mu) computes the eigenvalue of A closest to
%the given number (real or complex) mu  using the default error tolerance 
%(1e-6) the default maximum number of iterations (100) and the default initial
%vector.

%lambda = inverse_power_shift(A,mu,tol,nmax,x0) uses the absolute error
%tolerance tol, the maximum number of iterations nmax and the initial vector
%x0 provided by the user.

%[lambda, v, iter] = inverse_power_shift(A,mu,tol,nmax,x0) also returns the
%eigenvector v such that A*v = lambda*v and the iteration number at which v
%was computed.

[n,m]= size(A);
if n ~= m
    error ('Only for square matrices'); 
end

if nargin == 1
    x0 = rand (n ,1);
    nmax = 100;
    tol = 1e-06;
    mu = 0;
elseif nargin == 2
    x0 = rand(n,1);
    nmax = 100;
    tol = 1e-06;
end
[L,U]= lu (A - mu*eye(n));
% if norm ( x0 ) == 0
% x0 = rand (n ,1);
% end
x0 = x0 / norm ( x0 );
z0 = L \ x0 ;
pro = U \ z0 ;
lambda = x0'*pro ;
err = tol*abs(lambda)+1;
iter =0;

while (err > tol*abs(lambda)) && (iter <= nmax)
    x = pro;
    x = x/norm(x);
    z=L\x;
    pro = U\z;
    lambdanew = x'*pro;
    err = abs(lambdanew-lambda);
    lambda = lambdanew ;
    iter = iter + 1;
end

lambda = 1/lambda + mu ;

end