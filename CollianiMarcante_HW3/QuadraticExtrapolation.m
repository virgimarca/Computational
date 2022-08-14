function [xstar] = QuadraticExtrapolation(x1,x2,x3,x4) 
X = [x1 x2 x3 x4];
n = length(x1);
y = zeros(n,4); 
for j = 2:4 
    y(:,j)= X(:,j)-X(:,1); 
end 
Y = [y(:,2) y(:,3)]; 
[Q,R] = GramSchmidt(Y); % Using Gramh-Schmidt to get QR factorization
gamma = -R\(Q'*y(:,4)); % Using QR factorization to solve the system
beta0 = gamma(1)+gamma(2)+1; 
beta1 = gamma(2)+1; 
xstar = beta0*x2+beta1*x3+x4; 
end