clear all   
close all  
clc   
  
format long e  
%EXERCISE 3  
  
% LOADING OF THE SOURCE OF DATA  
load('X.mat')   
   
% INITIALIZATION OF PARAMETERS  
sigma = 1;   
s = @(xi,xj) exp(-(sqrt(sum((xi-xj).^2,2))).^2/sigma);   
K = 13;   
%K = 40;  
n = length(X);   
  
% CREATION OF MATRIX W USING THE KNN ALGORITHM  
W = spalloc(n,n,(K+1)*n);  
for i= 1:n   
    xi = X(i,:) .* ones(900, 1);  
    distances_xi= s(xi, X);  
      
    [B, I] = sort(distances_xi, "descend");  
    B = B(2:K+1); %values  
    I = I(2:K+1); %indexes  
    for k = 1:K 
       W(i,I(k))= B(k);   
    end
end
  
%W has to be symmetric, so I add what it needs to be symmetric  
for i = 1:n  
    for j = 1:n  
        if W(i,j) ~= 0 && W(j,i) == 0  
            W(j,i) = W(i,j);  
        end  
    end  
end  
 
figure(1)  
scatter(X(:, 1), X(:, 2), '.')  
figure (2)  
spy(W)  
figure (3)  
plot(graph(W))  
%%  
% EXERCISE 4  
  
%CALCULATING THE MATRIX L_sym  
d = sum(W)';  
D = spdiags(d, 0, n, n);  
D_12 = spdiags(1./sqrt(d), 0, n, n);  
B = D_12*W*D_12;  
L_sym = spdiags(ones(n, 1), 0, n, n) - B;  
  
%CALCULATING THE 5 SMALLEST EIGENVALUES OF L_sym  
eigenvalues = eigs(L_sym,5,'smallestreal');  

%%  
% EXERCISE 5   
% To apply the power deflation method we need a matrix that is SPD so that   
% the algorithm is stable, for this reason we create the matrix B_mod   
mu = max(sum(B, 2));   
B_mod = B + spdiags(mu*ones(n, 1), 0, n, n);   

% WE APPLY THE POWER DEFLATION METHOD TO B_mod  
tol = 1e-16;   
nmax = 100000;   
x0 = ones(n, 1);   
lambda = zeros(5,1);   
v = zeros(n, 5);  
[lambda(1), v_1, iter] = power_method(B_mod, tol, nmax, x0);  
v(:, 1) = v_1;  
tol = 1e-10;
[lambda(2:5)] = power_deflation_cycles(B_mod,v_1,tol,nmax, 4);  
  
for i = 1:5  
    x0 = rand(n, 1);   
    [lambda(i), v(:, i), iter] = inverse_power_shift(B_mod, lambda(i), tol, nmax, x0);   
end  

lambda = -lambda+mu+1;  
[lambda, indexes] = sort(lambda);  
v = v(:,indexes);    
%%  
% EXERCISE 6  
% Relative error  
norm(eigenvalues(1:5)-lambda(1:5))/norm(eigenvalues(1:5)) 
  
[eigenvalues, lambda] 
figure(4) 
hold on  
scatter(zeros(5,1),lambda, 'filled')  
scatter(zeros(5,1),eigenvalues,200,'x', 'Linewidth', 1) 
 
legend('lambda','eienvalues') 
hold off
m=3;
%m=2;
%% 
%EXERCISE 7 
[V,~] = eigs(L_sym,5,'smallestabs'); 
for i = 1: m
    V(:,i) = abs(V(:,i)/norm(V(:,i))); 
end

for i = 1:5  
    v(:,i) = abs(v(:,i)/norm(v(:,i)));  
end  

% Relative error between eigenvectors caluclated by Matlab and by us 
for i = 1: m
    norm(V(:,i)-v(:,i))/norm(V(:,i))
end
%% 
%EXERCISE 8 
U = zeros(900,m); 
for i = 1 : m
    U(:,i) = v(:,i); 
end

for i = 1:900 
    U(i,:) = abs(U(i,:)/norm(U(i,:))); 
    norm(U(i,:)) 
end 
%% 
%EXERCISE 9 
idu = kmeans(U,m); 
%% 
%EXERCISE 10 
figure(5)
gscatter(X(:,1),X(:,2),idu)
%%
%EXERCISE 11
[idx, ~, D] = spectralcluster(W, m,'Distance', 'precomputed','LaplacianNormalization','symmetric');
figure(6)
gscatter(X(:,1),X(:,2),idx)
spectral_eigenvalues = D(1:m, :);
spectral_eigenvalues = diag(diag(spectral_eigenvalues));
norm(spectral_eigenvalues- lambda(1:m))/norm(spectral_eigenvalues)
%%
%EXERCISE 12
idxx = kmeans(X, m);
figure(7)
gscatter(X(:,1),X(:,2),idxx)