clear all 
close all 
clc 
%% Getting the data
fileID = fopen('california.txt','r'); 
data = fscanf(fileID,'%d'); 
column= data(2:2:end); 
row= data(1:2:end); 
 
row=row+1; 
column=column+1; 
n = 9664; % Number of nodes

%% Creating the necessary matrix for the computation
W = sparse(row,column,1,n,n); %adjacency matrix 
v = 1/n*ones(n, 1); %giving to each node the same probability
deg = sum(W,2); % Outdegree of each node
d = zeros(n, 1);
for i=1:n
    if deg(i)==0 
        d(i, 1) = 1; %giving value 1 to the node with outdegree=0
    end
end 
% Normalizing matrix W to P
deg_1 = 1./deg; 
D_1 = spdiags(deg_1, 0, n, n);  
P = D_1*W; 

% Work with P_prime because P has row equal to all zero --> no sense
D = d*v';
P_prime = P + D;

%% General parameter for each method
beta = 0.15;  %value of constant beta 
%beta = 0.10;
%beta = 0.05;
%beta = 0.01;
mu = 1/n*ones(n,1);  %unitary input 
 
norma = 1000;  %
tol = 1e-3; % Tolerance
nmax = 1000; % Number of max iteration

Pcap = ((1-beta)*P_prime'+beta/n*ones(n,n)); 

%% Iterative method 
z_old_b = 1/n*ones(n,1);  %initialization of centrality vectors 
z_new_b = 1/n*ones(n,1); 

iterit=0;
tic
while norma > tol 
    z_new_b = (1-beta)*P_prime'*z_old_b + beta*mu;  %update
    norma = norm(z_new_b - z_old_b,1);  % check error 
    z_old_b = z_new_b; 
    iterit = iterit+1;
end 
z_new_b = z_new_b/norm(z_new_b); % normalization of the pagerank vector
Teiterative = toc; 

%% Quadratic Power Method 
% Choosing how each iteration to do the quadratic extrapolation
for p=1:30 
    [xq, iterq(p), errquad] =                          Quadratic_power_method(Pcap, tol ,nmax, v, p);
end
[itermin, ind] = min(iterq)

%% Comparison between Power method and Quadratic Power method 
tic
[xp, iterp, errpower] = power_method_markov(Pcap , tol ,nmax, v);
Tepower = toc

tic
[xq, iterq, errquad] = Quadratic_power_method(Pcap, tol ,nmax, v, 10);
Tequad = toc
%% Plot of the convergence for the Power method and the Quadratic Power method
hold on
figure(1)
plot((1:1:iterq),errquad)
plot((1:1:iterq),errpower(1:iterq))
legend('Quadratic Power', 'Power')
xlabel('Number of iterations') 
ylabel('L1 residual')
hold off