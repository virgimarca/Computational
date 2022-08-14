function [Q,R]=GramSchmidt(V) 
 
[m,n]=size(V); 
W=V; 
Q=zeros(m,n); 
R=zeros(n); 
 
for j=1:n 
    R(1:j-1,j)=Q(:,1:j-1)'*W(:,j); %when j=1 skip 
    W(:,j)=W(:,j)-Q(:,1:j-1)* R(1:j-1,j); 
    R(j,j)=norm(W(:,j),2); 
    Q(:,j)=W(:,j)/R(j,j); 
end