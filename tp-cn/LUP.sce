
function [L,U] = mylu3b(A)
    
[n,n]=size(A);
for k = 1 : n-1
for i = k + 1 : n
A(i,k) = A(i, k)/A(k, k);
end 
for i = k + 1 : n
for j = k + 1 : n
   A(i,j)= A(i,j)- A(i,k)*A(k,j);
end 
end 
end  
  U=triu(A);
// expression de L
L=tril(A)
L(1:n+1:$)=1;
endfunction


function [L,U]=mylu3b1boucle(A)
[n,n]=size(A);
// verification A is square
for k=1:n-1
// calcul des multiplicateurs
// ceux-ci sont mis dans la partie triangulaire inferieure de A
i=k+1:n;
A(i,k)=A(i,k)/A(k,k);
// mise `a jour de la partie superieure de A
j=k+1:n;
A(i,j)=A(i,j)-A(i,k)*A(k,j)
end
//expression de U
U=triu(A);
// expression de L
L=tril(A)
// on met la diagonale de L `a 1
L(1:n+1:$)=1;
endfunction





function [L,U,p]=mylu(A)
n=size(A,1);
q=zeros(1,n);
row=[1:n];

for k=1:n-1
[piv,ind]=max( abs(A(k:n,k)));

ind = k-1+ind;
q(1,k)=row(1,ind);

 if (ind~=k)
    new = A(ind,:);
    A(ind,:)=A(k,:);
    A(k,:)=new;
    row(1,ind)=row(1,k); 
    row(1,k)=q(1,k); 
end
for i =k+1:n
  A(i,k)=A(i,k)/A(k,k);
end
 for i=k+1:n;
    for j=k+1:n;
      A(i,j)=A(i,j)-A(i,k)*A(k,j)
    end
 end
end
Isp=speye(n,n);
p = Isp(row,:);
L=tril(A,-1)+speye(n,n);
U=triu(A);

endfunction
