function [matrice]=tridiag(A)//fonction pour rendre une matrice tridiagonale
    n=size(A,1);
  for i=1:n
     matrice(i,i)=A(i,i);
    end
    for (i=2:n)
            matrice(i-1,i)=A(i-1,i);
            matrice(i,i-1)=A(i,i-1);
    end



endfunction

function [L,U]=fact(matrice)//effectuer la factorisation LU 
    [matrice]=tridiag(A);
  for k=1:n-1
i=k+1:n;
matrice(i,k)=matrice(i,k)/matrice(k,k);
j=k+1:n;
matrice(i,j)=matrice(i,j)-matrice(i,k)*matrice(k,j)
end

U=triu(matrice);//matrice tridiagonale supérieure

L=tril(matrice)//matrice tridiagonale inférieure
endfunction

n=4
A=rand(n,n)
[matrice]=tridiag(A)
[L,U]=fact(matrice)

