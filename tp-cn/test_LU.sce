n=3
A=rand(n,n)
//[L,U]=mylu3b(A)
//[L,U]=mylu3b1boucle(A)
[L,U,p]=mylu(A)
erreur=A-L*U


