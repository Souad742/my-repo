s =100
rand("seed",s)

n=190;

Ad=rand(n,n);

A=triu(Ad);
xex = rand(n,1);

b = A*xex;

x = usolve(A,b);

fErrorB = norm(xex-x,2)/norm(xex,2)

bErrorB = norm(b-A*x,2)/norm(b,2)

