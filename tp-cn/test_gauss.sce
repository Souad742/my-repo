n=100;
A=rand(n,n);

xex = rand(n,1);

b = A*xex;

x = gausssolve(A,b);

fErrorB = norm(xex-x)/norm(xex)

bErrorB = norm(b-A*x)/norm(b)


