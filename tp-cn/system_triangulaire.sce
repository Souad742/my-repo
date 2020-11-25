function [x] = lsolve(A,b)
    n=size(A,1);
    m=size(b,2);
    
x = zeros(size(b,1),size(b,2));
x(1,1:m)=(b(1,1:m))./A(1,1);
for i=1:n
    x(i,1:m)=(b(i,1:m)-A(i,1:(i-1))*(x(1:(i-1),1:m)))./A(i,i);
end        
    
endfunction


function [x] = usolve(A,b)
    n=size(A,1);
    m=size(b,2);
    
x = zeros(size(b,1),size(b,2));
x(n,1:m)=(b(n,1:m))./A(n,n);
for i=n-1:-1:1
    x(i,1:m)=(b(i,1:m)-A(i,i+1:n))*(x((i+1):n,1:m))./A(i,i);
end        
    
endfunction
