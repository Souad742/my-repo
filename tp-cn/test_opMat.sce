n=3
 A=rand(n,n);
 B=rand(n,n);
t=zeros(1,20);
t1=zeros(1,20);
t2=zeros(1,20);
t3=zeros(1,20);
step=1;

for i=1:10;
    n=n+10;
    A=rand(n,n);
    B=rand(n,n);
for j=1:10;
tic();
C3= matmat3b(A,B);
t(step)=t(step)+toc();

tic();
C2= matmat2b(A,B);
t1(step)=t1(step)+toc();

tic;
C1= matmat1b(A,B);
t2(step)=t2(step)+toc();
tic()
A*B;
t3(step)=t3(step)+toc();
end
t(step)=t(step)/10;
t1(step)=t1(step)/10;
t2(step)=t2(step)/10;
t3(step)=t3(step)/10;
step=step+1;
end
t
t1
t2
t3
plot(t)
plot(t1)
plot(t2)
plot(t3)



