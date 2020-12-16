function[sol,niter,info]= myjacobi(A,b,nmaxit,tol)
    //vérification si aucun terme de la diagonal de A n'est nulle
    if ~and(diag(A)) then
        error('erreur:diagonale est nulle')
    end
    //décomposition de A: A=D-E-F
     D=diag(diag(A))
     E=-triu(A)+D
     F=-tril(A)+D
     x=inv(A)*b
   
     sol=b
     niter=0 
     info=0 
     err=[]     
           for k=1:nmaxit
               sol =(eye(n,n)-inv(D)*A)*sol+inv(D)*b
               err=[err,norm(x-sol)];
                if max(abs(A*sol-b))< tol 
                    info = 1;
                    niter= k;
                    break
                 end
            end
 xtitle('le graphe de convergence pour la méthode de jacobi')
 plot(1:niter,log(err),xtitle)    
endfunction
n=3
A=[2 -1 0;-1 2 -1;0 -1 2]
b=[1; 2; 3]
[sol,niter,info]= myjacobi(A,b,100,0.01)
x=inv(A)*b
b=A*x;
x=inv(A)*b
M=eye(n,n)-(inv(diag(diag(A)))*A); 
R_spectral=max(abs(spec(M)));

disp(R_spectral)




