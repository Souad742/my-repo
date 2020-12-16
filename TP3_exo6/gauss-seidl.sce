//méthode de jaccobi pour résoudre Ax=b

function[sol,niter,info]= gauss_seidl(A,b,nmaxit,tol)
    //vérification aucun terme de la diagonal de A n'est nul
    if ~and(diag(A)) then
        error('erreur:diagonale est nulle')
    end
    //décomposition de A=D-E-F
     D=diag(diag(A))
     E=-triu(A)+D
     F=-tril(A)+D
     x=inv(A)*b
   
     sol=b
     niter=0 
     info=0 
     err=[]     
           for k=1:nmaxit
               sol =inv(D-E)*((F*sol)+b)
               err=[err,norm(x-sol)];
                if max(abs(A*sol-b))< tol 
                    info = 1;
                    niter= k;
                    break
                 end
            end
 xtitle('le graphe de convergence pour la méthode de gauss-seidel')
plot(1:niter,log(err),xtitle)    
endfunction

n=3;
A=[2 -1 0;-1 2 -1;0 -1 2]
//D=diag(diag(A))
//E=-triu(A)+D
//F=-tril(A)+D
//A=D-E-F
b=[1; 2; 3]
[sol,niter,info]= gauss_seidl(A,b,100,0.01)

x=inv(A)*b
M=eye(n,n)-inv(D-E)*A; 
R_spectral=max(abs(spec(M)));
//disp(norm(M))
 //disp(M)
 //R_spectral=(cos(%pi/n+1))^2
 disp(R_spectral)




