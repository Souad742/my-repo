//méthode de jaccobi pour résoudre Ax=b

function[sol,niter,info]= myjacobi(A,b,nmaxit,tol)
    //vérification aucun terme de la diagonal de A n'est nul
    if ~and(diag(A)) then
        error('erreur:diagonale est nulle')
    end
    //décomposition de A=D-E-F
     D=diag(diag(A))
     E=-triu(A)+D
     F=-tril(A)+D
    
     
   
     sol=b
     niter=0 // aucune itération
     info=0 // pas de convergence
     //initialisation
     //boucle itérative de résolution 
  
     
           for k=1:nmaxit
               //sol =inv(D)*((E+F) *sol+b)
                 sol =(eye(n,n)-inv(D)*A)*sol+inv(D)*b
                  
                if max(abs(A*sol-b))< tol
                    info = 1
                    niter= k
                    break
                 end
            end
          
           
endfunction
n=3;
A=[2 -1 0;-1 2 -1;0 -1 2]

//D=diag(diag(A))
//E=-triu(A)+D
//F=-tril(A)+D
//A=D-E-F
b=[4; 3; 1]
[sol,niter,info]= myjacobi(A,b,100,0.01)
 M=eye(n,n)-inv(D)*A; 
 R_spectral=max(abs(spec(M)));
 disp(M)
  R_spectral=(cos(%pi/n+1))
 disp(R_spectral)
x=A\b


