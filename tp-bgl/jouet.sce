function jouet_variable(n, e)
    format ("e",e);
    A=rand(n,n);
    xex = rand(n, 1);
    b = A*xex;
    x = A\b;
    frelres= norm(x-xex)/norm (xex);
    brelres= norm(b-A*x)/norm (b);
    
    disp(frelres);
    disp(brelres);
    cap = cond(A);
    disp(cap);
    
    endfunction
    
    

    
 
