function jouet(m)
    
   mat=rand(m,m)
  
   //disp (mat)
   xex=rand(m,1)
   b=mat*xex;
   x=mat\b;
   frelres=norm(x-xex)/norm(xex)
   brelres=norm(b-mat*x)/norm(b)
   disp("frelres=",frelres);
   disp("brelres=",brelres);
endfunction
