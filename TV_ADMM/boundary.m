function u = boundary(u,Gx,Gy)
   [M,N] = size(u);
   for i = 2:257
       u(i,1) = u(i,2) - Gy(i-1,1);
       u(i,N) = u(i,N-1) - Gy(i-1,N-2);
       u(1,i) = u(2,i) - Gx(1,i-1);
       u(M,i) = u(M-1,i) - Gx(M-2,i-1);
   end
end

