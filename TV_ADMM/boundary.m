function u = boundary(u,Gx,Gy)
   [M,N] = size(u);
   for i = 2:257
       u(i,1) = u(i,2) - Gy(i-1,2);
       u(i,N) = u(i,N-1) - Gy(i-1,N-3);
       u(1,i) = u(2,i) - Gx(2,i-1);
       u(M,i) = u(M-1,i) - Gx(M-3,i-1);
   end
end

