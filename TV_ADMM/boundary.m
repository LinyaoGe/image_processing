function u = boundary(u,Gx,Gy)
   [M,N] = size(u);
   u(:,1) = u(:,2) - Gy(:,1);
   u(:,N) = u(:,N-1) - Gy(:,N);
   u(1,:) = u(2,:) - Gx(1,:);
   u(M,:) = u(M-1,:) - Gx(M,:);
end

