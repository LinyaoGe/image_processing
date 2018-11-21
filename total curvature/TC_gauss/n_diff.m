function [nx,ny] = n_diff(Nx,Ny)
    [M,N] = size(Nx);
    [n11,n12] = gradient(Nx);
    [n21,n22] = gradient(Ny);
    
    %处理n12
    A12_up = circshift(n12,[-1,0]);
    A12_up(M,:) = A12_up(M-1,:);
    r12 = A12_up - n12;
    R12 = circshift(r12,[0,1]);
    R12(:,1) = R12(:,2);
    
    %处理n22
    A22_left = circshift(n22,[0,-1]); 
    A22_left(:,N) = A22_left(:,N-1);
    A22_right = circshift(n22,[0,1]);
    A22_right(:,1) = A22_right(:,2);
    
    %处理n11
    A11_up = circshift(n11,[-1,0]);
    A11_up(M,:) = A11_up(M-1,:);
    A11_down = circshift(n11,[1,0]);
    A11_down(1,:) = A11_down(2,:);
    
    %处理n21
    A21_left = circshift(n21,[0,-1]);
    A21_left(:,M) = A21_left(:,M-1);
    r21 = A21_left - n21;
    R21 = circshift(r21,[-1,0]);
    R21(M,:) = R21(M-1,:);
    
    nx = A11_up + A11_down - 2*n11 + A21_left - n21 - R21;
    ny = A22_left + A22_right - 2*n22 + A12_up - n12 - R12;
    
    
    
end