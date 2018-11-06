function u = center_diff(u)
    [M,N] = size(u);
    left = circshift(u, [0,-1]);left(:,M) = left(:,M-1);
    right = circshift(u, [0,1]); right(:,1) = right(:,2);
    up = circshift(u, [1,0]);
    up(N,:) = up(N-1,:);
    down = circshift(u, [-1,0]);
    down(1,:) = down(2,:);
    
    u = left + right + up + down;
end

