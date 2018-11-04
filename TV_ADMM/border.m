function border_image = border(f)
    [M,N] = size(f);
    border_image = zeros(M+2, N+2);
    for i = 2:M+1
        for j = 2:N+1
            border_image(i,j) = f(i-1,j-1);
        end
        if i == 2
            border_image(1,:) = border_image(2,:);
        elseif i == M+1
            border_image(M+2,:) = border_image(M+1,:); 
        end
    end
    border_image(:,1) = border_image(:,2);
    border_image(:,N+2) = border_image(:,N+1);
end

