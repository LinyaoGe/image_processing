function final_image = get_image(u)
    final_image = zeros(256,256);
    for i = 1:256
        for j = 1:256
            final_image(i,j) = u(i+1,j+1);
        end
    end
end

