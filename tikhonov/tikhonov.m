%读取图像
Image = imread('lena.jpg');

%给图像加噪声
f = imnoise(Image,'gauss');
subplot(1,2,1)
imshow(f);
hold on;
f = double(f);

%加边界
% f2 = border(f);

%设置参数
u=f;
h = 0.5;
lamda = 0.8;

%遍历
for step = 1:200
    f2 = border(u);
    u = f2;
    [M,N] = size(u);
    for i = 2:M-1
        for j = 2:N-1
            u(i,j) = (f2(i,j) + lamda/h^2 * (u(i-1,j) + u(i,j-1) + u(i+1,j) + u(i,j+1))) / (1 + 4*lamda/h^2);
        end
    end
    u = get_image(u);
    [Fx,Fy] = gradient(u);
    E(step) = sum(0.5 * sum((u-f).^2) - 0.5*lamda*sum(Fx.^2 + Fy.^2));
    if step > 2
        if abs((E(step)-E(step-1))/E(step)) < 0.01
            break;
        end
    end
end

%取图
% final_image = zeros(256,256);
% for i = 1:256
%     for j = 1:256
%         final_image(i,j) = u(i+1,j+1);
%     end
% end
subplot(1,2,2)
imshow(uint8(u));







