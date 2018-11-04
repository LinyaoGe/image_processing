clc;clear;close all;
%读取图像
Image = imread('lena.jpg');

%给图像加噪声
f = imnoise(Image,'gauss');
subplot(1,2,1)
imshow(f);
hold on;
f = double(f);

%设置参数
[M,N] = size(f);
lamda = 35;
h = 0.5;
s = 256*256;arf = 0.00001;
mu = s*arf;
u = zeros(M,N);
[Wx,Wy] = gradient(u);
[Bx,By] = gradient(u);

for step = 1:10
    %第一步
    div_w = divergence(Wx,Wy);
    div_bita = divergence(Bx,By);
    F = f - (mu*div_w + div_bita);
    u = f;
    %处理边界
    u = border(u);
    Gx = mu*Wx + Bx;
    Gy = mu*Wy + Bx;
    u = boundary(u,Gx,Gy);
    %图像内部
    for i = 2:257
        for j = 2:257
             u(i,j) = (F(i-1,j-1) + mu * (u(i-1,j) + u(i,j-1) + u(i+1,j) + u(i,j+1))) / (1 + 4*mu);
        end
    end
    
    %取图
    u = get_image(u);
    
    %第二步迭代W
    [ux,uy] = gradient(u);
    abs_u = sqrt(ux.^2 + uy.^2);
    Ax = ux - Bx./mu;
    Ay = uy - By./mu;
    abs_A = sqrt(Ax.^2 + Ay.^2);
    Wx = max(abs_A - lamda./mu,0).*Ax./(abs_A+eps);
    Wy = max(abs_A - lamda./mu,0).*Ay./(abs_A+eps);
    %第三部update bita
    Bx = Bx + mu*(Wx - ux);
    By = By + mu*(Wy - uy);
    %判断终止条件
    E(step)=sum(sum(lamda.*abs_u+0.5.*(u-f).^2));
    if (step>2)
        if abs((E(step)-E(step-1))/E(step))<0.000000001
            fprintf('step:%d',step);
            break;
        end
    end
end

subplot(1,2,2);
imshow(uint8(u));






















