
clc;
clear;
image = imread('lena.jpg');
subplot(2,2,1);imshow(image);
f = imnoise(image, 'gaussian',0,0.01);
subplot(2,2,2);imshow(f);

%------------------初始化参数-------------------
f = double(f);
u = f;
[m,n] = size(f);

%bita的值
bita1 = zeros(m,n,2);
bita2 = zeros(m,n);
bita3 = zeros(m,n);
bita4 = zeros(m,n,2);

% mu的值,lamda的值
mu1 = 5;
mu2 = 0.00001;
mu3 = 100;
mu4 = 10;
lamda = 10000;

%w的值
w = zeros(m,n,2);

%q的值
q = zeros(m,n);

%n 和 m的值
N = zeros(m,n,2);
M = zeros(m,n,2); 

%------------------迭代-------------------
total = 200;
for step = 1:total
    
    %k+1步的u
    div_w = divergence(w(:,:,1),w(:,:,2));
    div_bita1 = divergence(bita1(:,:,1),bita1(:,:,2));
    for i = 1:2
        Au = center_diff(u);
        u = (f - mu1*div_w + mu1*Au - div_bita1) ./ (1+4*mu1);
    end
    
    %k+1步的w
    [ux,uy] = gradient(u);
    gra_u = cat(3,ux,uy);
    A = gra_u + ((mu2 + bita2).* M - bita1)./mu1 ;
    abs_A = sqrt(A(:,:,1).^2 + A(:,:,2).^2);
    w(:,:,1) = max(abs_A - (mu2 + bita2)./mu1,0) .* A(:,:,1) ./ (abs_A + eps);
    w(:,:,2) = max(abs_A - (mu2 + bita2)./mu1,0) .* A(:,:,2) ./ (abs_A + eps);

    %k+1步的q
    div_n = divergence(N(:,:,1),N(:,:,2));
    q = max(abs(div_n - bita3./mu3) - lamda./mu3,0) .* (div_n - bita3./mu3) ./ abs(div_n - bita3./mu3+ eps) ;
    
    %k+1步的m
    M = N + ((mu2 + bita2).*w + bita4)./mu4;
    for i =1:2
        abs_m = sqrt(M(:,:,1).^2 + M(:,:,2).^2);
        M(:,:,1) = M(:,:,1) ./ (max(abs_m,1) + eps);
        M(:,:,2) = M(:,:,2) ./ (max(abs_m,1) + eps);
    end
    %k+1步的n
    [qx,qy] = gradient(q);
    [bita3x,bita3y] = gradient(bita3);
    gra_q = cat(3,qx,qy);
    gra_b3 = cat(3,bita3x,bita3y);
    h = -gra_q - gra_b3./mu3 + mu4./mu3.*M - bita4./mu3;
    for i = 1:2
        [nx,ny] = n_diff(N(:,:,1),N(:,:,2));
        An = cat(3,nx,ny);
        N = mu3 .* (An + h)./(mu4+4*mu3*mu4);
    end
    
    %k+1步的bita
    [ux,uy] = gradient(u);
    gra_u = cat(3,ux,uy);
    bita1 = bita1 + mu1*(w - gra_u);
    bita2 = bita2 + mu2*(sqrt(w(:,:,1).^2 + w(:,:,2).^2) - N(:,:,1).*w(:,:,1)-N(:,:,2).*w(:,:,2));
    bita3 = bita3 + mu3*(q - N(:,:,1) - N(:,:,2));
    bita4 = bita4 + mu4*(M - N);
    
%   判断能量函数
    Ax = ux ./ sqrt(ux.^2 + uy.^2);
    Ay = uy ./ sqrt(ux.^2 + uy.^2);
    E(step) = sum(sum(lamda.*abs(divergence(Ax,Ay))+0.5.*((u-f).^2)));
%     if step >2
%         result(step) = abs((E(step) - E(step-1) )/ E(step));
%         if abs((E(step) - E(step-1) )/ E(step)) < 0.0000001
%             break;
%         end
%     end
end
subplot(2,2,3);imshow(uint8(u));title('result');
subplot(2,2,4);plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
PSNR = psnr(uint8(u),uint8(image))
PSNR2 = psnr(uint8(f),uint8(image))







