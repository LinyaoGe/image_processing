clc;
clear;
image = imread('lena.jpg');
subplot(2,2,1);imshow(image);title("原图");
f = imnoise(image, 'salt & pepper',0.1);
subplot(2,2,2);imshow(f);title("salt&pepper图");

%------------------初始化参数-------------------
[m,n] = size(f);
f = double(f);
u = f;

%bita的值
bita1 = zeros(m,n,2);
bita2 = zeros(m,n);
bita3 = zeros(m,n);
bita4 = zeros(m,n,2);
bita5 = zeros(m,n);

% mu的值,lamda的值
mu1 = 0.3;
mu2 = 0.01;
mu3 = 10;
mu4 = 10;
mu5 = 1;
lamda = 5;

%w的值
w = zeros(m,n,2);

%q的值
Q = zeros(m,n);

%n 和 m的值
N = zeros(m,n,2);
M = zeros(m,n,2);  

%v的值
v = zeros(m,n);

%------------------迭代-------------------
total = 100;
for step = 1:total
    
    %k+1步的u
    div_w = divergence(w(:,:,1),w(:,:,2));
    div_bita1 = divergence(bita1(:,:,1),bita1(:,:,2));
    for i = 1:2
        Au = center_diff(u);
        u = (v + f - mu1/mu5*div_w + mu1/mu5*Au - div_bita1./mu5 + bita5./mu5) ./ (mu5+mu5*4*mu1);
    end
    
    %k+1步的w
    [ux,uy] = gradient(u);
    gra_u = cat(3,ux,uy);
    A = gra_u + ((mu2 + bita2).*M - bita1)./mu1;
    abs_A = sqrt(A(:,:,1).^2 + A(:,:,2).^2);
    for i = 1:2
        w(:,:,1) = max(abs_A - (mu2 + bita2)./mu1,0) .* A(:,:,1) ./ (abs_A + eps);
        w(:,:,2) = max(abs_A - (mu2 + bita2)./mu1,0) .* A(:,:,2) ./ (abs_A + eps);
    end
    
    %k+1步的q
    div_n = divergence(N(:,:,1),N(:,:,2));
    Q = max(abs(div_n - bita3./mu3) - 1./mu3,0) .* (div_n - bita3./mu3) ./ abs(div_n - bita3./mu3 + eps) ;
    
    %k+1步的m
    M = N + ((mu2 + bita2).*w + bita4)./mu4;
    abs_m = sqrt(M(:,:,1).^2 + M(:,:,2).^2);
    M(:,:,1) = M(:,:,1) ./ (max(abs_m,1) + eps);
    M(:,:,2) = M(:,:,2) ./ (max(abs_m,1) + eps);
    
    %k+1步的n
    [qx,qy] = gradient(Q);
    Aq = cat(3,qx,qy);
    [bita3x,bita3y] = gradient(bita3);
    Ab3 = cat(3,bita3x,bita3y);
    H = -Aq - Ab3./mu3 + mu4./mu3.*M - bita4./mu3;
    for i = 1:2
        [nx,ny] = n_diff(N(:,:,1),N(:,:,2));
        An = cat(3,nx,ny);
        N = mu3.*(An + H)./(mu4+mu4*2*mu3);
    end
    
    %k+1步的v
    v = max(abs(u - f -bita5./mu5) - lamda./mu5,0).* (u - f -bita5./mu5) ./ abs(u - f -bita5./mu5 + eps) ;
    
    %k+1步的bita
    gra_u = cat(3,ux,uy);
    bita1 = bita1 + mu1*(w - gra_u);
    bita2 = bita2 + mu2*(sqrt(w(:,:,1).^2 + w(:,:,2).^2) - N(:,:,1).*w(:,:,1)-N(:,:,2).*w(:,:,2));
    bita3 = bita3 + mu3*(Q - N(:,:,1) - N(:,:,2));
    bita4 = bita4 + mu4*(M - N);
    bita5 = bita5 + mu5*(v - u + f);
    
%   判断能量函数
    Ax = ux ./ sqrt(ux.^2 + uy.^2);
    Ay = uy ./ sqrt(ux.^2 + uy.^2);
    E(step) = sum(sum(lamda.*abs(divergence(Ax,Ay)) + 0.5*(u-f)));
    if step >2
        result(step) = abs((E(step) - E(step-1) )/ E(step-1));
        if abs((E(step) - E(step-1) )/ E(step)) < 0.0001
            break;
        end
    end
end
subplot(2,2,3);imshow(uint8(u));title('result');
subplot(2,2,4);plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
PSNR = psnr(uint8(u),uint8(image))
PSNR2 = psnr(uint8(f),uint8(image))











