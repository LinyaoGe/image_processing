clc;
clear;
image = imread('lena.jpg');
subplot(2,2,1);imshow(image);
f = imnoise(image, 'gaussian',0,0.005);
subplot(2,2,2);imshow(f);

%------------------初始化参数-------------------
f = double(f);
u = f;
[M,N] = size(f);
%bita的值
bita1x = zeros(M,N); 
bita1y = zeros(M,N);
bita2 = zeros(M,N);
bita3 = zeros(M,N);
bita4x = zeros(M,N);
bita4y = zeros(M,N);

% mu的值,lamda的值
mu1 = 0.02;
mu2 = 0.3;
mu3 = 0.0001;
mu4 = 0.0001;
lamda = 1000;
h = 5;
%w的值
Wx = zeros(M,N);
Wy = zeros(M,N); 

%q的值
q = zeros(M,N);

%n 和 m的值
Nx = zeros(M,N); 
Ny = zeros(M,N); 
Mx = zeros(M,N); 
My = zeros(M,N); 

%------------------迭代-------------------
total = 100;
for step = 1:total
    %k+1步的bita
    [ux,uy] = gradient(u);
    m = divergence(Mx,My);
    bita1x = bita1x + mu1*(Wx - ux);
    bita1y = bita1y + mu1*(Wy - uy);
    bita2 = bita2 + mu2*(sqrt(Wx.^2 + Wy.^2) - Nx.*Wx-Ny.*Wy);
    bita3 = bita3 + mu3*(q - Nx - Ny);
    bita4x = bita4x + mu4*(Mx - Nx);
    bita4y = bita4y + mu4*(My - Ny);
    
    %k+1步的u
    div_w = divergence(Wx,Wy);
    div_bita1 = divergence(bita1x,bita1y);
    for i = 1:4
        Au = center_diff(u);
        u = (f - mu1*div_w + mu1./h^2*Au - div_bita1) ./ (1+4*mu1/h^2);
    end
    
    %k+1步的w
    [ux,uy] = gradient(u);
    Ax = ux + ((mu2 + bita2).* Mx - bita1x)./mu1 ;
    Ay = uy + ((mu2 + bita2).* My - bita1y)./mu1 ;
    abs_A = sqrt(Ax.^2 + Ay.^2);
    Wx = max(abs_A - (mu2 + bita2)./mu1,0) .* abs_A ./ abs(abs_A + eps);
    Wy = max(abs_A - (mu2 + bita2)./mu1,0) .* abs_A ./ abs(abs_A + eps);
    
    %k+1步的q
    q = max(abs(divergence(Nx,Ny) - bita3./mu3) - lamda./mu3,0) .* (divergence(Nx,Ny) - bita3./mu3) ./ abs(divergence(Nx,Ny) - bita3./mu3+ eps) ;
    
%     k+1步的m
    Mx = Nx + ((mu2 + bita2).*Wx + bita4x)./mu4;
    My = Ny + ((mu2 + bita2).*Wy + bita4y)./mu4;
    abs_m = sqrt(Mx.^2 + My.^2);
    Mx = Mx ./ max(abs_m,1);
    My = My ./ max(abs_m,1);
%     Mx = Mx./max(abs(Mx),1);
%     My = My./max(abs(My),1);
    
    %k+1步的n
    [qx,qy] = gradient(q);
    [bita3x,bita3y] = gradient(bita3);
    h1 = -qx - bita3x./mu3 + mu4./mu3.*Mx - bita4x./mu3;
    h2 = -qy - bita3y./mu3 + mu4./mu3.*My - bita4x./mu3;
    for i = 1:4
        [nx,ny] = n_diff(Nx,Ny);
        Nx = mu3./h.^2 .* (nx + h1)./(mu4+2*mu3./h^2);
        Ny = mu3./h.^2 .* (ny + h2)./(mu4+2*mu3./h^2);
%         Nx = mu3./h.^2 .* (center_diff(Nx) + h1)./(mu4+2*mu3./h^2);
%         Ny = mu3./h.^2 .* (center_diff(Ny) + h2)./(mu4+2*mu3./h^2);
    end
    
    
%   判断能量函数
    Ax = ux ./ sqrt(ux.^2 + uy.^2);
    Ay = uy ./ sqrt(ux.^2 + uy.^2);
%     abs_A = sqrt(Ax.^2 + Ay.^2);
%     E(step) = sum(lamda.*sum(abs_A) + 0.5*sum((u-f).^2));
    E(step) = sum(lamda.*sum(abs(divergence(Ax,Ay))) + 0.5.*sum((u-f).^2));
    if step >2
        result(step) = abs((E(step) - E(step-1) )/ E(step-1));
        if abs((E(step) - E(step-1) )/ E(step)) < 0.001
            break;
        end
    end
end
subplot(2,2,3);imshow(uint8(u));title('result');
subplot(2,2,4);plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
PSNR = psnr(uint8(u),uint8(image))
PSNR2 = psnr(uint8(f),uint8(image))














