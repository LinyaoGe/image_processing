clc;
clear;
f = imread('lena.jpg');
subplot(2,2,1);imshow(f);
f = imnoise(f, 'gaussian',0,0.01);
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

%mu的值,lamda的值
mu1 = 0.02;
mu2 = 0.3;
mu3 = 0.002;
mu4 = 2;
lamda = 25;

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
    bita3 = bita3 + mu3*(q - m);
    bita4x = bita4x + mu4*(Mx - Nx);
    bita4y = bita4y + mu4*(My - Ny);
    
    %k+1步的u
    div_w = divergence(Wx,Wy);
    div_bita1 = divergence(bita1x,bita1y);
    Au = center_diff(u);
    u = (f - mu1*div_w + mu1*Au - div_bita1) ./ (1+4*mu1);
    
    %k+1步的w
    [ux,uy] = gradient(u);
    A = ux - ((mu2 + bita2).* Mx - bita1x)./mu1 ;
    B = uy - ((mu2 + bita2).* My - bita1y)./mu1 ;
    Wx = max(abs(A) - (mu2 + bita2)./mu1,0) .* A ./ abs(A + eps);
    Wy = max(abs(B) - (mu2 + bita2)./mu1,0) .* B ./ abs(B + eps);
    
    %k+1步的q
    q = max(abs(divergence(Nx,Ny) - bita3./mu3) - lamda./mu3,0) .* (divergence(Nx,Ny) - bita3./mu3) ./ abs(divergence(Nx,Ny) - bita3./mu3+ eps) ;
    
    %k+1步的m
    Mx = Nx + ((mu2 + bita2).*Wx + bita4x)./mu4;
    My = Ny + ((mu2 + bita2).*Wy + bita4y)./mu4;
    
    Mx = Mx./max(abs(Mx),1);
    My = My./max(abs(My),1);
    
    %k+1步的n
    [qx,qy] = gradient(q);
    [bita3x,bita3y] = gradient(bita3);
    h1 = -qx - bita3x./mu3 + mu4./mu3.*Mx - bita4x./mu3;
    h2 = -qy - bita3y./mu3 + mu4./mu3.*My - bita4x./mu3;
    [nxx,nxy] = gradient(Nx);
    [nyx,nyy] = gradient(Ny);
    nxx = padarray(nxx,[1,1],'replicate');
    nyy = padarray(nyy,[1,1],'replicate');
    nyx = padarray(nyx,[1,1],'replicate');
    nxy = padarray(nxy,[1,1],'replicate');
    for i = 2:M+1
        for j = 2:N+1
            nxx(i,j) = nxx(i+1,j) + nxx(i-1,j) - 2*nxx(i,j);
            nyy(i,j) = nyy(i+1,j) + nyy(i-1,j) - 2*nyy(i,j);
            nyx(i,j) = nyx(i,j+1) - nyx(i,j) - nyx(i-1,j+1) + nyx(i-1,j);
            nxy(i,j) = nxy(i+1,j) - nxy(i,j) - nxy(i+1,j-1) + nxy(i,j-1);
        end
    end
    Nx = nxx(2:M+1,2:N+1) + nyx(2:M+1,2:N+1);
    Ny = nyy(2:M+1,2:N+1) + nxy(2:M+1,2:N+1);
    
%   判断能量函数
    A = ux ./ sqrt(ux.^2 + uy.^2);
    B = uy ./ sqrt(ux.^2 + uy.^2);
    E(step) = sum(lamda.*sum(abs(divergence(A,B))) + 0.5.*sum(u-f));
    if step >2
        if abs((E(step) - E(step-1) )/ E(step-1)) < 0.01
            break;
        end
    end
end
subplot(2,2,3);imshow(uint8(u));title('result');
subplot(2,2,4);plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
PSNR = psnr(u,f)














