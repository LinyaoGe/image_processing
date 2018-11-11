image = imread('lena.jpg');
image = imnoise(image,'gauss');
subplot(1,2,1);
imshow(image);
image = double(image);

%参数设置
lamda = 25;
[M,N] = size(image);
bita1X = zeros(M,N);
bita1Y = zeros(M,N);
bita2 = zeros(M,N);
Wx = zeros(M,N);
Wy = zeros(M,N);
V = zeros(M,N);
% fai = v.^2;
% div_fai = 2.*v;
mu = 0.5;
u = image;


for step = 1:100
    %处理u的
    W = divergence(Wx,Wy);
    div_bita1 = divergence(bita1X,bita1Y);
    [Ux,Uy]=gradient(u);
    
    %更新bita
    bita1X = bita1X + mu*(W - Ux);
    bita1Y = bita1Y + mu*(W - Uy);
    bita2 = bita2 + mu*(V - W);
    
    F1 = image - mu*W - div_bita1;
    u = center_diff(u);
    u = (F1 + mu*u)./(1+4*mu);
    
    %对新的u求梯度
    [Ux,Uy] = gradient(u);
    [Vx,Vy] = gradient(V);
    [B2x,B2y] = gradient(bita2);
    F2x = -mu*Ux + bita1X + mu*Vx + B2x;
    F2y = -mu*Uy + bita1Y + mu*Vy + B2y;
    
    %对W
    AW1 = center_diff(Wx);
    AW2 = center_diff(Wy);
    Wx = (mu*AW1 - F2x) ./ (5*mu);
    Wy = (mu*AW2 - F2y) ./ (5*mu);
    
    %更新V
    W = divergence(Wx,Wy);
    V = max(abs(W - bita2/mu) - lamda/mu,0) .* (W - bita2*mu) ./ abs(W - bita2*mu);
%     V =  max(abs(W - bita2/mu) - lamda.*div_fai/mu,0) .* (W - bita2*mu) ./ abs(W - bita2*mu);
%     [f1 f2 Af2]=grad(fai);
%     abs_fai=sqrt(f1.^2+f2.^2);
%     E(step)=sum(sum(lamda.*abs_fai+0.5.*(u-image).^2));    
    u2 = divergence(Ux,Uy);
%     E(step)=0.5*sum(sum((u-image).^2))+lamda*sum(sum(Ux.^2+Uy.^2));
    E(step) = sum(lamda*sum(abs(u2)) + 0.5*sum(u-image).^2);
    if step > 2
        abs(E(step) - E(step-1) / E(step)) < 0.001;
        break;
    end
end
subplot(1,2,2);
imshow(uint8(u));











