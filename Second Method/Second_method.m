image = imread('lena.jpg');
subplot(2,2,1);imshow(image);
image = imnoise(image,'gaussian',0,0.01);
subplot(2,2,2);imshow(image);

%初始化参数
image = double(image);
[M,N] = size(image);
u = image;
Wx = zeros(M,N);
Wy = zeros(M,N);
V = zeros(M,N);
bita1X = zeros(M,N);
bita1Y = zeros(M,N);
bita2 = zeros(M,N);
mu = 0.5;lamda = 25;
total = 100;

%迭代
for step = 1:total
    [ux,uy] = gradient(u);
    Au = center_diff(u);
    div_w  = divergence(Wx,Wy);
    div_bita1 = divergence(bita1X,bita1Y);
    %更新bita
    bita1X = bita1X + mu*(Wx-ux);
    bita1Y = bita1Y + mu*(Wy-uy);
    bita2 = bita2 + mu*(V - div_w);
    u = (image + mu*Au - mu*div_w - div_bita1) ./ (1 + 4*mu);
    [ux,uy] = gradient(u);
    [b2x,b2y] = gradient(bita2);
    [Vx,Vy] = gradient(V);
    %迭代Wx,Wy
    AW1 = center_diff(Wx);
    AW2 = center_diff(Wy);
    Wx = (mu*ux - bita1X - b2x - mu*Vx + mu*AW1)./(5*mu);
    Wy = (mu*uy - bita1Y - b2y - mu*Vy + mu*AW2)./(5*mu);
    div_w = divergence(Wx,Wy);
    V = max(div_w - bita2/mu,0) .* (div_w - bita2/mu) ./ abs(div_w - bita2/mu + eps);
    E(step) = sum(lamda*sum(abs(divergence(ux,uy))) + 0.5*sum((u-image).^2));
    if (step>2)
        if abs((E(step)-E(step-1))/E(step))<0.001
            break;
        end
    end
end


subplot(2,2,3);imshow(uint8(u));title('result');
subplot(2,2,4);plot(E);xlabel('Iterations');ylabel('Energy');legend('Energy/Iterations');
PSNR=psnr(u,image)









