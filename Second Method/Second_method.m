image = imread('lena.jpg');
image = imnoise(image,'gauss');
subplot(1,2,1);
imshow(image);
image = double(image);

%参数设置
lamda = 30;
[M,N] = size(image);
bita1X = zeros(M,N);
bita1Y = zeros(M,N);
bita2 = zeros(M,N);
Wx = zeros(M,N);
Wy = zeros(M,N);
V = zeros(M,N);
s = 256*256;arf = 0.000003;arf2 = 0.01;
mu1 = s*arf;
mu2 = s*arf2;
u = image;


for step = 1:100
    %处理u的
    W = divergence(Wx,Wy);
    thita1 = divergence(bita1X,bita1Y);
    F1 = image - mu1*W - thita1;
    u = center_diff(u);
    u = (F1 + mu1*u)./(1+4*mu1);
    
    %对新的u求梯度
    [Ux,Uy] = gradient(u);
    [Vx,Vy] = gradient(V);
    [B2x,B2y] = gradient(bita2);
    F2x = -mu1*Ux + bita1X + mu2*Vx + B2x;
    F2y = -mu2*Uy + bita1Y + mu2*Vy + B2y;
    
    Wx = padarray(Wx,[1,1],'replicate');
    Wy = padarray(Wy,[1,1],'replicate');
    %更新W
    for i = 2:257
        for j = 2:257
            Wx(i,j) = mu2 * (Wx(i+1,j) - Wx(i,j) - Wx(i+1,j+1) - Wx(i,j-1) + Wy(i,j+1) + Wy(i,j-1) - 2*Wy(i,j));
            Wy(i,j) = mu2 * (Wy(i,j+1) - Wy(i,j) - Wy(i-1,j+1) - Wy(i-1,j) + Wx(i+1,j) + Wx(i-1,j) - 2*Wx(i,j));
        end
    end
    Wx = Wx(2:257,2:257);
    Wy = Wy(2:257,2:257);
    
    %更新V
    W = divergence(Wx,Wy);
    V = max(abs(W - bita2/mu2) - ((-bita2)/mu2 - (V - W)*abs(V)/V),0)*(W - bita2*mu2) / abs(W - bita2*mu2);
    
    %更新bita
    bita1X = bita1X + mu1*(W - Ux);
    bita1Y = bita1Y + mu1*(W - Uy);
    bita2 = bita2 + mu2*(V - W);
    u2 = divergence(Ux,Uy);
    E(step) = sum(lamda*sum(abs(u2)) + 0.5*sum(u-image).^2);
    if step > 2
        abs(E(step) - E(step-1) / E(step-1)) < 0.0000000001;
        break;
    end
end
subplot(1,2,2);
imshow(uint8(u));











