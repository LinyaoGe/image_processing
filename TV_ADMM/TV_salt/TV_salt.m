%����ͼ��
image = imread('lena.jpg');
subplot(2,2,1);imshow(image);title('ԭͼ');

%��ͼ�������
f = imnoise(image,'salt & pepper');subplot(2,2,2);imshow(f);title('����ͼ��');

%��ʼ������
f = double(f);
u = f;
[m,n] = size(f);
w = zeros(m,n,2);
v = zeros(m,n);
bita1 = zeros(m,n,2);
bita2 = zeros(m,n);
mu1 = 0.05;
mu2 = 1;
lamda = 0.1;
iteration = 10000;

for step = 1:iteration
    %k+1����u
    div_w = divergence(w(:,:,1),w(:,:,2));
    div_bita1 = divergence(bita1(:,:,1),bita1(:,:,2));
    u = (f + v + (bita2 + div_bita1)./mu2 - mu1/mu2.*(div_w - center_diff(u))) ./ (mu2 + 4*mu1*mu2);
    
    %k+1����w
    [ux,uy] = gradient(u);
    grad_u = cat(3,ux,uy);
    q = grad_u - bita1./mu1;
    abs_q = sqrt(q(:,:,1).^2 + q(:,:,2).^2);
    for i = 1:2
        w(:,:,1) = max(abs_q - lamda./mu1,0).*q(:,:,1)./(abs_q + eps);
        w(:,:,2) = max(abs_q - lamda./mu1,0).*q(:,:,2)./(abs_q + eps);
    end
    %k+1����v
    l = u -f - bita2./mu2;
    v = max(abs(l)-0.5*mu2,0).*l./(abs(l)+eps);
    
    %k+1����bita
    bita1 = bita1 + mu1*(w - grad_u);
    bita2 = bita2 + mu2*(v - u + f);
    
    %��������
    abs_u = sqrt(grad_u(:,:,1).^2 + grad_u(:,:,2).^2);
    E(step) = sum(sum(lamda*abs_u + 0.5*(u-f)));
    
    %stop iteration
    if step > 1
        result(step) = abs((E(step)-E(step-1))/E(step));
        if abs((E(step)-E(step-1))/E(step)) < 0.00001
            break;
        end
    end
end

subplot(2,2,3);imshow(uint8(u));title('�����');
subplot(2,2,4);plot(E);
PSNR = psnr(uint8(u),image)







