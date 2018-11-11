%¶ÁÈ¡Í¼Ïñ
Image = imread('lena.jpg');

%¸øÍ¼Ïñ¼ÓÔëÉù
f = imnoise(Image,'gauss');
subplot(1,2,1)
imshow(f);
hold on;
f = double(f);

%ÉèÖÃ²ÎÊı
u=f;
h = 0.5;
lamda = 0.8;

%±éÀú
for step = 1:200
    u = center_diff(u);
    u = (f + lamda/h^2*u) / (1+4*lamda/h^2);    
    [Fx,Fy] = gradient(u);
    E(step) = sum(0.5 * sum((u-f).^2) + 0.5*lamda*sum(Fx.^2 + Fy.^2));
    if step > 2
        if abs((E(step)-E(step-1))/E(step)) < 0.01
            break;
        end
    end
end
subplot(1,2,2)
imshow(uint8(u));







