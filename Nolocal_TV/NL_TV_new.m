clc;clear;close all;
%read the image
image = imread('lena.jpg');

%add the noise
% image_noise = imnoise(image,'gaussian',0.01);
image_noise = imnoise(image,'poisson');
figure;
imshow(image_noise);
u = double(image_noise);
f = u;
[m,n] = size(u);

%init output image
output = zeros(m,n);

%init para
h = 12;
gamma = 700;
mu = 0.5;
sigma = 3.5;

%init window&patch size
window_radius = 4;
patch_radius = 3;
patch_size = 2*patch_radius + 1;
total_radius = window_radius + patch_radius;

%init w&lamda
w = cell(m,n,2);
lamda = cell(m,n,2);
for i = 1:m
    for j = 1:n
        w{i,j,1} = zeros(2*window_radius+1,2*window_radius+1);
        w{i,j,2} = zeros(2*window_radius+1,2*window_radius+1);
        lamda{i,j,1} = zeros(2*window_radius+1,2*window_radius+1);
        lamda{i,j,2} = zeros(2*window_radius+1,2*window_radius+1);
    end
end

%set Gaussian convolution
kernal = fspecial('gaussian',patch_size,sigma);

%Nolocal method
for step = 1:2
    u2 = padarray(u,[total_radius,total_radius],'symmetric');
    %calculate k
    for i =1: m
        for j =1: n
            i1 = i + total_radius;
            j1 = j + total_radius;
            %init sum
            sum_k = 0;sum_y_k = 0;
            sum_w_k = 0;sum_lamda_k = 0;

            %chose w&lamda
            w11 = w{i,j,1};w12 = w{i,j,2};
            lamda11 = lamda{i,j,1};lamda12 = lamda{i,j,2};
            %x patch
            X = u2(i1-patch_radius:i1+patch_radius,j1-patch_radius:j1+patch_radius);
            p=0;
            for r = i1-window_radius:i1+window_radius
                p = p+1;
                q = 0;
                for s = j1-window_radius:j1+window_radius
                    q = q+1;
                    %y patch
                    Y = u2(r-patch_radius:r+patch_radius,s-patch_radius:s+patch_radius);
                    d = sum(sum(kernal.*(X - Y).^2));
                    k = exp(-d./(h.^2));
                    sum_k = sum_k + k;
                    sum_y_k = sum_y_k + u2(r,s)*k;
                    sum_w_k = sum_w_k + (w11(p,q) - w12(p,q)).*sqrt(k);
                    sum_lamda_k = sum_lamda_k + (lamda11(p,q) - lamda12(p,q)).*sqrt(k);
                end
            end
            output(i,j) = (f(i,j) + 2*mu*sum_y_k - mu*sum_w_k - sum_lamda_k)/(1+2*mu*sum_k);
            
            %update w&lamda
            div_u = zeros(2*window_radius+1,2*window_radius+1,2);
            p=0;
            for r = i1-window_radius:i1+window_radius
                p = p+1;
                q = 0;
                for s = j1+window_radius:j1+window_radius
                    q = q+1;
                    Y = u2(r-patch_radius:r+patch_radius,s-patch_radius:s+patch_radius);
                    d = sum(sum(kernal.*(X - Y).^2));
                    d2 = sum(sum(kernal.*(Y - X).^2));
                    k = exp(-d./(h.^2));
                    k2 = exp(-d2/(h.^2));
                    y = u2(r,s);
                    div_u(p,q,1) = (y - output(i,j)).*sqrt(k);
                    div_u(p,q,2) = (output(i,j) - y).*sqrt(k2);
                end
            end
            abs = sqrt((div_u(:,:,1)-lamda{i,j,1}).^2 + (div_u(:,:,2)-lamda{i,j,2}).^2);
            w{i,j,1} = max(abs - gamma/mu,0).*(div_u(:,:,1)-lamda11)./(abs+eps);
            w{i,j,2} = max(abs - gamma/mu,0).*(div_u(:,:,2)-lamda12)./(abs+eps);
            
            lamda{i,j,1} = lamda{i,j,1} + mu*(w{i,j,1} - div_u(:,:,1));
            lamda{i,j,2} = lamda{i,j,2} + mu*(w{i,j,1} - div_u(:,:,2));
        end
    end
    u = output;
end
figure;imshow(uint8(u));
PSNR1 = psnr(uint8(u),image)
PSNR2 = psnr(uint8(image_noise),image)


















