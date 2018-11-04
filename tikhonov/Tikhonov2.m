clc
p1=imread('lena.jpg');
imshow(p1);
p=imnoise(p1,'gauss');
figure;
imshow(p);
[m,n]=size(p);
p=double(p);
%u=zeros(m,n);
u=p;
lamda=0.8;
h=0.5;
for k=1:8
    for i=2:255;
        for j=2:255;
            u(i,j)=(p(i,j)+lamda/h^2*(u(i-1,j)+u(i,j-1)+u(i,j+1)+u(i+1,j)))/(1+4*lamda/h^2);
        end
    end
end
figure;
imshow(uint8(u));
            