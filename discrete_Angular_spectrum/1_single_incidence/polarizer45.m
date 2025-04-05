%获得透射矩阵
I2=zeros(2172,2172);
I_check=zeros(4879,4879);%3450sqrt(2)
x1=-140*17+70;
x2=140*17-70;
x=x1:140:x2;
for i=1:4879
    for j=1:4879
        for num=1:34
            if abs(i-x(num)-2439.5)<=35
                I_check(i,j)=1;
            end
        end
    end
end
for i=1:4879
    for j=1:4879
        if -i+j-2128.5>0
            I_check(i,j)=0;
        end
        if -i+j+2128.5<0
            I_check(i,j)=0;
        end
        if i+j-7007.5>0%4879-311+4879/2
            I_check(i,j)=0;
        end
        if i+j-2750.5<0%311+4879/2
            I_check(i,j)=0;
        end
    end
end
%imwrite(I_check,"mn.png");
K=4879/724;
for i=1:724
    for j=1:724
        ii=round(K*(i-1)+1);
        jj=round(K*(j-1)+1);
        if I_check(ii,jj)==1
            I2(i+724,j+724)=1;
        end
    end
end
I_temp2=imrotate(I2,45,'bilinear','crop');
imwrite(I_temp2,"mn_scaled.png");

U1=fft2(fftshift(I2));
SL=3.45*3*sqrt(2.0);% Side length um
lambda=0.555;% wavelength
k=2*pi/lambda;% wavevector
N=2172;% samples for side length
d=0.7;% The distance between the input and output planes
dx=SL/N;
fx=-1/(2*dx):1/SL:1/(2*dx)-1/SL;  %freq coords
[FX,FY]=meshgrid(fx,fx);
% FF(aa<0)=0;
aa=1-(lambda*FX).^2-(lambda*FY).^2;
H=exp(1i*k*d*sqrt(aa)); 
H1=H;
H1(aa<0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H1).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I1=I_temp(319:1854,319:1854);
get_img(I1,"Icx");
H2=H;
H2(aa>0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H2).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I2=I_temp(319:1854,319:1854);
get_img(I2,"Isx");
% FF=ones(1536,1536);%观察频域
% hsv=zeros(1536,1536,3);
% FF(aa<0)=0;
% hsv(:,:,3)=FF;
% rgb=hsv2rgb(hsv);
% imwrite(rgb,'FF.png');%观察频域
bb=sqrt(aa);
H=exp(1i*k*d*bb)*lambda.*FY./bb; 
H3=H;
H3(aa<0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H3).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I3=I_temp(319:1854,319:1854);
get_img(I3,"Icz");
H4=H;
H4(aa>0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H4).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I4=I_temp(319:1854,319:1854);
get_img(I4,"Isz");
Isum=I1+I2+I3+I4;
get_img(Isum,"Isum");
cc=[sum(sum(I1)) sum(sum(I2)) sum(sum(I3)) sum(sum(I4))]./sum(sum(Isum));
dd=[get_ctsum(I1) get_ctsum(I2) get_ctsum(I3) get_ctsum(I4)];
ee=dd./sum(dd);
k1=sum(sum(Isum(1:512,1:512)))/sum(sum(Isum));
k2=sum(sum(Isum(1:512,513:1024)))/sum(sum(Isum));
k3=sum(sum(Isum(1:512,1025:1536)))/sum(sum(Isum));
ksum=get_ctsum(Isum)/sum(sum(Isum));
function y=get_ctsum(I1)
all=sum(sum(I1));
center=sum(sum(I1(513:1024,513:1024)));
y=all-center;
end
