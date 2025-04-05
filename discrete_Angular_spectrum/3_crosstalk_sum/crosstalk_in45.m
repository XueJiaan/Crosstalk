%0像素，只考虑邻边像素，包含2个45，2个135
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
%imwrite(I2,"mn_scaled.png");

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
H2=H;
H2(aa>0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H2).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I2=I_temp(319:1854,319:1854);
bb=sqrt(aa);
H=exp(1i*k*d*bb)*lambda.*FY./bb; 
H3=H;
H3(aa<0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H3).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I3=I_temp(319:1854,319:1854);
H4=H;
H4(aa>0)=0;
I_rotated=abs(fftshift(ifft2(fftshift(H4).*U1))).^2;
I_temp=imrotate(I_rotated,45,'bilinear','crop');
I4=I_temp(319:1854,319:1854);
Isum45=I1+I2+I3+I4;

center=zeros(512,512);
center=center+0.707*Isum0(513:1024,1:512);
center=center+0.707*Isum0(513:1024,1025:1536);
center=center+0.707*flip(Isum0(1:512,513:1024),2);
center=center+0.707*flip(Isum0(1025:1536,513:1024),2);
;
Imax=Isum45(513:1024,513:1024)+center;
Imin=center+Isum45(1:512,1:512)+Isum45(1:512,1025:1536)+Isum45(1025:1536,1:512)+Isum45(1025:1536,1025:1536);

I_max=max(max(Imin));
I_min=min(min(Imin));
fprintf('%f %f\n',I_min,I_max);
e=(center-I_min)/(0.5*I_max-I_min);
e(e>1)=1;
e(e<0)=0;
h=abs(e);
hsv=ones(512,512,3);
hsv(:,:,1)=0.6-0.6*h;
hsv(:,:,3)=0.6*h+0.4;
rgb=hsv2rgb(hsv);
imwrite(rgb,"crosstalk_in0.png");

e=sum(sum(Imax))/sum(sum(Imin));
fprintf('e= %f\n',e);