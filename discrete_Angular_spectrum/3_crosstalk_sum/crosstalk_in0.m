I2=zeros(1536,1536);
I_check=zeros(3450,3450);
x1=-140*11+70;
x2=140*11-70;
x=x1:140:x2;
for i=1:3450
    for j=1:3450
        for num=1:22
            if abs(i-x(num)-1725.5)<=35&&abs(j-1725.5)<=1470
                I_check(i,j)=1;
            end
        end
    end
end
K=3450/512;
for i=1:512
    for j=1:512
        ii=round(K*(i-1)+1);
        jj=round(K*(j-1)+1);
        if I_check(ii,jj)==1
            I2(i+512,j+512)=1;
        end
    end
end
%imwrite(I2,"result1.png");

U1=fft2(fftshift(I2));
SL=3.45*3;% Side length um
lambda=0.555;% wavelength
k=2*pi/lambda;% wavevector
N=1536;% samples for side length
d=0.7;% The distance between the input and output planes
dx=SL/N;
fx=-1/(2*dx):1/SL:1/(2*dx)-1/SL;  %freq coords
[FX,FY]=meshgrid(fx,fx);
aa=1-(lambda*FX).^2-(lambda*FY).^2;
H=exp(1i*k*d*sqrt(aa)); 
H1=H;
H1(aa<0)=0;
I1=abs(fftshift(ifft2(fftshift(H1).*U1))).^2;
H2=H;
H2(aa>0)=0;
I2=abs(fftshift(ifft2(fftshift(H2).*U1))).^2;
bb=sqrt(aa);
H=exp(1i*k*d*bb)*lambda.*FY./bb; 
H3=H;
H3(aa<0)=0;
I3=abs(fftshift(ifft2(fftshift(H3).*U1))).^2;
H4=H;
H4(aa>0)=0;
I4=abs(fftshift(ifft2(fftshift(H4).*U1))).^2;
Isum0=I1+I2+I3+I4;

center=zeros(512,512);
center=center+0.707*Isum45(513:1024,1:512);
center=center+0.707*Isum45(513:1024,1025:1536);
center=center+0.707*flip(Isum45(1:512,513:1024),2);
center=center+0.707*flip(Isum45(1025:1536,513:1024),2);

corner=Isum0(1:512,1:512)+Isum0(1:512,1025:1536)+Isum0(1025:1536,1:512)+Isum0(1025:1536,1025:1536);
Imax=Isum0(513:1024,513:1024)+center;
Imin=center+corner;
e=sum(sum(Imax))/sum(sum(Imin));
fprintf('e= %f\n',e);
fprintf('%f %f\n',sum(sum(corner)),sum(sum(center)));

I_max=max(max(center));
I_min=min(min(center));
fprintf('%f %f\n',I_min,I_max);
e=(center-I_min)/(0.5*I_max-I_min);
e(e>1)=1;
e(e<0)=0;
h=abs(e);
hsv=ones(512,512,3);
hsv(:,:,1)=0.6-0.6*h;
hsv(:,:,3)=0.6*h+0.4;
rgb=hsv2rgb(hsv);
imwrite(rgb,"crosstalk_in45.png");