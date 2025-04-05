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
% FF(aa<0)=0;
tic
aa=1-(lambda*FX).^2-(lambda*FY).^2;
H=exp(1i*k*d*sqrt(aa)); 
H1=H;
H1(aa<0)=0;
I1=abs(fftshift(ifft2(fftshift(H1).*U1))).^2;
%get_img(I1,"Icx");
H2=H;
H2(aa>0)=0;
I2=abs(fftshift(ifft2(fftshift(H2).*U1))).^2;
%get_img(I2,"Isx");
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
I3=abs(fftshift(ifft2(fftshift(H3).*U1))).^2;
%get_img(I3,"Icz");
H4=H;
H4(aa>0)=0;
I4=abs(fftshift(ifft2(fftshift(H4).*U1))).^2;
toc
% get_img(I4,"Isz");
Isum=I1+I2+I3+I4;
all=sum(sum(I2));
get_img(Isum,"Isum");%end

%statistics
cc=[sum(sum(I1)) sum(sum(I2)) sum(sum(I3)) sum(sum(I4))]./sum(sum(Isum));
dd=[get_ctsum(I1) get_ctsum(I2) get_ctsum(I3) get_ctsum(I4)];
ee=dd./sum(dd);
ksum=get_ctsum(Isum)/sum(sum(Isum));
function y=get_ctsum(I1)
all=sum(sum(I1));
center=sum(sum(I1(513:1024,513:1024)));
y=all-center;
end
