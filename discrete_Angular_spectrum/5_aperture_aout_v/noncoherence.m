Ixy=zeros(1536,1536);
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
            Ixy(i+512,j+512)=1;
        end
    end
end

fx=1:1536;  %freq coords
fx2=0.001*K*fx;
[FX,FY]=meshgrid(fx2,fx2);
SL=3.45*3;% Side length um
lambda=0.555;% wavelength
k=2*pi/lambda;% wavevector
N=1536;% samples for side length
distance=0.7;% The distance between the input and output planes
dx=SL/N;
fx=-1/(2*dx):1/SL:1/(2*dx)-1/SL;  %freq coords
[FX2,FY2]=meshgrid(fx,fx);
aa=1-(lambda*FX2).^2-(lambda*FY2).^2;
H12=exp(1i*k*distance*sqrt(aa));

H1=H12;
H1(aa<0)=0;
H2=H12;
H2(aa>0)=0;
bb=sqrt(aa);
H34=exp(1i*k*distance*bb)*lambda.*FY2./bb;
H3=H34;
H3(aa<0)=0;
H4=H34;
H4(aa>0)=0;
%---------------------------------------------光圈计算

f=35;
r_array=[1.1 1.6 2.2 3.1 4.4 6.25 8.75 10.95];
for v=35:35
    Isum=zeros(1536,1536);
    for number=1:8
        r=r_array(number);
        x=-r:0.05:r;
        [X,Y]=meshgrid(x,x);
        d=X.^2+Y.^2;
        d(d>r^2)=0;
        if number>1
            r0=r_array(number-1);
            d(d<=r0^2)=0;
        end
        numZeros =length(x)^2-numel(find(d == 0));

        %---------------------------------------------光圈计算
        t=0;
        for t1=1:length(x)
            for t2=1:length(x)
                if d(t1,t2)~=0
                    tanb1=x(t1)/v;
                    tanb2=x(t2)/v;
                    gama=sqrt(tanb1^2+tanb2^2+1);
                    cosxi1=tanb1/gama;
                    cosxi2=tanb2/gama;
                    I_xy2=Ixy.*exp(1i*2*pi*(cosxi1*FX+cosxi2*FY)/0.555);
                    U1=fft2(fftshift(I_xy2));
                    %----------------------------------------------------------------


                    I1=abs(fftshift(ifft2(fftshift(H1).*U1))).^2;
                    I2=abs(fftshift(ifft2(fftshift(H2).*U1))).^2;
                    I3=abs(fftshift(ifft2(fftshift(H3).*U1))).^2;
                    I4=abs(fftshift(ifft2(fftshift(H4).*U1))).^2;

                    Isum=Isum+I1+I2+I3+I4;
                    fprintf('%d/%d----%d\n',t,numZeros,number);
                    t=t+1;
                end
            end
        end
        ename=strcat('Isum',num2str(v),"_",num2str(number),'.mat');
        save(ename,'Isum');
    end
end