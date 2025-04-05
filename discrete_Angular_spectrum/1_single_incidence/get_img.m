function y=get_img(I,name) 
    sizeI=size(I,1);
    I_max=max(max(I));
    I_min=min(min(I));
    fprintf('%f %f\n',I_min,I_max);
    e=(I-I_min)/(0.05*I_max-I_min);
    e(e>1)=1;
    e(e<0)=0;
    h=abs(e);
    hsv=ones(sizeI,sizeI,3);
    hsv(:,:,1)=0.6-0.6*h;
    hsv(:,:,3)=0.6*h+0.4;
    onethird=512;
    twothird=1024;
    hsv(onethird:onethird+1,:,3)=1;
    hsv(twothird:twothird+1,:,3)=1;
    hsv(:,onethird:onethird+1,3)=1;
    hsv(:,twothird:twothird+1,3)=1;
    hsv(onethird:onethird+1,:,2)=0;
    hsv(twothird:twothird+1,:,2)=0;
    hsv(:,onethird:onethird+1,2)=0;
    hsv(:,twothird:twothird+1,2)=0;
    rgb=hsv2rgb(hsv);
    imwrite(rgb,name+".png");
end