
%-----------------------------------------
% Zero - Fill K-Space
%-----------------------------------------

function KZ = Kzerofill_isotropic(K,ZF)

[x,y,z] = size(K);          % these values should always be odd (center of k-space @ (x+1)/2)
KZ = complex(zeros(ZF,ZF,ZF,'single'),zeros(ZF,ZF,ZF,'single'));

%KZ(1:x,1:y,1:z) = K;       % get same 'abs' image just doing this, but 'phase' image will be toast (i.e. ringing at pixel frequency) 

KZ(1:(x+1)/2,1:(y+1)/2,1:(z+1)/2) = K(1:(x+1)/2,1:(y+1)/2,1:(z+1)/2);
KZ(ZF-(x+1)/2+2:ZF,ZF-(y+1)/2+2:ZF,ZF-(z+1)/2+2:ZF) = K((x+1)/2+1:x,(y+1)/2+1:y,(z+1)/2+1:z);
KZ(1:(x+1)/2,ZF-(y+1)/2+2:ZF,ZF-(z+1)/2+2:ZF) = K(1:(x+1)/2,(y+1)/2+1:y,(z+1)/2+1:z);
KZ(ZF-(x+1)/2+2:ZF,1:(y+1)/2,ZF-(z+1)/2+2:ZF) = K((x+1)/2+1:x,1:(y+1)/2,(z+1)/2+1:z);
KZ(ZF-(x+1)/2+2:ZF,ZF-(y+1)/2+2:ZF,1:(z+1)/2) = K((x+1)/2+1:x,(y+1)/2+1:y,1:(z+1)/2);
KZ(1:(x+1)/2,1:(y+1)/2,ZF-(z+1)/2+2:ZF) = K(1:(x+1)/2,1:(y+1)/2,(z+1)/2+1:z);
KZ(ZF-(x+1)/2+2:ZF,1:(y+1)/2,1:(z+1)/2) = K((x+1)/2+1:x,1:(y+1)/2,1:(z+1)/2);
KZ(1:(x+1)/2,ZF-(y+1)/2+2:ZF,1:(z+1)/2) = K(1:(x+1)/2,(y+1)/2+1:y,1:(z+1)/2);

%figure(2);
%t = [1:ZF];
%plot(t,KZ(:,1,1));
%test = 0;

