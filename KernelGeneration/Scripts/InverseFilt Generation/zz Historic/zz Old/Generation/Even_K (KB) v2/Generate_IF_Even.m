%=========================================================
% Generate Inverse Filter for Kaiser_Bessel Regridding
%=========================================================

function Generate_IF_Even

%==================
W = 6;                  
beta = 10.5;              
order = 0;
ZF = 168;
%==================
if beta > 10
    name = strcat('W',num2str(W,2),'_B',num2str(beta,3),'_ZF',num2str(ZF));
else
    name = strcat('W',num2str(W,2),'_B',num2str(beta,2),'_ZF',num2str(ZF));
end
name = regexprep(name,'\.', 'p');

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
%res = 0.25;
res = 0.125;
%res = 0.0625;
KB = KaiBesImg_Full(W,beta,res,order);

%-------------------------------------------------
% Shape of Filter in Image-Space (full width)
%-------------------------------------------------  
zf1 = 256;
KB = ifftshift(KB);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(KB,zf1);
V = fftshift(fftn(V));                              

%figure(2);
%FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
%t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1);
%plot(t,imag(V(:,zf1/2+1,zf1/2+1)));                                   % imaginary component should be ~ zero
%drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = floor(zf1/2*(1-res)+1);
c = floor(zf1/2*(1+res));
np2 = zf1*res;
V2 = V(b:c,b:c,b:c);

figure(3);
t = [1:np2];
plot(t,real(V2(:,np2/2+1,np2/2+1)));
%plot(t,log(abs(V2(:,np2/2+1,np2/2+1))));
hold on
title('Convolution Roll-off Within Image Space');
drawnow;

if rem(np2,2)
    error
    return
end
%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
V2 = ifftshift(V2);
V2 = fftn(V2);

V = complex(zeros(ZF,ZF,ZF,'single'),zeros(ZF,ZF,ZF,'single'));
V(1:np2/2,1:np2/2,1:np2/2) = V2(1:np2/2,1:np2/2,1:np2/2);
V(ZF-np2/2+1:ZF,ZF-np2/2+1:ZF,ZF-np2/2+1:ZF) = V2(np2/2+1:np2,np2/2+1:np2,np2/2+1:np2);
V(1:np2/2,ZF-np2/2+1:ZF,ZF-np2/2+1:ZF) = V2(1:np2/2,np2/2+1:np2,np2/2+1:np2);
V(ZF-np2/2+1:ZF,1:np2/2,ZF-np2/2+1:ZF) = V2(np2/2+1:np2,1:np2/2,np2/2+1:np2);
V(ZF-np2/2+1:ZF,ZF-np2/2+1:ZF,1:np2/2) = V2(np2/2+1:np2,np2/2+1:np2,1:np2/2);
V(1:np2/2,1:np2/2,ZF-np2/2+1:ZF) = V2(1:np2/2,1:np2/2,np2/2+1:np2);
V(ZF-np2/2+1:ZF,1:np2/2,1:np2/2) = V2(np2/2+1:np2,1:np2/2,1:np2/2);
V(1:np2/2,ZF-np2/2+1:ZF,1:np2/2) = V2(1:np2/2,np2/2+1:np2,1:np2/2);
clear V2;

V = ifftn(V);
V = fftshift(V);
V = V/max(max(max(V)));
V = abs(V);                                                         % imaginary component ~ zero

%-----------------------------------------------
% filter manipulation
%V = 1./V;
%V2 = V(ZF/4+1:3*ZF/4,ZF/4+1:3*ZF/4,ZF/4+1:3*ZF/4); 
%vmax = max(max(max(V2)));
%V = ones(ZF,ZF,ZF,'single')*vmax;
%V(ZF/4+1:3*ZF/4,ZF/4+1:3*ZF/4,ZF/4+1:3*ZF/4) = V2;

%figure(4);
%t = [1:ZF];
%plot(t,10*log(V(:,ZF/2+1,ZF/2+1)));
%hold on
%title('Log Convolution Roll-off Across Full Image Space');
%drawnow;

figure(5);
t = [1:ZF];
plot(t,V(:,ZF/2+1,ZF/2+1),'c');
set(gca,'ylim',[0 1]);
hold on
title('Convolution Roll-off Across Full Image Space');
drawnow;

button = questdlg('Save?');
if strcmp(button,'Yes')
    file = strcat('D:\0 Programs\A0 Inverse Filters\Current Files\',name);
    save(file,'V','ZF');
end
