%=========================================================
% Generate Inverse Filter for Kaiser_Bessel Regridding
%=========================================================

function Generate_IF

%==================
W = 4;
beta = 6.5;
order = 0;
ZF = 312;
%==================
name = strcat('W',num2str(W,2),'_B',num2str(beta,2),'_ZF',num2str(ZF));
name = regexprep(name,'\.', 'p');

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
res = 0.25;
KB = KaiBesImg_Full(W,beta,res,order);
np = W/res + 1;

%-------------------------------------------------
% Shape of Filter in Image-Space (full width)
%-------------------------------------------------  
KB = ifftshift(KB);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(KB,ZF);
V = fftshift(fftn(V));                              
figure(2);
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = [-FOV:(1/res)/ZF:FOV-(1/res)/ZF];
plot(t,real(V(:,ZF/2+1,ZF/2+1)));                                   % imaginary component should be ~ zero
%plot(t,log(abs(V(:,ZF/2+1,ZF/2+1))));
title('Log Convolution Roll-off');
drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = floor(ZF/2*(1-res)+1);
c = floor(ZF/2*(1+res));
np2 = ZF*res;
V2 = V(b:c,b:c,b:c);
%figure(3);
%t = [1:np2];
%V2 = 100*V2/max(max(max(V2)));
%plot(t,real(V2(:,np2/2+1,np2/2+1)));
%plot(t,log(abs(V2(:,np2/2+1,np2/2+1))));
%title('Convolution Roll-off Within Image Space');
%drawnow;

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

figure(4);
t = [1:ZF];
plot(t,10*log(V(:,ZF/2+1,ZF/2+1)));
title('Log Convolution Roll-off Across Full Image Space');
drawnow;

figure(5);
t = [1:ZF];
plot(t,V(:,ZF/2+1,ZF/2+1));
title('Convolution Roll-off Across Full Image Space');
drawnow;

button = questdlg('Save?');
if strcmp(button,'Yes')
    file = strcat('D:\8 Programs\2 NL-PRODS\Inversion Filters\',name);
    save(file,'V','ZF');
end
