%=========================================================
% Generate Inverse Filter for Kaiser_Bessel Regridding
%=========================================================

function Generate_IF_Odd

%==================
W = 4;
res = W/16;
beta = 6.5;
order = 0;
ZF = 129;
%==================
name = strcat('W',num2str(W,2),'_B',num2str(beta,2),'_ZF',num2str(ZF));
name = regexprep(name,'\.', 'p');

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
KB = KaiBesImg_Full(W,beta,res,order);
np = W/res + 1;

%-------------------------------------------------
% Shape of Filter in Image-Space (full width)
%-------------------------------------------------  
zf1 = 257;
KB = ifftshift(KB);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(KB,zf1);
V = fftshift(fftn(V));                              
figure(2);
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = [-FOV+(1/res)/(2*zf1):(1/res)/zf1:FOV-(1/res)/(2*zf1)];
plot(t,real(V(:,(zf1+1)/2,(zf1+1)/2)));                                   % imaginary component should be ~ zero
%plot(t,log(abs(V(:,zf1/2+1,zf1/2+1))));
hold on
title('Convolution Roll-off');
drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = floor(zf1/2*(1-res)+1);
c = floor(zf1/2*(1+res)+1);
V2 = V(b:c,b:c,b:c);
np2 = c-b+1;
figure(3);
t = (1:np2);
plot(t,real(V2(:,(np2+1)/2,(np2+1)/2)),'c');
%plot(t,log(abs(V2(:,np2/2+1,np2/2+1))));
hold on
title('Convolution Roll-off Within Image Space');
drawnow;

if rem(np2+1,2)
    return
end
%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
V2 = ifftshift(V2);
V2 = fftn(V2);

if length(V2(:,1,1)) > ZF
    error
end
V = Kzerofill_isotropic_Kodd(V2,ZF);
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
%plot(t,10*log(V(:,(ZF+1)/2,(ZF+1)/2)));
%hold on
%title('Log Convolution Roll-off Across Full Image Space');
%drawnow;

figure(5);
t = [1:ZF];
plot(t,V(:,(ZF+1)/2,(ZF+1)/2),'c');
set(gca,'ylim',[0 1]);
hold on
title('Convolution Roll-off Across Full Image Space');
drawnow;

button = questdlg('Save?');
if strcmp(button,'Yes')
    file = strcat('D:\8 Programs\2 NL-PRODS\Inversion Filters\',name);
    [file,path] = uiputfile(file);
    save(strcat(path,file),'V','ZF');
end
