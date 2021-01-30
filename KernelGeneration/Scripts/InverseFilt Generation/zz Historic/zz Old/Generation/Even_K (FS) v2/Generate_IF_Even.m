%=========================================================
% Generate Inverse Filter
%=========================================================

function Generate_IF_Even

%==================
W = 12;                  
BW = 0.9;
beta = 6;              
ZF = 96;
subsamp = 1.2;
%==================
if beta > 10
    name = strcat('W',num2str(W,2),'_BW',num2str(BW,3),'_B',num2str(beta,3),'_ZF',num2str(ZF));
else
    name = strcat('W',num2str(W,2),'_BW',num2str(BW,3),'_B',num2str(beta,2),'_ZF',num2str(ZF));
end
name = regexprep(name,'\.', 'p');

%-------------------------------------------------
% Get Filter Shape
%------------------------------------------------- 
%res = 0.25;
res = 0.125;
%res = 0.0625;
[FS] = FilSincImg_Full(W,BW,beta,res);
%FS = ones(size(FS));

%-------------------------------------------------
% Shape of Filter in Image-Space (full width)
%-------------------------------------------------  
zf1 = 384;
FS = ifftshift(FS);                                                 % ifftshift because odd (gives max on left)
V = Kzerofill_isotropic(FS,zf1);
V = fftshift(fftn(V));                              
V = V/max(V(:));

figure(2); hold on;
FOV = 1/(2*res);                                                    % Filter V2 = 1/res times larger than image
t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1);
plot(t,real(V(:,zf1/2+1,zf1/2+1)),'r');                                   
plot(t,imag(V(:,zf1/2+1,zf1/2+1)),'b');                                   % imaginary component should be ~ zero
plot([-0.5 -0.5],[0 1],'k:'); plot([0.5 0.5],[0 1],'k:');
imext = 0.5/subsamp;
plot([-imext -imext],[0 1],'k:'); plot([imext imext],[0 1],'k:');
axis([-1 1 -0.1 1.1]);
title('Convolution Roll-off Across Centre of Image Space');
drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = floor(zf1/2*(1-res)+1);
c = floor(zf1/2*(1+res)+1);
np2 = zf1*res+1;
V2 = V(b:c,b:c,b:c);

figure(3); hold on;
t = (1:np2);
plot(t,real(V2(:,(np2+1)/2,(np2+1)/2)),'b-*');
%plot(t,log(abs(V2(:,np2/2+1,np2/2+1))));
ylim([0 1.1]);
title('Convolution Roll-off Across Centre of Image Space');
drawnow;

if rem(np2+1,2)
    error
    return
end
%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
V2 = fftn(ifftshift(V2));

doshift = 1;
if doshift == 1
    for a = -(np2-1)/2:(np2-1)/2
        for b = -(np2-1)/2:(np2-1)/2
            for c = -(np2-1)/2:(np2-1)/2
                Shift(a+(np2-1)/2+1,b+(np2-1)/2+1,c+(np2-1)/2+1) = exp(1i*pi*0.5*a/np2 + 1i*pi*0.5*b/np2 + 1i*pi*0.5*c/np2);
            end
        end
    end
    Shift = ifftshift(Shift);
    V2 = V2.*Shift;
end

V = Kzerofill_isotropic(V2,ZF);
%figure; hold on; plot(real(squeeze(V(1,1,:))),'r'); plot(imag(squeeze(V(1,1,:))),'b');

V = fftshift(ifftn(V));
V = V/max(max(max(V)));
V = abs(V);                                                         % imaginary component ~ zero

figure(5); hold on;
t = (1:ZF);
plot(t,squeeze(V(:,ZF/2+1,ZF/2+1)),'c'); plot(t,squeeze(V(ZF/2+1,:,ZF/2+1)),'b'); plot(t,squeeze(V(ZF/2+1,ZF/2+1,:)),'r');
ylim([0 1.1]);
title('Interpolated Convolution Roll-off Across Centre of Image Space');
drawnow;

figure(6); hold on;
t = (1:ZF);
edge = round((ZF*(1-1/subsamp))/2+1);
plot(t,squeeze(V(:,edge,edge)),'c'); plot(t,squeeze(V(edge,:,edge)),'b'); plot(t,squeeze(V(edge,edge,:)),'r');
ylim([0 1.1]);
title('Interpolated Convolution Roll-off Across Edge of Image Space');
drawnow;

button = questdlg('Save?');
if strcmp(button,'Yes')
    file = strcat('D:\0 Programs\A0 Inverse Filters\Current Files\',name);
    save(file,'V','ZF');
end
