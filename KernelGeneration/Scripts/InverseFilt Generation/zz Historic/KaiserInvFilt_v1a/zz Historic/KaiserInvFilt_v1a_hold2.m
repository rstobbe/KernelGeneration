%==================================================
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = KaiserInvFilt_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

W = str2double(SCRPTipt(strcmp('Width',{SCRPTipt.labelstr})).entrystr);
beta = str2double(SCRPTipt(strcmp('Beta',{SCRPTipt.labelstr})).entrystr); 
ZF = str2double(SCRPTipt(strcmp('ZF',{SCRPTipt.labelstr})).entrystr); 
res = str2double(SCRPTipt(strcmp('Res',{SCRPTipt.labelstr})).entrystr); 
zf1 = str2double(SCRPTipt(strcmp('ZFinter',{SCRPTipt.labelstr})).entrystr); 

Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

if rem(ZF,2)
    err.flag = 1;
    err.msg = 'ZF must be even';
    return
end
if rem((W/2),res)
    err.flag = 1;
    err.msg = 'Width/2 not a multiple of Res';
    return
end
%if rem((zf1*res),2)
%    err.flag = 1;
%    err.msg = 'ZFinter*Res must be odd';
%    return
%end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',1);
u = abs(-(W/2):res:(W/2)-res);                                              % defined over |u| = W/2;     
M = beta * sqrt(1 - (2*u/W).^2);
KB_1D = besseli(0,M);

N = length(KB_1D);
KB = zeros(N,N,N);
for a = 1:N
    for b = 1:N
        for c = 1:N
            %KB(a,b,c) = KB_1D(a)*KB_1D(b)*KB_1D(c);
            KB(a,b,c) = 1;
        end
    end
    Status2('busy',num2str(N-a),2); 
end
Status2('done','',2);
KB = KB / max(max(max(KB)));

figure(1);
cen = length(KB)/2+1;
plot(squeeze(KB(cen,cen,:)));                                   
drawnow;

%-------------------------------------------------
% Shift and Zero-Fill
%-------------------------------------------------  
KB = ifftshift(KB); 
V = zerofill_isotropic_even_doubles(KB,zf1);

figure(2);
plot(squeeze(V(1,1,:)));                                   
drawnow;

%-------------------------------------------------
% FT
%-------------------------------------------------  
V = fftshift(ifftn(V));                              

figure(3); hold on
FOV = 1/(2*res);                                                        % Filter V2 = 1/res times larger than image
t = (-FOV:(1/res)/zf1:FOV-(1/res)/zf1);
plot(t,imag(V(:,zf1/2+1,zf1/2+1)),'b');                                 % imaginary component should be ~ zero
plot(t,real(V(:,zf1/2+1,zf1/2+1)),'r');                                 % real component - 'minimal fold-over at edges'
plot(t,squeeze(real(V(zf1/2+1,:,zf1/2+1))),'r'); 
plot(t,squeeze(real(V(zf1/2+1,zf1/2+1,:))),'k-*'); 
drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = (zf1/2)*(1-res)+1.5;
c = (zf1/2)*(1+res)+0.5;
V2 = V(b:c,b:c,b:c);
V2 = V2/max(max(max(V2)));
np2 = length(V2);

figure(4); hold on
t2 = t(b:c);
%plot(t2,real(V2(:,np2/2+1,np2/2+1)),'r-*');
plot(t2,real(V2(:,(np2+1)/2,(np2+1)/2)),'r-*');
hold on
title('Convolution Roll-off Within Image Space');
drawnow;

%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
V2 = fftshift(fftn(ifftshift(V2)));
np3 = length(V2);

figure(5);
%plot(1:np3,real(V2(:,np3/2+1,np3/2+1)),'r-');
plot(1:np3,real(V2(:,(np3+1)/2,(np3+1)/2)),'r-');
hold on
title('Kernel inter');
drawnow;

V = zerofill_isotropic_odd_doubles(ifftshift(V2),ZF);

figure(6); hold on
plot(1:ZF,real(V(:,1,1)),'r-*');
hold on
title('Zerofill Kernel');
drawnow;

V = fftshift(ifftn(V));
V = V/max(max(max(V)));

figure(4);
t = (-0.5:1/ZF:0.5-1/ZF);
%t = 1:ZF
plot(t,V(:,ZF/2+1,ZF/2+1),'c');

button = questdlg('Save?');
if strcmp(button,'Yes')
    file = strcat('D:\0 Programs\A0 Inverse Filters\Current Files\',name);
    save(file,'V','ZF');
end



%==================
W = 6;                  
beta = 10.5;              
ZF = 168;
%==================
if beta > 10
    name = strcat('W',num2str(W,2),'_B',num2str(beta,3),'_ZF',num2str(ZF));
else
    name = strcat('W',num2str(W,2),'_B',num2str(beta,2),'_ZF',num2str(ZF));
end
name = regexprep(name,'\.', 'p');

