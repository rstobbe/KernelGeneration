%==================================================
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = KaiserInvFilt_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

KernWidth = str2double(SCRPTipt(strcmp('KernWidth',{SCRPTipt.labelstr})).entrystr);
KernBeta = str2double(SCRPTipt(strcmp('KernBeta',{SCRPTipt.labelstr})).entrystr); 
KernRes = str2double(SCRPTipt(strcmp('KernRes',{SCRPTipt.labelstr})).entrystr); 
SS = str2double(SCRPTipt(strcmp('SubSamp',{SCRPTipt.labelstr})).entrystr); 
ZF = str2double(SCRPTipt(strcmp('ZF',{SCRPTipt.labelstr})).entrystr); 
res1 = str2double(SCRPTipt(strcmp('Res_dev',{SCRPTipt.labelstr})).entrystr); 
zf1 = str2double(SCRPTipt(strcmp('ZF_dev',{SCRPTipt.labelstr})).entrystr); 

Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

if rem(ZF,2)
    err.flag = 1;
    err.msg = 'ZF must be even';
    return
end
if rem(1,KernRes)
    err.flag = 1;
    err.msg = '1/KernRes is not an integer';
    return
end
if rem((KernWidth/2),KernRes)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(round(1e9*(1/(KernRes*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end
if rem(((KernWidth*SS)/2),res1)
    err.flag = 1;
    err.msg = 'KernWidth*SS/2 not a multiple of Res_dev';
    return
end
if rem((zf1*res1)+1,2)
    err.flag = 1;
    err.msg = 'ZF_dev*Res_dev must be odd';
    return
end

%-------------------------------------------------
% Build Kernel
%-------------------------------------------------  
Status2('busy','Generate Convolution Kernel',1);
W = KernWidth*SS;
u = abs(-(W/2):res1:(W/2)-res1);                                              % defined over |u| = W/2;     
M = KernBeta * sqrt(1 - (2*u/W).^2);
KB_1D = besseli(0,M);

N = length(KB_1D);
KB = zeros(N,N,N);
for a = 1:N
    for b = 1:N
        for c = 1:N
            KB(a,b,c) = KB_1D(a)*KB_1D(b)*KB_1D(c);
            %KB(a,b,c) = 1;
        end
    end
    Status2('busy',num2str(N-a),2); 
end
Status2('done','',2);
KB = KB / max(max(max(KB)));

%figure(1);
%cen = length(KB)/2+1;
%plot(squeeze(KB(cen,cen,:)));                                   
%drawnow;

%-------------------------------------------------
% Shift and Zero-Fill
%-------------------------------------------------  
KB = ifftshift(KB); 
V = zerofill_isotropic_even_doubles(KB,zf1);

%figure(2);
%plot(squeeze(V(1,1,:)));                                   
%drawnow;

%-------------------------------------------------
% FT
%-------------------------------------------------  
Status2('busy','FT',1);
V = fftshift(ifftn(V));                              
FOV = SS/(2*res1);                                                        % Filter V2 = 1/res1 times larger than image
t = (-FOV:(SS/res1)/zf1:FOV-(SS/res1)/zf1);

%figure(3); hold on
%plot(t,imag(V(:,zf1/2+1,zf1/2+1)),'b');                                 % imaginary component should be ~ zero
%plot(t,real(V(:,zf1/2+1,zf1/2+1)),'r');                                 % real component - 'minimal fold-over at edges'
%plot(t,squeeze(real(V(zf1/2+1,:,zf1/2+1))),'r'); 
%plot(t,squeeze(real(V(zf1/2+1,zf1/2+1,:))),'k-*'); 
%drawnow;

%-------------------------------------------------
% Shape of Filter in Image-Space (image width)
%-------------------------------------------------
b = (zf1/2)*(1-res1)+1.5;
c = (zf1/2)*(1+res1)+0.5;
V2 = V(b:c,b:c,b:c);
V2 = V2/max(max(max(V2)));
np2 = length(V2);

figure(4); hold on
t2 = t(b:c);
plot(t2,real(V2(:,(np2+1)/2,(np2+1)/2)),'r-*');
hold on
title('Convolution Roll-off Within Image Space');
drawnow;

%-------------------------------------------------
% Interpolate to Desired Resolution
%-------------------------------------------------
Status2('busy','Interpolate to ZF',1);
V2 = fftshift(fftn(ifftshift(V2)));

%figure(5);
%np3 = length(V2);
%plot(1:np3,real(V2(:,(np3+1)/2,(np3+1)/2)),'r-');
%hold on
%title('Kernel inter');
%drawnow;

V = zerofill_isotropic_odd_doubles(ifftshift(V2),ZF);

%figure(6); hold on
%plot(1:ZF,real(V(:,1,1)),'r-*');
%hold on
%title('Zerofill Kernel');
%drawnow;

V = fftshift(ifftn(V));
V = V/max(max(max(V)));
V = abs(V);

figure(4);
t = (-0.5*SS:SS/ZF:0.5*SS-SS/ZF);
plot(t,V(:,ZF/2+1,ZF/2+1),'c');
plot([-0.5 -0.5],[-0.2 1],'k:'); 
plot([0.5 0.5],[-0.2 1],'k:'); 
plot([-0.5*SS 0.5*SS],[0 0],'k:');
axis([-0.5*SS 0.5*SS -0.1 1]);

%-------------------------------------------------
% Save
%-------------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    IFparams.W = W;
    IFparams.ZF = ZF;
    IFparams.KernWidth = KernWidth;
    IFparams.KernBeta = KernBeta;
    IFparams.KernRes = KernRes;
    IFparams.SubSamp = SS;
    IFparams.Res_dev = res1;
    IFparams.ZF_dev = zf1;

    name = ['IFKBv2_W',num2str(W,3),'_B',num2str(KernBeta,3),'_ZF',num2str(ZF)];
    name = regexprep(name,'\.', 'p');
    
    path = ['D:\1 Scripts\zs Shared\zx Inverse Filters\Kaiser_v2\',name];
    mkdir(path);
    save([path,'\',name],'V','IFparams');
    [Out,err] = ExternalSave(SCRPTGBL.scrptnum,SCRPTGBL.scrpt,path);
    SCRPTipt = Out.saveSCRPTipt;
    saveSCRPTipt = Out.saveSCRPTipt;
    saveSCRPTGBL = Out.saveSCRPTGBL;
    saveSCRPTIPTGBL = Out.saveSCRPTIPTGBL;
    saveSCRPTPATHS = Out.saveSCRPTPATHS;
    saveScript = Out.saveScript;
    save([path,'\Scrpt_',name],'saveSCRPTipt','saveSCRPTGBL','saveSCRPTIPTGBL','saveSCRPTPATHS','saveScript');
end





