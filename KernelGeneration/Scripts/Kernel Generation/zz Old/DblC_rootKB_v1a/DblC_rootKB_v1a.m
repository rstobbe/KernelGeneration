%==================================================
% (v1a) 
%   
%==================================================

function [SCRPTipt,SCRPTGBL,err] = DblC_rootKB_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

KB_W = str2double(SCRPTipt(strcmp('KB_Width',{SCRPTipt.labelstr})).entrystr);
KB_beta = str2double(SCRPTipt(strcmp('KB_Beta',{SCRPTipt.labelstr})).entrystr);
KB_res = str2double(SCRPTipt(strcmp('KB_Res',{SCRPTipt.labelstr})).entrystr); 

Kern_zWadd = str2double(SCRPTipt(strcmp('Kern_zWadd',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('TestForSS',{SCRPTipt.labelstr})).entrystr); 
Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

Kern_W = KB_W;
Kern_res = KB_res;
Kern_zW = KB_W + Kern_zWadd; 
KRNprms.prec = 'singles';

if rem((KB_W/2),KB_res)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(1,KB_res)
    err.flag = 1;
    err.msg = '1/KernRes is not an integer';
    return
end
if rem(round(1e9*(1/(KB_res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end
zWtest = 2*(ceil((KB_W*SS-2)/2)+1)/SS;
if zWtest > Kern_zW
    err.flag = 1;
    err.msg = 'Kern Zero-Fill Too Small';
    return
end

%--------------------------------------------
% Generate KB Kernel Shape
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);
res0 = 0.001;
u0 = (-KB_W/2:res0:KB_W/2-res0);  
M = KB_beta * sqrt(1 - (2*u0/KB_W).^2);
KB_1D = besseli(0,M);
KB_1D = KB_1D/max(KB_1D);
figure(100); hold on; plot(u0,KB_1D,'r');

%--------------------------------------------
% Plot Full Kernel
%--------------------------------------------
FTKB_1D = fft(fftshift(KB_1D));

figure(101); hold on; 
f = (-1/(2*res0):1/KB_W:1/(2*res0)-1/KB_W);
plot(f,real(fftshift(FTKB_1D)/max(abs(FTKB_1D(:)))),'r');
plot(f,imag(fftshift(FTKB_1D)/max(abs(FTKB_1D(:)))),'r:');
plot([-0.5 -0.5],[0 1],'k:'); plot([0.5 0.5],[0 1],'k:'); xlim([-4 4]);
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

%--------------------------------------------
% Create Root Kernel
%--------------------------------------------
FTKB_1D = real(FTKB_1D);
FTKB_1D(FTKB_1D<0) = 0;
rootFTKB_1D = FTKB_1D.^(0.5);
rootFTKB_1D(FTKB_1D<=5e-4) = 0;

plot(f,real(fftshift(rootFTKB_1D)/max(abs(rootFTKB_1D(:)))),'g');
plot(f,imag(fftshift(rootFTKB_1D)/max(abs(rootFTKB_1D(:)))),'g:');

rootKB_1D = ifftshift(ifft(rootFTKB_1D));
rootKB_1D = rootKB_1D/max(rootKB_1D);

figure(100); hold on; 
plot(u0,rootKB_1D,'g');

rootFTKB_1D = fft(fftshift(rootKB_1D));
figure(101); hold on; 
plot(f,real(fftshift(rootFTKB_1D)/max(abs(rootFTKB_1D(:)))),'g');
plot(f,imag(fftshift(rootFTKB_1D)/max(abs(rootFTKB_1D(:)))),'g:');

%--------------------------------------------
% Generate Full Kaiser Kernel
%--------------------------------------------
Status2('busy','Generate Full Kaiser Kernel',1);
u = (0:Kern_res:(Kern_W/2));                                                 
KB_1D = interp1(u0,KB_1D,u,'linear','extrap');
figure(100); plot(u,KB_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    KBKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    KBKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'double');
end
N = length(KB_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            KBKern(a,b,c) = KB_1D(a)*KB_1D(b)*KB_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
KBKern = KBKern / max(max(max(KBKern)));

%--------------------------------------------
% Generate Root Kaiser Kernel
%--------------------------------------------
Status2('busy','Generate Root Kaiser Kernel',1);
u = (0:Kern_res:(Kern_W/2));                                                 
rootKB_1D = interp1(u0,rootKB_1D,u,'linear','extrap');
figure(100); plot(u,rootKB_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    rootKBKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    rootKBKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'double');
end
N = length(rootKB_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rootKBKern(a,b,c) = rootKB_1D(a)*rootKB_1D(b)*rootKB_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
rootKBKern = rootKBKern / max(max(max(rootKBKern)));

%--------------------------------------------
% Save
%--------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'DblC_rootKB_v1';
    KRNprms.DblKern.convscaleval = 1;
    KRNprms.DblKern.OKforSS = SS;
    KRNprms.DblKern.W = Kern_W;
    KRNprms.DblKern.iKern = 1/Kern_res;
    KRNprms.DblKern.res = Kern_res;
    KRNprms.DblKern.zW = Kern_zW;
    KRNprms.DblKern.Kern = KBKern;
    
    KRNprms.FwdKern.W = KB_W;
    KRNprms.FwdKern.beta = KB_beta;
    KRNprms.FwdKern.iKern = 1/KB_res;
    KRNprms.FwdKern.res = KB_res;
    KRNprms.FwdKern.zW = Kern_zW;
    KRNprms.FwdKern.Kern = rootKBKern;
    
    KRNprms.RvsKern.W = KB_W;
    KRNprms.RvsKern.beta = KB_beta;
    KRNprms.RvsKern.iKern = 1/KB_res;
    KRNprms.RvsKern.res = KB_res;
    KRNprms.RvsKern.zW = Kern_zW;
    KRNprms.RvsKern.Kern = rootKBKern;    
  
    [file,path] = uiputfile('*.mat','Save Convolution Kernel','D:\1 Scripts\zs Shared\zy Convolution Kernels\');  
    KRNprms.name = file;

    save([path,file],'KRNprms');
    [Out,err] = ExternalSave(SCRPTGBL.scrptnum,SCRPTGBL.scrpt,path);
    SCRPTipt = Out.saveSCRPTipt;
    saveSCRPTipt = Out.saveSCRPTipt;
    saveSCRPTGBL = Out.saveSCRPTGBL;
    saveSCRPTIPTGBL = Out.saveSCRPTIPTGBL;
    saveSCRPTPATHS = Out.saveSCRPTPATHS;
    saveScript = Out.saveScript;
    save([path,'\Scrpt_',KRNprms.name],'saveSCRPTipt','saveSCRPTGBL','saveSCRPTIPTGBL','saveSCRPTPATHS','saveScript');    
end