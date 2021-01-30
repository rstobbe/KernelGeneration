%==================================================
% (v1a) 
%   
%==================================================

function [SCRPTipt,SCRPTGBL,err] = DblC_FSFS_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

FS1_W = str2double(SCRPTipt(strcmp('FS1_Width',{SCRPTipt.labelstr})).entrystr);
FS1_BW = str2double(SCRPTipt(strcmp('FS1_BW',{SCRPTipt.labelstr})).entrystr); 
FS1_beta = str2double(SCRPTipt(strcmp('FS1_Beta',{SCRPTipt.labelstr})).entrystr);
FS1_res = str2double(SCRPTipt(strcmp('FS1_Res',{SCRPTipt.labelstr})).entrystr); 

FS2_W = str2double(SCRPTipt(strcmp('FS2_Width',{SCRPTipt.labelstr})).entrystr);
FS2_BW = str2double(SCRPTipt(strcmp('FS2_BW',{SCRPTipt.labelstr})).entrystr); 
FS2_beta = str2double(SCRPTipt(strcmp('FS2_Beta',{SCRPTipt.labelstr})).entrystr);
FS2_res = str2double(SCRPTipt(strcmp('FS2_Res',{SCRPTipt.labelstr})).entrystr); 

Kern_zWadd = str2double(SCRPTipt(strcmp('Kern_zWadd',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('TestForSS',{SCRPTipt.labelstr})).entrystr); 
Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

if FS1_W >= FS2_W
    Kern_W = FS1_W;
    Kern_zW = FS1_W + Kern_zWadd;
    FS1_zW = FS1_W + Kern_zWadd;
    FS2_zW = FS2_W + Kern_zWadd;
else
    err.flag = 1;
    err.msg = 'FS2 width should not be larger than FS1 width';
    return
end    
KRNprms.prec = 'singles';

if rem((FS1_W/2),FS1_res) || rem((FS2_W/2),FS2_res)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(1,FS1_res) || rem(1,FS2_res)
    err.flag = 1;
    err.msg = '1/KernRes is not an integer';
    return
end
if rem(round(1e9*(1/(FS1_res*SS)))/1e9,1) || rem(round(1e9*(1/(FS2_res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end
zWtest = 2*(ceil((FS2_W*SS-2)/2)+1)/SS;
if zWtest > FS2_zW
    err.flag = 1;
    err.msg = 'FS2 Zero-Fill Too Small';
    return
end
zWtest = 2*(ceil((FS1_W*SS-2)/2)+1)/SS;
if zWtest > FS1_zW
    err.flag = 1;
    err.msg = 'FS1 Zero-Fill Too Small';
    return
end
if FS2_zW/FS2_res > 1100
    err.flag = 1;
    err.msg = ['FS2 kernel size: ',num2str(FS2_zW/FS1_res)];
    return
end
if FS1_zW/FS1_res > 1100
    err.flag = 1;
    err.msg = ['FS1 kernel size: ',num2str(FS1_zW/FS1_res)];
    return
end

%--------------------------------------------
% Generate FS1 and FS2 Kernel Shapes
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);
res0 = 0.001;
u1 = (-FS1_W/2:res0:FS1_W/2-res0);   
Sinc_1D = sinc(u1*FS1_BW);
M = FS1_beta * sqrt(1 - (2*u1/FS1_W).^2);
FS1KB_1D = besseli(0,M);
FS1_1D = Sinc_1D.*FS1KB_1D;
FS1_1D = FS1_1D/max(FS1_1D);
figure(100); hold on; plot(u1,FS1_1D,'b');

u2 = (-FS2_W/2:res0:FS2_W/2-res0);   
Sinc_1D = sinc(u2*FS2_BW);
M = FS2_beta * sqrt(1 - (2*u2/FS2_W).^2);
FS2KB_1D = besseli(0,M);
FS2_1D = Sinc_1D.*FS2KB_1D;
FS2_1D = FS2_1D/max(FS2_1D);

lendif = length(FS1_1D) - length(FS2_1D);
ZFFS2_1D = [zeros(1,lendif/2) FS2_1D zeros(1,lendif/2)];
figure(100); hold on; plot(u1,ZFFS2_1D,'r');

%--------------------------------------------
% Plot Kernel Frequency Profiles
%--------------------------------------------
FTFS1_1D = fft(fftshift(FS1_1D));
FTZFFS2_1D = fft(fftshift(ZFFS2_1D));

figure(101); hold on; 
f = (-1/(2*res0):1/FS1_W:1/(2*res0)-1/FS1_W);
plot(f,real(fftshift(FTFS1_1D)/max(abs(FTFS1_1D(:)))),'b');
plot(f,imag(fftshift(FTFS1_1D)/max(abs(FTFS1_1D(:)))),'b');
plot(f,real(fftshift(FTZFFS2_1D)/max(abs(FTZFFS2_1D(:)))),'r');
plot(f,imag(fftshift(FTZFFS2_1D)/max(abs(FTZFFS2_1D(:)))),'r');
plot([-0.5 -0.5],[0 1],'k:'); plot([0.5 0.5],[0 1],'k:'); xlim([-4 4]);
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

%--------------------------------------------
% Create Combined Kernel
%--------------------------------------------
FTDBL_1D = real(FTFS1_1D .* FTZFFS2_1D);
DBL_1D = ifftshift(ifft(FTDBL_1D));
DBL_1D = DBL_1D/max(DBL_1D);

figure(100); hold on; plot(u1,DBL_1D,'g');
FTDBL_1D = fft(fftshift(DBL_1D));
figure(101); hold on; 
plot(f,real(fftshift(FTDBL_1D)/max(abs(FTDBL_1D(:)))),'g');
plot(f,imag(fftshift(FTDBL_1D)/max(abs(FTDBL_1D(:)))),'g');

%--------------------------------------------
% Generate Combined Kernel
%--------------------------------------------
Status2('busy','Generate Combined Convolution Kernel',1);
u = (0:FS1_res:(Kern_W/2));                                                 
DBL_1D = interp1(u1,DBL_1D,u,'linear','extrap');
figure(100); plot(u,DBL_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    DblKern = zeros((Kern_zW/2)*round(1/FS1_res)+1,(Kern_zW/2)*round(1/FS1_res)+1,(Kern_zW/2)*round(1/FS1_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    DblKern = zeros((Kern_zW/2)*round(1/FS1_res)+1,(Kern_zW/2)*round(1/FS1_res)+1,(Kern_zW/2)*round(1/FS1_res)+1,'double');
end
N = length(DBL_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            DblKern(a,b,c) = DBL_1D(a)*DBL_1D(b)*DBL_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
DblKern = DblKern / max(max(max(DblKern)));

%--------------------------------------------
% Generate FS1 Kernel
%--------------------------------------------
Status2('busy','Generate Filtered-Sinc1 Kernel',1);
u = (0:FS1_res:(Kern_W/2));                                                 
FS1_1D = interp1(u1,FS1_1D,u,'linear','extrap');
figure(100); plot(u,FS1_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    FS1Kern = zeros((FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    FS1Kern = zeros((FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,'double');
end
N = length(FS1_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            FS1Kern(a,b,c) = FS1_1D(a)*FS1_1D(b)*FS1_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
FS1Kern = FS1Kern / max(max(max(FS1Kern)));

%--------------------------------------------
% Generate FS2 Kernel
%--------------------------------------------
Status2('busy','Generate Filtered-Sinc2 Kernel',1);
u = (0:FS2_res:(Kern_W/2));                                                 
FS2_1D = interp1(u2,FS2_1D,u,'linear','extrap');
figure(100); plot(u,FS2_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    FS2Kern = zeros((FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    FS2Kern = zeros((FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,'double');
end
N = length(FS2_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            FS2Kern(a,b,c) = FS2_1D(a)*FS2_1D(b)*FS2_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
FS2Kern = FS2Kern / max(max(max(FS2Kern)));

button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'DblC_FSFS_v1';
    KRNprms.DblKern.convscaleval = 1;
    KRNprms.DblKern.OKforSS = SS;
    KRNprms.DblKern.W = Kern_W;
    KRNprms.DblKern.iKern = 1/FS1_res;
    KRNprms.DblKern.res = FS1_res;
    KRNprms.DblKern.zW = Kern_zW;
    KRNprms.DblKern.Kern = DblKern;
    
    KRNprms.FwdKern.W = FS1_W;
    KRNprms.FwdKern.BW = FS1_BW;
    KRNprms.FwdKern.beta = FS1_beta;
    KRNprms.FwdKern.iKern = 1/FS1_res;
    KRNprms.FwdKern.res = FS1_res;
    KRNprms.FwdKern.zW = FS1_zW;
    KRNprms.FwdKern.Kern = FS1Kern;
    
    KRNprms.RvsKern.W = FS2_W;
    KRNprms.RvsKern.BW = FS2_BW;
    KRNprms.RvsKern.beta = FS2_beta;
    KRNprms.RvsKern.iKern = 1/FS2_res;
    KRNprms.RvsKern.res = FS2_res;
    KRNprms.RvsKern.zW = FS2_zW;
    KRNprms.RvsKern.Kern = FS2Kern;   
  
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