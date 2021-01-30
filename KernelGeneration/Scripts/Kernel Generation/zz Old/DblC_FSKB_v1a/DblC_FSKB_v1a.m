%==================================================
% 
%==================================================

function [SCRPTipt,SCRPTGBL,err] = DblC_FSKB_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

FS_W = str2double(SCRPTipt(strcmp('FS_Width',{SCRPTipt.labelstr})).entrystr);
FS_BW = str2double(SCRPTipt(strcmp('FS_BW',{SCRPTipt.labelstr})).entrystr); 
FS_beta = str2double(SCRPTipt(strcmp('FS_Beta',{SCRPTipt.labelstr})).entrystr);

KB_W = str2double(SCRPTipt(strcmp('KB_Width',{SCRPTipt.labelstr})).entrystr);
KB_beta = str2double(SCRPTipt(strcmp('KB_Beta',{SCRPTipt.labelstr})).entrystr);

Kern_res = str2double(SCRPTipt(strcmp('Kern_Res',{SCRPTipt.labelstr})).entrystr); 
Kern_zWadd = str2double(SCRPTipt(strcmp('Kern_zWadd',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('TestForSS',{SCRPTipt.labelstr})).entrystr); 
Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

if FS_W >= KB_W
    Kern_W = FS_W;
    Kern_zW = FS_W + Kern_zWadd;
else
    err.flag = 1;
    err.msg = 'KB width should not be larger than FS width';
    return
end    
KRNprms.prec = 'singles';

if rem((FS_W/2),Kern_res) || rem((KB_W/2),Kern_res)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(1,Kern_res)
    err.flag = 1;
    err.msg = '1/KernRes is not an integer';
    return
end
if rem(round(1e9*(1/(Kern_res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end

%--------------------------------------------
% Generate FS and KB Kernel Shapes
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);
res0 = 0.001;
u0 = (-FS_W/2:res0:FS_W/2-res0);   
Sinc_1D = sinc(u0*FS_BW);
M = FS_beta * sqrt(1 - (2*u0/FS_W).^2);
FSKB_1D = besseli(0,M);
FS_1D = Sinc_1D.*FSKB_1D;
FS_1D = FS_1D/max(FS_1D);
figure(100); hold on; plot(u0,FS_1D,'b');

KBu0 = (-KB_W/2:res0:KB_W/2-res0);  
M = KB_beta * sqrt(1 - (2*KBu0/KB_W).^2);
KB_1D = besseli(0,M);
KB_1D = KB_1D/max(KB_1D);

%T_W = 1;
%KB_1D = [(0:res0/(T_W/2):1) (1-res0/(T_W/2):-res0/(T_W/2):res0/(T_W/2))];

lendif = length(FS_1D) - length(KB_1D);
ZFKB_1D = [zeros(1,lendif/2) KB_1D zeros(1,lendif/2)];
figure(100); hold on; plot(u0,ZFKB_1D,'r');

%--------------------------------------------
% Plot Kernel Frequency Profiles
%--------------------------------------------
FTFS_1D = fft(fftshift(FS_1D));
FTZFKB_1D = fft(fftshift(ZFKB_1D));

figure(101); hold on; 
f = (-1/(2*res0):1/FS_W:1/(2*res0)-1/FS_W);
plot(f,real(fftshift(FTFS_1D)/max(abs(FTFS_1D(:)))),'b');
plot(f,imag(fftshift(FTFS_1D)/max(abs(FTFS_1D(:)))),'b');
plot(f,real(fftshift(FTZFKB_1D)/max(abs(FTZFKB_1D(:)))),'r');
plot(f,imag(fftshift(FTZFKB_1D)/max(abs(FTZFKB_1D(:)))),'r');
plot([-0.5 -0.5],[0 1],'k:'); plot([0.5 0.5],[0 1],'k:'); xlim([-4 4]);
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

%--------------------------------------------
% Create Combined Kernel
%--------------------------------------------
FSDBL_1D = real(FTFS_1D .* FTZFKB_1D);
DBL_1D = ifftshift(ifft(FSDBL_1D));
DBL_1D = DBL_1D/max(DBL_1D);

figure(100); hold on; plot(u0,DBL_1D,'g');
FTDBL_1D = fft(fftshift(DBL_1D));
figure(101); hold on; 
plot(f,real(fftshift(FTDBL_1D)/max(abs(FTDBL_1D(:)))),'g');
plot(f,imag(fftshift(FTDBL_1D)/max(abs(FTDBL_1D(:)))),'g');

%--------------------------------------------
% Generate Combined Kernel
%--------------------------------------------
Status2('busy','Generate Combined Convolution Kernel',1);
u = (0:Kern_res:(Kern_W/2));                                                 
DBL_1D = interp1(u0,DBL_1D,u,'linear','extrap');
figure(100); plot(u,DBL_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    DblKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    DblKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'double');
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
% Generate FS Kernel
%--------------------------------------------
Status2('busy','Generate Filtered-Sinc Kernel',1);
u = (0:Kern_res:(Kern_W/2));                                                 
FS_1D = interp1(u0,FS_1D,u,'linear','extrap');
figure(100); plot(u,FS_1D,'k:');
if strcmp(KRNprms.prec,'singles')
    FSKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    FSKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'double');
end
N = length(FS_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            FSKern(a,b,c) = FS_1D(a)*FS_1D(b)*FS_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
FSKern = FSKern / max(max(max(FSKern)));

%--------------------------------------------
% Generate Kaiser Kernel
%--------------------------------------------
Status2('busy','Generate Kaiser Kernel',1);
u = (0:Kern_res:(KB_W/2));                                                 
KB_1D = interp1(u0,ZFKB_1D,u,'linear','extrap');
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

button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'DblC_FSKB_v1';
    KRNprms.DblKern.convscaleval = 1;
    KRNprms.DblKern.OKforSS = SS;
    KRNprms.DblKern.W = Kern_W;
    KRNprms.DblKern.iKern = 1/Kern_res;
    KRNprms.DblKern.res = Kern_res;
    KRNprms.DblKern.zW = Kern_zW;
    KRNprms.DblKern.Kern = DblKern;
    
    KRNprms.FwdKern.W = Kern_W;
    KRNprms.FwdKern.BW = FS_BW;
    KRNprms.FwdKern.beta = FS_beta;
    KRNprms.FwdKern.iKern = 1/Kern_res;
    KRNprms.FwdKern.res = Kern_res;
    KRNprms.FwdKern.zW = Kern_zW;
    KRNprms.FwdKern.Kern = FSKern;
    
    KRNprms.RvsKern.W = Kern_W;
    KRNprms.RvsKern.beta = KB_beta;
    KRNprms.RvsKern.iKern = 1/Kern_res;
    KRNprms.RvsKern.res = Kern_res;
    KRNprms.RvsKern.zW = Kern_zW;
    KRNprms.RvsKern.Kern = KBKern;    
  
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