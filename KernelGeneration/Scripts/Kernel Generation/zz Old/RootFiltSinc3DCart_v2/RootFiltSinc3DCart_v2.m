%==================================================
% (v2) 
%   - Singles/Doubles 
%   - zW (zero-filled kernel width)
%==================================================

function [SCRPTipt,SCRPTGBL,err] = RootFiltSinc3DCart_v2(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

W = str2double(SCRPTipt(strcmp('KernWidth',{SCRPTipt.labelstr})).entrystr);
res = str2double(SCRPTipt(strcmp('KernRes',{SCRPTipt.labelstr})).entrystr); 
BW = str2double(SCRPTipt(strcmp('KernBW',{SCRPTipt.labelstr})).entrystr); 
beta = str2double(SCRPTipt(strcmp('KernBeta',{SCRPTipt.labelstr})).entrystr); 
expenv = str2double(SCRPTipt(strcmp('KernExpEnv',{SCRPTipt.labelstr})).entrystr); 
zW_add = str2double(SCRPTipt(strcmp('zW_add',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('FacilSubSamp',{SCRPTipt.labelstr})).entrystr); 
Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

zW = W + zW_add;
KRNprms.prec = 'singles';

if rem((W/2),res)
    err.flag = 1;
    err.msg = 'KernWidth/2 not a multiple of KernRes';
    return
end
if rem(1,res)
    err.flag = 1;
    err.msg = '1/KernRes is not an integer';
    return
end
if rem(round(1e9*(1/(res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(KernRes*SS) not an integer';
    return
end

%--------------------------------------------
% Generate 'Full' Kernel Shape
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);
res0 = 0.001;
u0 = (-W/2:res0:W/2-res0);   
Sinc_1D = sinc(u0*BW);
M = beta * sqrt(1 - (2*u0/W).^2);
KB_1D = besseli(0,M);
E_1D = exp(-abs(u0)/expenv);
Kern_1D = Sinc_1D.*KB_1D.*E_1D;
Kern_1D = fftshift(Kern_1D);
Kern_1D = Kern_1D/Kern_1D(1);
Kern_1D = ifftshift(Kern_1D);
%figure; plot(Kern_1D);

%--------------------------------------------
% Generate Root Kernel Shape
%--------------------------------------------
FTKern_1D = fft(fftshift(Kern_1D));
%figure; plot(real(fftshift(FTKern_1D)));
%figure; plot(imag(fftshift(FTKern_1D)));
test = min(real(FTKern_1D));
if test < 0
    err.flag = 1;
    err.msg = 'decrease ''expenv''';
    return
end
rootFTKern_1D = (real(FTKern_1D)).^(0.5);
rootKern_1D = ifftshift(ifft(rootFTKern_1D));
%figure; plot(real(rootKern_1D));
%figure; plot(imag(rootKern_1D));
%figure; plot(real(fftshift(fft(fftshift(rootKern_1D)))).^2);
rootKern_1D = real(rootKern_1D)/max(real(rootKern_1D));

%--------------------------------------------
% Generate 'Full' Kernel
%--------------------------------------------
Status2('busy','Generate Non-Root Convolution Kernel',1);
u = (0:res:(W/2));                                                 
Kern_1D = interp1(u0,Kern_1D,u,'linear','extrap');
%figure; plot(u,Kern_1D);

if strcmp(KRNprms.prec,'singles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'double');
end

N = length(Kern_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = Kern_1D(a)*Kern_1D(b)*Kern_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
Kern = Kern / max(max(max(Kern)));

%--------------------------------------------
% Generate Root Kernel
%--------------------------------------------
Status2('busy','Generate Root Convolution Kernel',1);
u = (0:res:(W/2));                                                 
rootKern_1D = interp1(u0,rootKern_1D,u,'linear','extrap');
%figure; plot(u,rootKern_1D);

if strcmp(KRNprms.prec,'singles')
    rootKern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    rootKern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'double');
end

N = length(rootKern_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rootKern(a,b,c) = rootKern_1D(a)*rootKern_1D(b)*rootKern_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
rootKern = rootKern / max(max(max(rootKern)));

%--------------------------------------------
% Visuals
%--------------------------------------------
if strcmp(Vis,'On')
    figure(1100); hold on; 
    plot((1:length(Kern)),squeeze(Kern(1,1,:)),'r'); xlim([0 length(Kern)]); xlabel('kernel size');
    plot((1:length(rootKern)),squeeze(rootKern(1,1,:)),'b');
    figure(1101); hold on; 
    plot((0:length(Kern)-1)*res,squeeze(Kern(1,1,:)),'r'); xlim([0 zW/2]); xlabel('kernel width');
    plot((0:length(rootKern)-1)*res,squeeze(rootKern(1,1,:)),'b'); 
    
    KernProf = [Kern_1D zeros(1,2000) flipdim(Kern_1D(2:length(Kern_1D)),2)];
    rootKernProf = [rootKern_1D zeros(1,2000) flipdim(rootKern_1D(2:length(rootKern_1D)),2)];
    OS = length(KernProf)/(2*length(Kern_1D));
    T1 = ifftshift(ifft(rootKernProf));
    T1 = T1/max(T1);
    sqT1 = T1.^2;
    sqT1 = sqT1/max(sqT1);
    T2 = ifftshift(ifft(KernProf));
    T2 = T2/max(T2);
    
    bot = (length(T1)+1)/2 +1 - (W/2*OS);
    top = (length(T1)+1)/2 +1 + (W/2*OS);
    figure(1102); hold on; plot((1:length(T1)),T1,'r'); plot([bot bot],[0 1],'k:'); plot([top top],[0 1],'k:'); xlim([0.98*bot top/0.98]);
    figure(1103); hold on; 
    plot((1:length(T1)),sqT1,'r'); plot([bot bot],[0 1],'k:'); plot([top top],[0 1],'k:'); xlim([0.98*bot top/0.98]);
    plot((1:length(T1)),T2,'b');
end
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'RootFiltSinc_v2';
    KRNprms.iKern = 1/res;
    KRNprms.W = W;
    KRNprms.res = res;
    KRNprms.zW = zW;
    KRNprms.beta = beta;
    KRNprms.BW = BW;
    KRNprms.ExpEnv = expenv;
    KRNprms.convscaleval = 1;
    KRNprms.SubSamp = SS;
    
    name = ['RFSv2_W',num2str(W,2),'_E',num2str(expenv,3),'_BW',num2str(BW,3),'_B',num2str(beta,3)];
    KRNprms.name = regexprep(name,'\.', 'p');

    path = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\RootFiltSinc_v2\',KRNprms.name];
    mkdir(path);
    save([path,'\',KRNprms.name],'rootKern','Kern','KRNprms');
    [Out,err] = ExternalSave(SCRPTGBL.scrptnum,SCRPTGBL.scrpt,path);
    SCRPTipt = Out.saveSCRPTipt;
    saveSCRPTipt = Out.saveSCRPTipt;
    saveSCRPTGBL = Out.saveSCRPTGBL;
    saveSCRPTIPTGBL = Out.saveSCRPTIPTGBL;
    saveSCRPTPATHS = Out.saveSCRPTPATHS;
    saveScript = Out.saveScript;
    save([path,'\Scrpt_',KRNprms.name],'saveSCRPTipt','saveSCRPTGBL','saveSCRPTIPTGBL','saveSCRPTPATHS','saveScript');    
end
