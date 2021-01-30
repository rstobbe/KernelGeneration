%==================================================
% (v1a) 
%  
%==================================================

function [SCRPTipt,SCRPTGBL,err] = FiltBox3DCart_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

W = str2double(SCRPTipt(strcmp('KernWidth',{SCRPTipt.labelstr})).entrystr);
res = str2double(SCRPTipt(strcmp('KernRes',{SCRPTipt.labelstr})).entrystr); 
BW = str2double(SCRPTipt(strcmp('KernBW',{SCRPTipt.labelstr})).entrystr); 
beta = str2double(SCRPTipt(strcmp('KernBeta',{SCRPTipt.labelstr})).entrystr); 
zW_add = str2double(SCRPTipt(strcmp('zW_add',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('SubSamp',{SCRPTipt.labelstr})).entrystr); 
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
% Generate FS1 and FS2 Kernel Shapes
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);                                                
res0 = 0.001;
u1 = (-W/2:res0:W/2-res0);   
Sinc_1D = sinc(u1*BW);
M = beta * sqrt(1 - (2*u1/W).^2);
KB_1D = besseli(0,M);
Kern_1D = Sinc_1D.*KB_1D;
Kern_1D = Kern_1D/max(Kern_1D);
figure(100); hold on; plot(u1,Kern_1D,'b');

%--------------------------------------------
% Plot Kernel Frequency Profiles
%--------------------------------------------
FTKern_1D = fft(fftshift(Kern_1D));
figure(101); hold on; 
f = (-1/(2*res0):1/W:1/(2*res0)-1/W);
plot(f,real(fftshift(FTKern_1D)/max(abs(FTKern_1D(:)))),'b');
plot(f,imag(fftshift(FTKern_1D)/max(abs(FTKern_1D(:)))),'b');
plot([-0.5 -0.5],[0 1],'k:'); plot([0.5 0.5],[0 1],'k:'); xlim([-4 4]);
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

%--------------------------------------------
% Test ConvScaleVal
%--------------------------------------------
Status2('busy','ConvScaleVal',1);
CSVres = 0.1;
u = (-(W/2):CSVres:(W/2));                                                 
CSVtest = interp1(u1,Kern_1D,u,'linear','extrap');
figure(100); plot(u,CSVtest,'r*');
N = length(CSVtest);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = CSVtest(a)*CSVtest(b)*CSVtest(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
Kern = Kern / max(max(max(Kern)));
ConvScaleVal = sum(Kern(:))*CSVres^3

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
u = (0:res:(W/2));                                                 
Kern_1D = interp1(u1,Kern_1D,u,'linear','extrap');
figure(100); plot(u,Kern_1D,'k*');
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
% Kernel Vis
%--------------------------------------------
if strcmp(Vis,'On')
    figure(1100); hold on; plot((1:length(Kern)),squeeze(Kern(1,1,:))); xlim([0 length(Kern)]); xlabel('kernel size');
    figure(1101); hold on; plot((0:length(Kern)-1)*res,squeeze(Kern(1,1,:))); xlim([0 zW/2]); xlabel('kernel width');
    
    KernProf = [Kern_1D zeros(1,2000) flipdim(Kern_1D(2:length(Kern_1D)),2)];
    OS = length(KernProf)/(2*length(Kern_1D));
    T1 = ifftshift(ifft(KernProf));
    T1 = T1/max(T1);
    
    bot = (length(T1)+1)/2 +1 - (W/2*OS);
    top = (length(T1)+1)/2 +1 + (W/2*OS);
    figure(1102); hold on; plot((1:length(T1)),T1,'r'); plot([bot bot],[0 1],'k:'); plot([top top],[0 1],'k:'); xlim([0.98*bot top/0.98]);
end

%--------------------------------------------
% Save
%--------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'FiltSinc_v2';
    KRNprms.iKern = 1/res;
    KRNprms.W = W;
    KRNprms.res = res;
    KRNprms.zW = zW;
    KRNprms.beta = beta;
    KRNprms.BW = BW;
    KRNprms.convscaleval = ConvScaleVal;
    KRNprms.SubSamp = SS;
    
    name = ['FSv2_W',num2str(W,2),'_R',num2str(res,3),'_BW',num2str(BW,3),'_B',num2str(beta,3),'_zW',num2str(zW,3)];
    KRNprms.name = regexprep(name,'\.', 'p');

    path = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\FiltSinc_v2\',KRNprms.name];
    mkdir(path);
    save([path,'\',KRNprms.name],'Kern','KRNprms');
    [Out,err] = ExternalSave(SCRPTGBL.scrptnum,SCRPTGBL.scrpt,path);
    SCRPTipt = Out.saveSCRPTipt;
    saveSCRPTipt = Out.saveSCRPTipt;
    saveSCRPTGBL = Out.saveSCRPTGBL;
    saveSCRPTIPTGBL = Out.saveSCRPTIPTGBL;
    saveSCRPTPATHS = Out.saveSCRPTPATHS;
    saveScript = Out.saveScript;
    save([path,'\Scrpt_',KRNprms.name],'saveSCRPTipt','saveSCRPTGBL','saveSCRPTIPTGBL','saveSCRPTPATHS','saveScript');    
end
