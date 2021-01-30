%==================================================
% (v1a) 
%   
%==================================================

function [SCRPTipt,SCRPTGBL,err] = FiltSinc3DRad_v1a(SCRPTipt,SCRPTGBL)

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
% Generate FS and Kaiser Kernel Shapes for Testing
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);                                                
res0 = 0.001;
u0 = (-W/2:res0:W/2-res0);   
Sinc_1D = sinc(u0*BW);
M = beta * sqrt(1 - (2*u0/W).^2);
KB_1D = besseli(0,M);
Kern_1D = Sinc_1D.*KB_1D;
Kern_1D = Kern_1D/max(Kern_1D);
figure(100); hold on; plot(u0,Kern_1D,'b');

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
% Generate FS and Kaiser Kernel Shapes for Building
%--------------------------------------------
Status2('busy','Generate Kernel Shape',1);                                                
res0 = 0.00001;
u0 = (0:res0:W/2);   
Sinc_1D = sinc(u0*BW);
M = beta * sqrt(1 - (2*u0/W).^2);
KB_1D = besseli(0,M);
Kern_1D = Sinc_1D.*KB_1D;
Kern_1D = Kern_1D/max(Kern_1D);

%--------------------------------------------
% Test ConvScaleVal
%--------------------------------------------
Status2('busy','ConvScaleVal',1);
CSVres = 0.1;
u = (-(W/2):CSVres:(W/2));                                                 
Kern = zeros(length(u),length(u),length(u),'single');
Rad = ones(length(u),length(u),length(u),'single');
N = length(u);
L = length(u0);
pts = 0;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                Kern(a,b,c) = lin_interp4(Kern_1D,rad,L);
                %V = lin_interp4(Kern_1D,rad,L);
                pts = pts+1;
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
figure(100); plot(Rad(:)*(W/2),Kern(:),'k*');
ConvScaleVal = sum(Kern(:))*CSVres^3

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
u = (0:res:(W/2)); 
Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
Rad = ones((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
N = length(u);
L = length(u0);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                Kern(a,b,c) = lin_interp4(Kern_1D,rad,L);
                %V = lin_interp4(Kern_1D,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
%figure(100); plot(Rad(:)*(W/2),Kern(:),'r*');

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
if strcmp(Vis,'On')    
    rows = floor(sqrt(N));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1; IMSTRCT.stop = N; 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 1]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 900; 
    AxialMontage_v2a(Kern,IMSTRCT);   
end

%--------------------------------------------
% Save
%--------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'FiltSincRad';
    KRNprms.iKern = 1/res;
    KRNprms.W = W;
    KRNprms.res = res;
    KRNprms.zW = zW;
    KRNprms.beta = beta;
    KRNprms.BW = BW;
    KRNprms.convscaleval = ConvScaleVal;
    KRNprms.SubSamp = SS;
    
    name = ['FSR_W',num2str(W,2),'_R',num2str(res,3),'_BW',num2str(BW,3),'_B',num2str(beta,3),'_zW',num2str(zW,3)];
    KRNprms.name = regexprep(name,'\.', 'p');

    path = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\FiltSincRad\',KRNprms.name];
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
