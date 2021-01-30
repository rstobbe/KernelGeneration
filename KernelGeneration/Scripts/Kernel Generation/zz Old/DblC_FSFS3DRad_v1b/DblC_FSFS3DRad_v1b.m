%==================================================
% (v1b) 
%       - calculate convscaleval for each kernel
%==================================================

function [SCRPTipt,SCRPTGBL,err] = DblC_FSFS3DRad_v1b(SCRPTipt,SCRPTGBL)

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
    Kern_res = FS1_res;
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
res0 = 0.0001;
u0 = (0:res0:(FS1_W/2));   

FTSphereProf1 = 3*(sin(u0*pi*FS1_BW) - pi*FS1_BW*u0.*cos(u0*pi*FS1_BW))./((u0*pi*FS1_BW).^3);
FTSphereProf1(1) = 1;                                              
M = FS1_beta * sqrt(1 - (2*u0/FS1_W).^2);
KBProf1 = besseli(0,M);
KBProf1 = KBProf1/max(KBProf1);
KernProf1 = FTSphereProf1.*KBProf1;
figure(200); hold on;
plot(u0,FTSphereProf1,'b');
plot(u0,KBProf1,'r');
plot(u0,KernProf1,'m');

FTSphereProf2 = 3*(sin(u0*pi*FS2_BW) - pi*FS2_BW*u0.*cos(u0*pi*FS2_BW))./((u0*pi*FS2_BW).^3);
FTSphereProf2(1) = 1;                                              
M = FS2_beta * sqrt(1 - (2*u0/FS2_W).^2);
KBProf2 = besseli(0,M);
KBProf2 = KBProf2/max(KBProf2);
KernProf2 = FTSphereProf2.*KBProf2;
figure(200); hold on;
plot(u0,FTSphereProf2,'b:');
plot(u0,KBProf2,'r:');
plot(u0,KernProf2,'m:');

%--------------------------------------------
% Test
%--------------------------------------------
Status2('busy','Testing',1);
CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5; 
uT = (-(Kern_W/2):CSVres:(Kern_W/2)); 
Rad = ones(length(uT),length(uT),length(uT));
tKern1 = zeros(length(uT),length(uT),length(uT));
tKern2 = zeros(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(Kern_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                tKern1(a,b,c) = lin_interp4(KernProf1,rad,L);
                tKern2(a,b,c) = lin_interp4(KernProf2,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
FwdKern_CSV = sum(tKern1(:))*CSVres^3;
RvsKern_CSV = sum(tKern2(:))*CSVres^3;

interptest = 0;
if interptest == 1
    figure(200); plot(uT,tKern1(:,(N+1)/2,(N+1)/2),'k*');
    %figure(200); plot(Rad(:)*(Kern_W/2),tKern1(:),'k*');
    %figure(200); plot(Rad(:)*(Kern_W/2),tKern2(:),'k*');
end

Status2('busy','Fourier Transform',2);
ZF = 320;
C = ZF/2 + 1;
TF1 = fftshift(fftn(ifftshift(tKern1)));
TF2 = fftshift(fftn(ifftshift(tKern2)));
%TF = TF1;          % testing                                          
TF = TF1.*TF2;
ZFTF = zerofill_isotropic_odd_doubles(ifftshift(TF),ZF);
DblKern = fftshift(ifftn(ZFTF));
clear ZFTF
DblKernProf = DblKern(:,C,C)/abs(DblKern(C,C,C));
clear DblKern
iKern_W = Kern_W + CSVres;
tu0 = (-(iKern_W/2):iKern_W/ZF:(iKern_W/2)-iKern_W/ZF);
iDblKernProf = interp1(tu0,DblKernProf,u0);
figure(200); hold on;
plot(u0,iDblKernProf,'k:');
DblKernProf = real(iDblKernProf);

ZFtKern1 = zerofill_isotropic_odd_doubles(ifftshift(tKern1),ZF);
ZFtKern2 = zerofill_isotropic_odd_doubles(ifftshift(tKern2),ZF);
TF1 = fftshift(fftn(ZFtKern1));
TF2 = fftshift(fftn(ZFtKern2));
clear ZFtKern1
clear ZFtKern2
TF = TF1.*TF2;
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
TF2Prof = squeeze(TF2(C,C,:))/abs(TF2(C,C,C));
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF1
clear TF2
iCSVres = 1/((1/CSVres)+(1/Kern_W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));
figure(300); hold on; title('Image Profile');
plot(f,real(TF1Prof),'b');
plot(f,imag(TF1Prof),'b:');
plot(f,real(TF2Prof),'r');
plot(f,imag(TF2Prof),'r:');
plot(f,real(TFProf),'g');
plot(f,imag(TFProf),'g:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

imageshow = 0;
if imageshow == 1
    tTF = TF(ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6));
    rows = floor(sqrt(length(tTF)));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1;  IMSTRCT.stop = length(tTF); 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 max(abs(tTF(:)))*1.05]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; 
    IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 400; 
    AxialMontage_v2a(tTF,IMSTRCT);
end

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Testing',1);
Rad = ones(length(uT),length(uT),length(uT));
tKern = zeros(length(uT),length(uT),length(uT));
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(Kern_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                tKern(a,b,c) = lin_interp4(DblKernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
DblKern_CSV = sum(tKern(:))*CSVres^3;

Status2('busy','Fourier Transform',2);
ZFtKern = zerofill_isotropic_odd_doubles(ifftshift(tKern),ZF);
TF = fftshift(fftn(ZFtKern));
clear ZFtKern
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF
figure(300); hold on; title('Image Profile');
plot(f,real(TFProf),'k:');
ind1 = find(f<-2,1,'last');
ind2 = find(f>2,1,'first');
xlim([f(ind1) f(ind2)]);
clear tKern
clear tKern1
clear tKern2

Status2('done','',2);
button = questdlg('Continue?','Convolution Kernel Generation','Yes','No','Yes');
if strcmp(button,'No')
    return
end

%==============================================================================

%--------------------------------------------
% Generate Double Convolution Kernel Approximation
%--------------------------------------------
Status2('busy','Generate Double Convolution Kernel Approximation',1);
u = (0:Kern_res:(Kern_W/2)); 
DblKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
Rad = ones((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(Kern_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                DblKern(a,b,c) = lin_interp4(DblKernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(Kern_W/2),DblKern(:),'k*');
end

%--------------------------------------------
% Generate FS1 Kernel
%--------------------------------------------
Status2('busy','Generate FS1 Kernel',1);
u = (0:FS1_res:(FS1_W/2)); 
FS1Kern = zeros((FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,'single');
Rad = ones((FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,(FS1_zW/2)*round(1/FS1_res)+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(FS1_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                FS1Kern(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(FS1_W/2),FS1Kern(:),'k*');
end

%--------------------------------------------
% Generate FS2 Kernel
%--------------------------------------------
Status2('busy','Generate FS2 Kernel',1);
u = (0:FS2_res:(FS2_W/2)); 
FS2Kern = zeros((FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,'single');
Rad = ones((FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,(FS2_zW/2)*round(1/FS2_res)+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(FS2_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                FS2Kern(a,b,c) = lin_interp4(KernProf2,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(FS2_W/2),FS2Kern(:),'k*');
end

%--------------------------------------------
% Save
%--------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'DblC_FSFS3DRad_v1b';
    KRNprms.DblKern.convscaleval = DblKern_CSV;
    KRNprms.DblKern.OKforSS = SS;
    KRNprms.DblKern.W = Kern_W;
    KRNprms.DblKern.iKern = 1/Kern_res;
    KRNprms.DblKern.res = Kern_res;
    KRNprms.DblKern.zW = Kern_zW;
    KRNprms.DblKern.Kern = DblKern;
    
    KRNprms.FwdKern.convscaleval = FwdKern_CSV;
    KRNprms.FwdKern.W = FS1_W;
    KRNprms.FwdKern.BW = FS1_BW;
    KRNprms.FwdKern.beta = FS1_beta;
    KRNprms.FwdKern.iKern = 1/FS1_res;
    KRNprms.FwdKern.res = FS1_res;
    KRNprms.FwdKern.zW = FS1_zW;
    KRNprms.FwdKern.Kern = FS1Kern;
    
    KRNprms.RvsKern.convscaleval = RvsKern_CSV;
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
