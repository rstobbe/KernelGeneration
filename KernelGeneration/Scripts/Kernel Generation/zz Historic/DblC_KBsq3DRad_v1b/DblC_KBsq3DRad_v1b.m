%==================================================
% (v1a) 
%     
%==================================================

function [SCRPTipt,SCRPTGBL,err] = DblC_KBsq3DRad_v1a(SCRPTipt,SCRPTGBL)

SCRPTGBL.TextBox = '';
SCRPTGBL.Figs = [];
SCRPTGBL.Data = [];

err.flag = 0;
err.msg = '';

Kern_W = str2double(SCRPTipt(strcmp('Width',{SCRPTipt.labelstr})).entrystr);
Kern_res = str2double(SCRPTipt(strcmp('Res',{SCRPTipt.labelstr})).entrystr); 
Kern_beta = str2double(SCRPTipt(strcmp('Beta',{SCRPTipt.labelstr})).entrystr); 
Kern_zWadd = str2double(SCRPTipt(strcmp('zWadd',{SCRPTipt.labelstr})).entrystr);
SS = str2double(SCRPTipt(strcmp('TestForSS',{SCRPTipt.labelstr})).entrystr); 
Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr;
if iscell(Vis)
    Vis = SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entrystr{SCRPTipt(strcmp('Visuals',{SCRPTipt.labelstr})).entryvalue};
end

Kern_zW = Kern_W + Kern_zWadd; 
KRNprms.prec = 'singles';

if rem((Kern_W/2),Kern_res)
    err.flag = 1;
    err.msg = 'Width/2 not a multiple of Res';
    return
end
if rem(1,Kern_res) 
    err.flag = 1;
    err.msg = '1/Res is not an integer';
    return
end
if rem(round(1e9*(1/(Kern_res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(Res*TestForSS) not an integer';
    return
end
zWtest = 2*(ceil((Kern_W*SS-2)/2)+1)/SS;
if zWtest > Kern_zW
    err.flag = 1;
    err.msg = 'zWadd Too Small';
    return
end
if Kern_zW/Kern_res > 1100
    err.flag = 1;
    err.msg = ['Kern kernel size: ',num2str(Kern_zW/Kern_res)];
    return
end

%--------------------------------------------
% Generate Kaiser profile
%--------------------------------------------
Status2('busy','Generate FT-Sphere Profile',1);
res0 = 0.00001;
u0 = (0:res0:(Kern_W/2));                                               
M = Kern_beta * sqrt(1 - (2*u0/Kern_W).^2);
KBProf = besseli(0,M);
KBProf = KBProf/KBProf(1);
figure(200); hold on;
plot(u0,KBProf,'k');

%--------------------------------------------
% Test
%--------------------------------------------
Status2('busy','Testing',1);
CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5; 
uT = (-(Kern_W/2):CSVres:(Kern_W/2)); 
Rad = ones(length(uT),length(uT),length(uT));
tKern = zeros(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(Kern_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                tKern(a,b,c) = lin_interp4(KBProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
FwdKern_CSV = sum(tKern(:))*CSVres^3;
RvsKern_CSV = FwdKern_CSV;

interptest = 0;
if interptest == 1
    figure(200); plot(uT,tKern(:,(N+1)/2,(N+1)/2),'k*');
    %figure(200); plot(Rad(:)*(Kern_W/2),tKern(:),'k*');
end

Status2('busy','Fourier Transform',2);
ZF = 320;
C = ZF/2 + 1;
TF1 = fftshift(fftn(ifftshift(tKern)));                                       
TF = TF1.^2;
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

ZFtKern = zerofill_isotropic_odd_doubles(ifftshift(tKern),ZF);
TF1 = fftshift(fftn(ZFtKern));
clear ZFtKern
TF = TF1.^2;
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF1
iCSVres = 1/((1/CSVres)+(1/Kern_W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));
figure(300); hold on; title('Image Profile');
plot(f,real(TF1Prof),'b');
plot(f,imag(TF1Prof),'b:');
plot(f,real(TFProf),'g');
plot(f,imag(TFProf),'g:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');

imageshow = 1;
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
interptest = 1;
if interptest == 1
    figure(200); plot(Rad(:)*(Kern_W/2),DblKern(:),'k*');
end

%--------------------------------------------
% Generate KB Kernel
%--------------------------------------------
Status2('busy','Generate KB Kernel',1);
u = (0:Kern_res:(Kern_W/2)); 
KBKern = zeros((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
Rad = ones((Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,(Kern_zW/2)*round(1/Kern_res)+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(Kern_W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                KBKern(a,b,c) = lin_interp4(KBProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 1;
if interptest == 1
    figure(200); plot(Rad(:)*(Kern_W/2),KBKern(:),'k*');
end

%--------------------------------------------
% Save
%--------------------------------------------
button = questdlg('Save?');
if strcmp(button,'Yes')
    KRNprms.type = 'DblC_KBsq3DRad_v1a';
    KRNprms.DblKern.convscaleval = DblKern_CSV;
    KRNprms.DblKern.OKforSS = SS;
    KRNprms.DblKern.W = Kern_W;
    KRNprms.DblKern.beta = Kern_beta;
    KRNprms.DblKern.iKern = 1/Kern_res;
    KRNprms.DblKern.res = Kern_res;
    KRNprms.DblKern.zW = Kern_zW;
    KRNprms.DblKern.Kern = DblKern;
    
    KRNprms.FwdKern.convscaleval = FwdKern_CSV;
    KRNprms.FwdKern.W = Kern_W;
    KRNprms.FwdKern.beta = Kern_beta;
    KRNprms.FwdKern.iKern = 1/Kern_res;
    KRNprms.FwdKern.res = Kern_res;
    KRNprms.FwdKern.zW = Kern_zW;
    KRNprms.FwdKern.Kern = KBKern;
    
    KRNprms.RvsKern.convscaleval = RvsKern_CSV;
    KRNprms.RvsKern.W = Kern_W;
    KRNprms.RvsKern.beta = Kern_beta;
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
