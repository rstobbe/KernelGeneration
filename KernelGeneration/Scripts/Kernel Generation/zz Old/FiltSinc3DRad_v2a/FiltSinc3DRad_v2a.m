%==================================================
% (v2a) 
%   
%==================================================

function [SCRPTipt,SCRPTGBL,err] = FiltSphere_v1a(SCRPTipt,SCRPTGBL)

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
% Generate fSphere profile
%--------------------------------------------
res0 = 0.0001;
u0 = (0:res0:(W/2));   
FTSphereProf = 3*(sin(u0*pi*BW) - pi*BW*u0.*cos(u0*pi*BW))./((u0*pi*BW).^3);
FTSphereProf(1) = 1;
figure(200); hold on;
plot(u0,FTSphereProf,'k');

%--------------------------------------------
% fSphere profile test
%--------------------------------------------
test = 0;
if test == 1
    Status2('busy','Generate FT-Sphere Profile',1);
    Status2('busy','Generate Sphere',2);
    N = 320;
    Sphere_diam = 32;
    R = Sphere_diam/2;
    Sphere = zeros(N,N,N); 
    C = (N/2)+1;
    Tot = 0;
    for x = C-R:C+R
        for y = C-R:C+R
            for z = C-R:C+R
                r = (sqrt((x-C)^2 + (y-C)^2 + (z-C)^2))/R;
                if r <= 1
                    Sphere(x,y,z) = 1;
                    Tot = Tot + 1;
                end
            end
        end
        Status2('busy',num2str(x),3);    
    end
    Status2('done','',3);
    compV = pi*(4/3)*R.^3;
    Dif = compV/Tot;
    if abs((1-Dif)*100) > 2
        err.flag = 1;
        err.msg = 'Sphere Creation More than 2% in error (increase matrix dimension)';
        return
    end
    Status2('busy','Fourier Transform',2);
    FTsphere = fftshift(ifftn(ifftshift(Sphere)));
    FTsphere = FTsphere/max(abs(FTsphere(:)));
    figure(100); hold on;
    t = (-Sphere_diam/2:Sphere_diam/N:Sphere_diam/2-Sphere_diam/N);
    plot(t,squeeze(real(FTsphere(C,C,:))),'b');
    plot(t,squeeze(imag(FTsphere(C,C,:))),'b:');
    tFTSphereProf = squeeze(real(FTsphere(C,C,:)));
    res0 = 0.01;
    tu0 = (0:res0:(W/2)*BW);   
    tFTSphereProf = interp1(t,tFTSphereProf,tu0);
    figure(100); hold on;
    plot(tu0,tFTSphereProf,'k*');
    tu0 = tu0/BW;
    figure(200); hold on;
    plot(tu0,tFTSphereProf,'b:');
    Status2('done','',2);
end

%--------------------------------------------
% Generate Kaiser Profile
%--------------------------------------------
Status2('busy','Generate Kaiser Profile',1);                                                
M = beta * sqrt(1 - (2*u0/W).^2);
KBProf = besseli(0,M);
KBProf = KBProf/max(KBProf);
figure(200); hold on; plot(u0,KBProf,'r');

%--------------------------------------------
% Create Combined Profile
%--------------------------------------------
KernProf = FTSphereProf.*KBProf;
figure(200); hold on; plot(u0,KernProf,'g');

%--------------------------------------------
% Test and Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Testing',1);
CSVres = 0.25;                                             
uT = (-(W/2):CSVres:(W/2)); 
tKern = zeros(length(uT),length(uT),length(uT),'single');
Rad = ones(length(uT),length(uT),length(uT),'single');
N = length(uT);
L = length(u0);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad < 1
                Rad(a,b,c) = rad;
                tKern(a,b,c) = lin_interp4(KernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
ConvScaleVal = sum(tKern(:))*CSVres^3;

interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),tKern(:),'k*');
end

ZF = 300;
tKern = zerofill_isotropic_odd_doubles(ifftshift(tKern),ZF);
TF = fftshift(fftn(tKern));
C = ZF/2 + 1;
TFProf1 = squeeze(TF(C,C,:));
TFProf2 = squeeze(TF(:,C,C));
TFProf3 = squeeze(TF(C,:,C));
f = (-1/(2*CSVres):1/(ZF*CSVres):1/(2*CSVres) - 1/(ZF*CSVres));
figure(300); hold on; title('Image Profile');
plot(f,real(TFProf1));
plot(f,imag(TFProf1));
plot(f,real(TFProf2));
plot(f,real(TFProf3));
TFProf = real(TFProf1);
plot([-0.5 -0.5],[min(TFProf) max(TFProf)],'k:');
plot([0.5 0.5],[min(TFProf) max(TFProf)],'k:');

imageshow = 1;
if imageshow == 1
    TF = TF(ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6),ZF/2-floor(ZF/6):ZF/2+ceil(ZF/6));
    rows = floor(sqrt(length(TF)));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1;  IMSTRCT.stop = length(TF); 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 max(TFProf)*1.05]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 400; 
    AxialMontage_v2a(TF,IMSTRCT);
end

button = questdlg('Continue?','Convolution Kernel Generation','Yes','No','Yes');
if strcmp(button,'No')
    return
end

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
                Kern(a,b,c) = lin_interp4(KernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),Kern(:),'k*');
end

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
    KRNprms.type = 'FiltSphere';
    KRNprms.iKern = 1/res;
    KRNprms.W = W;
    KRNprms.res = res;
    KRNprms.zW = zW;
    KRNprms.beta = beta;
    KRNprms.BW = BW;
    KRNprms.convscaleval = ConvScaleVal;
    KRNprms.SubSamp = SS;
    
    name = ['FSph_W',num2str(W,2),'_R',num2str(res,3),'_BW',num2str(BW,3),'_B',num2str(beta,3),'_zW',num2str(zW,3)];
    KRNprms.name = regexprep(name,'\.', 'p');

    path = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\FiltSphere\',KRNprms.name];
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
