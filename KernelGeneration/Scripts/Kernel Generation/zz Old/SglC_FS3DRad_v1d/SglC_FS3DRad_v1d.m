%==================================================
% (v1d) 
%   - Update generation to '<='
%   - search for possible SS
%==================================================

function [SCRPTipt,SCRPTGBL,err] = SglC_FS3DRad_v1d(SCRPTipt,SCRPTGBL)

Status('busy','Create Convolution Kernel');
Status2('done','',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Load Input
%---------------------------------------------
W = str2double(SCRPTGBL.CurrentScript.Width);
beta = str2double(SCRPTGBL.CurrentScript.Beta);
BW = str2double(SCRPTGBL.CurrentScript.BW);
res = str2double(SCRPTGBL.CurrentScript.Res);
SS = str2double(SCRPTGBL.CurrentScript.DesForSS); 

%---------------------------------------------
% Tests
%---------------------------------------------
if rem((W/2),res)
    err.flag = 1;
    err.msg = 'Width/2 not a multiple of Res';
    return
end
if rem(round(1e9*(1/(res*SS)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(Res*DesForSS) not an integer';
    return
end

%---------------------------------------------
% Find all SS values that are compatible with res
%---------------------------------------------
m = 1;
for n = 1:0.0005:6
    test = round(1e9/(res*n))/1e9;
    if not(rem(test,1)) && not(rem(W*n,1))
        posSS(m) = n;
        m = m+1;
    end    
end

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
for n = 1:length(posSS)
    zWtest = 2*(ceil((W*posSS(n)-2)/2)+1)/posSS(n);
    zW = W + 0.1;
    while true
        if zW > zWtest
            break
        end
        zW = zW + 0.1;
    end
    zWarr(n) = round(zW*1e9)/1e9;
end
zW = max(zWarr);

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
kernelsize = ceil((zW/2)*(1/res))+1;
button = questdlg(['Kernel Dimensions: ',num2str(kernelsize),'.  Continue?']);
if strcmp(button,'No')
    err.flag = 4;
    err.msg = '';
    return
end
if kernelsize > 550
    err.flag = 1;
    err.msg = 'excessive kernlsize';
    return
end

%--------------------------------------------
% Generate fSphere profile
%--------------------------------------------
Status2('busy','Generate Filtered Sphere Sphape',1); 
res0 = 0.0001;
u0 = (0:res0:(W/2));   
FTSphereProf = 3*(sin(u0*pi*BW) - pi*BW*u0.*cos(u0*pi*BW))./((u0*pi*BW).^3);
FTSphereProf(1) = 1;

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

%--------------------------------------------
% Create Combined Profile
%--------------------------------------------
KernProf = FTSphereProf.*KBProf;

%--------------------------------------------
% Visualize Kernel
%--------------------------------------------
figure(200); hold on;
%plot(u0,FTSphereProf,'m');
%plot(u0,KBProf,'g');
plot(u0,KernProf,'b');
%plot([0 zW/2],[0 0],'k:');
%xlim([0 zW/2]); 
plot([0 W/2],[0 0],'k:');
xlim([0 W/2]); 
ylim([-0.2 1.1])
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');
OutStyle(1);

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value',1);
CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5;                                              
uT = (-(W/2):CSVres:(W/2)); 
tKern = zeros(length(uT),length(uT),length(uT));
Rad = ones(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                tKern(a,b,c) = lin_interp4(KernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
CSV = sum(tKern(:))*CSVres^3

%--------------------------------------------
% Interpolation Test
%--------------------------------------------
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),tKern(:),'k*');
end

%--------------------------------------------
% Image Profile Test
%--------------------------------------------
testIP = 1;
if testIP == 1
    Status2('busy','Test Image Space Profile',1); 
    ZF = 300;
    tKern = zerofill_isotropic_odd_doubles(ifftshift(tKern),ZF);
    TF = fftshift(fftn(tKern))*CSVres^3/CSV;
    C = ZF/2 + 1;
    TFProf1 = squeeze(TF(C,C,:));
    TFProf2 = squeeze(TF(:,C,C));
    TFProf3 = squeeze(TF(C,:,C));
    f = (-1/(2*CSVres):1/(ZF*CSVres):1/(2*CSVres) - 1/(ZF*CSVres));
    figure(300); hold on;
    plot(f,real(TFProf1),'r');
    plot(f,imag(TFProf1),'b');
    plot(f,real(TFProf2),'r');
    plot(f,real(TFProf3),'r');
    TFProf = real(TFProf1);
    plot([-0.5 -0.5],[min(TFProf) max(TFProf)],'k:');
    plot([0.5 0.5],[min(TFProf) max(TFProf)],'k:');
    ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
    xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
    xlim([-1.2 1.2]); ylim([-0.1 1.1]);
    OutStyle(1);
    
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
Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
Rad = ones(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(u);
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                Kern(a,b,c) = lin_interp4(KernProf,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%--------------------------------------------
% Interpolation Test
%--------------------------------------------
interptest = 0;
if interptest == 1
    figure(200); plot(Rad(:)*(W/2),Kern(:),'k*');
end

%--------------------------------------------
% Kernal Visualization
%--------------------------------------------
Vis = 1;
if Vis == 1;
    rows = floor(sqrt(N));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1; IMSTRCT.stop = N; 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 1]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 900; 
    AxialMontage_v2a(Kern,IMSTRCT);   
end

%--------------------------------------------
% Output
%--------------------------------------------
KRNprms.type = 'SglC_FS3DRad_v1d';
KRNprms.prec = 'singles';
KRNprms.posSS = posSS;

KRNprms.convscaleval = CSV;
KRNprms.DesforSS = SS;
KRNprms.iKern = 1/res;
KRNprms.W = W;
KRNprms.BW = BW;
KRNprms.beta = beta;
KRNprms.res = res;
KRNprms.zW = zW;
KRNprms.Kern = Kern;

name = ['Kern_FSRw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3)];
KRNprms.name = regexprep(name,'\.', 'p');

%--------------------------------------------
% Return
%--------------------------------------------
name = inputdlg('Name Kern:','Name Kernel',1,{KRNprms.name});
if isempty(name)
    SCRPTGBL.RWSUI.SaveGlobal = 'no';
    return
end
SCRPTipt(strcmp('Kern_Name',{SCRPTipt.labelstr})).entrystr = cell2mat(name);

SCRPTGBL.RWSUI.SaveVariables = {KRNprms};
SCRPTGBL.RWSUI.SaveVariableNames = {'KRNprms'};

SCRPTGBL.RWSUI.SaveGlobal = 'yes';
SCRPTGBL.RWSUI.SaveGlobalNames = name;
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\',KRNprms.name];
SCRPTGBL.RWSUI.SaveScriptName = KRNprms.name;

Status('done','');
Status2('done','',2);
Status2('done','',3);
