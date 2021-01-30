%==================================================
% (v1b) 
%       - backup posterity / update for RWSUI_B9
%==================================================

function [SCRPTipt,SCRPTGBL,err] = SglC_FB3DCart_v1b(SCRPTipt,SCRPTGBL)

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
    err.msg = '1/(Res*TestForSS) not an integer';
    return
end

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
zWtest = 2*(ceil((W*SS-2)/2)+1)/SS;
zW = W + 0.1;
while true
    if zW > zWtest
        break
    end
    zW = zW + 0.1;
end
zW = round(zW*1e9)/1e9;

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
if (ceil((zW/2)*(1/res))+1) > 500
    err.flag = 1;
    err.msg = ['kernel size: ',num2str(ceil((zW/2)*(1/res))+1)];
    return
end

%--------------------------------------------
% Generate Kernel Shape
%--------------------------------------------
Status2('busy','Generate Filtered Box Shape',1);                                                
res0 = 0.0001;
u0 = (-(W/2):res0:(W/2)); 
SincProf = sinc(u0*BW);
M = beta * sqrt(1 - (2*u0/W).^2);
KBProf = besseli(0,M);
KernProf = SincProf.*KBProf;
KernProf = KernProf/max(KernProf);

%--------------------------------------------
% Visualize Kernel
%--------------------------------------------
Status2('busy','Visualize Kernel',1);  
figure(200); hold on;
%plot(u0,SincProf,'m');
%plot(u0,KBProf,'g');
plot(u0,KernProf,'b');
xlim([0 zW/2]); 
ylabel('Kernel k-Space Shape','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value',1);
CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5;                                              
uCSV = (-(W/2):CSVres:(W/2)); 
KernProfCSV = interp1(u0,KernProf,uCSV,'linear','extrap');
tKern = zeros(length(uCSV),length(uCSV),length(uCSV));
N = length(uCSV);
for a = 1:N
    for b = 1:N
        for c = 1:N
            tKern(a,b,c) = KernProfCSV(a)*KernProfCSV(b)*KernProfCSV(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
CSV = sum(tKern(:))*CSVres^3;

%--------------------------------------------
% Plot Image Profile
%--------------------------------------------
Status2('busy','Plot Image Space Profile',1);  
figure(300); hold on;
anlzzf = 64000;
IPres = 0.1;
uIP = (0:IPres:(W/2));
KernProfIP = interp1(u0,KernProf,uIP,'linear','extrap');
ZFKernProfIP = [KernProfIP zeros(1,anlzzf-2*length(KernProfIP)+1) flipdim(KernProfIP(2:length(KernProfIP)),2)];
FTProf = ifftshift(ifft(ZFKernProfIP));
FTProf = FTProf/max(FTProf);
FOV = 1/(2*IPres);                                                        
t = (-FOV:(1/IPres)/anlzzf:FOV-(1/IPres)/anlzzf); 
plot(t,real(FTProf),'r');
plot(t,imag(FTProf),'b'); 
plot([0.5 0.5],[-0.2 1.2],'k:');
plot([-0.5 -0.5],[-0.2 1.2],'k:');
plot([SS/2 SS/2],[-0.2 1.2],'k:');
plot([-SS/2 -SS/2],[-0.2 1.2],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.4) (SS/2+0.4)]);

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
    plot(f,real(TFProf1),'k:');
    plot(f,imag(TFProf1),'k:');
    plot(f,real(TFProf2),'k:');
    plot(f,real(TFProf3),'k:');
    TFProf = real(TFProf1);
    plot([-0.5 -0.5],[min(TFProf) max(TFProf)],'k:');
    plot([0.5 0.5],[min(TFProf) max(TFProf)],'k:');
    ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
    xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

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

%==============================================================

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
u = (0:res:(W/2)); 
KernProfGen = interp1(u0,KernProf,u,'linear','extrap');
figure(200); hold on;
plot(u,KernProfGen,'k*');
Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(KernProfGen);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = KernProfGen(a)*KernProfGen(b)*KernProfGen(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%--------------------------------------------
% Kernal Visualization
%--------------------------------------------
Vis = 0;
if Vis == 1;
    rows = floor(sqrt(N));
    IMSTRCT.type = 'abs'; IMSTRCT.start = 1; IMSTRCT.step = 1; IMSTRCT.stop = N; 
    IMSTRCT.rows = rows; IMSTRCT.lvl = [0 1]; IMSTRCT.docolor = 1; IMSTRCT.ColorMap = 'ColorMap4'; IMSTRCT.SLab = 1; IMSTRCT.figsize = []; IMSTRCT.figno = 900; 
    AxialMontage_v2a(Kern,IMSTRCT);   
end

%--------------------------------------------
% Display
%--------------------------------------------
SCRPTGBL.RWSUI.LocalOutput(1).label = 'ConvScaleVal';
SCRPTGBL.RWSUI.LocalOutput(1).value = num2str(CSV);

%--------------------------------------------
% Output
%--------------------------------------------
KRNprms.type = 'SglC_FB3DCart_v1b';
KRNprms.prec = 'singles';

KRNprms.convscaleval = CSV;
KRNprms.DesforSS = SS;
KRNprms.iKern = 1/res;
KRNprms.W = W;
KRNprms.BW = BW;
KRNprms.beta = beta;
KRNprms.res = res;
KRNprms.zW = zW;
KRNprms.Kern = Kern;

name = ['Kern_FBCw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3)];
KRNprms.name = regexprep(name,'\.', 'p');

SCRPTGBL.RWSUI.SaveGlobal = 'yes';
SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\',KRNprms.name];
SCRPTGBL.RWSUI.SaveVariables = {KRNprms};
SCRPTGBL.RWSUI.SaveVariableNames = {'KRNprms'};

