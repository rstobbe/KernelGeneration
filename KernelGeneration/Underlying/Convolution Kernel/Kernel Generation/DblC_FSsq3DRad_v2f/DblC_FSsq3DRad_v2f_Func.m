%=========================================================
% 
%=========================================================

function [KDES,err] = DblC_FSsq3DRad_v2f_Func(KDES,INPUT)

Status2('done','Create Kernel',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Get Input
%---------------------------------------------
W = KDES.W;
BW = KDES.BW;
res = KDES.res;
beta = KDES.beta;
SS = KDES.DesforSS;
clear INPUT;

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
if rem(W*SS,1)
    err.flag = 1;
    err.msg = 'Width*DesForSS not an integer';
    return
end

%---------------------------------------------
% Find all SS values that are compatible with res
%---------------------------------------------
m = 1;
for n = 1:0.0005:10
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
if kernelsize > 1200
    err.flag = 1;
    err.msg = 'excessive kernlsize';
    return
end

%--------------------------------------------
% Generate FilteredSinc Kernel Profile
%--------------------------------------------
res0 = 0.0001;
u0 = (0:res0:(W/2));   
FTSphereProf1 = 3*(sin(u0*pi*BW) - pi*BW*u0.*cos(u0*pi*BW))./((u0*pi*BW).^3);
FTSphereProf1(1) = 1;                                              
M = beta * sqrt(1 - (2*u0/W).^2);
KBProf1 = besseli(0,M);
KBProf1 = KBProf1/max(KBProf1);
KernProf1 = FTSphereProf1.*KBProf1;
figure(200); clf; hold on;
plot(u0,KernProf1,'b');
plot(u0,KBProf1,'g');
plot(u0,FTSphereProf1,'c');
box on;
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

%--------------------------------------------
% Convolution Scale Value
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value of Single Kernel',1);
CSVres = 0.05;                                      % size of this value doesn't really matter -> same result just more accurate calculation
uT = (-(W/2):CSVres:(W/2)); 
Rad = ones(length(uT),length(uT),length(uT));
tKern1 = zeros(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                tKern1(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
CSV = sum(tKern1(:))*CSVres^3;

%--------------------------------------------
% Generate Double Kernel Profile (for testing)
%--------------------------------------------
Status2('busy','Determine Double Kern k-Space Shape',1);
%ZF = 320;
ZF = 600;
C = ZF/2 + 1;
TF1 = fftshift(fftn(ifftshift(tKern1)));                                     
TF = TF1.^2;
ZFTF = zerofill_isotropic_odd_doubles(ifftshift(TF),ZF);
DblKern = fftshift(ifftn(ZFTF));
clear ZFTF
DblKernProf = DblKern(:,C,C)/abs(DblKern(C,C,C));
clear DblKern
iW = W + CSVres;
tu0 = (-(iW/2):iW/ZF:(iW/2)-iW/ZF);
iDblKernProf = interp1(tu0,DblKernProf,u0);
hFig = figure(200); hold on;
plot(u0,iDblKernProf,'r');
KDES.Figure(1).Name = 'KernShape';
KDES.Figure(1).Type = 'Graph';
KDES.Figure(1).hFig = hFig;
KDES.Figure(1).hAx = gca;

%--------------------------------------------
% Image Space Profiles
%--------------------------------------------
Status2('busy','Determine Double Kern Image Shape',1);
ZFtKern1 = zerofill_isotropic_odd_doubles(ifftshift(tKern1),ZF);
TF1 = fftshift(fftn(ZFtKern1));
clear ZFtKern1
TF = TF1.^2;
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
TFProf = squeeze(TF(C,C,:))/abs(TF(C,C,C));
clear TF1
iCSVres = 1/((1/CSVres)+(1/W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));
hFig = figure(300); clf; hold on;
plot(f,real(TF1Prof),'b');
plot(f,imag(TF1Prof),'b:');
plot(f,real(TFProf),'r','linewidth',2);
plot(f,imag(TFProf),'r:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
%plot([-SS+0.5 -SS+0.5],[0 1],'k:'); plot([SS-0.5 SS-0.5],[0 1],'k:');
xlim([-1.5 1.5]);
box on;
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
KDES.Figure(2).Name = 'KernImSpace';
KDES.Figure(2).Type = 'Graph';
KDES.Figure(2).hFig = hFig;
KDES.Figure(2).hAx = gca;

%--------------------------------------------
% Continue
%--------------------------------------------
Status2('done','',2);
button = questdlg('Continue?','Convolution Kernel Generation','Yes','No','Yes');
if strcmp(button,'No')
    err.flag = 4;
    return
end

%==============================================================================

%--------------------------------------------
% Generate FilteredSinc Kernel
%--------------------------------------------
Status2('busy','Generate FilteredSinc Kernel',1);
u = (0:res:(W/2)); 
FS1Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
Rad = ones(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(u);
Status2('busy','Build 3D Kernel',2);  
parfor a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Rad(a,b,c) = rad;
                FS1Kern(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    %Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%----------------------------------------------------
% Panel Output
%----------------------------------------------------
Panel(1,:) = {'',KDES.method,'Output'};
Panel(2,:) = {'Width',W,'Output'};
Panel(3,:) = {'Beta',beta,'Output'};
Panel(4,:) = {'Res',res,'Output'};
Panel(5,:) = {'DesForSS',KDES.DesforSS,'Output'};
Panel(6,:) = {'ConvScalVal',CSV,'Output'};
KDES.Panel = Panel;
KDES.PanelOutput = cell2struct(Panel,{'label','value','type'},2);

%--------------------------------------------
% Output
%--------------------------------------------
KDES.posSS = posSS;

KDES.DblKern.OKforSS = SS;
KDES.DblKern.W = W;
KDES.DblKern.iKern = 1/res;
KDES.DblKern.res = res;
KDES.DblKern.zW = zW;

KDES.FwdKern.convscaleval = CSV;
KDES.FwdKern.W = W;
KDES.FwdKern.BW = BW;
KDES.FwdKern.beta = beta;
KDES.FwdKern.iKern = 1/res;
KDES.FwdKern.res = res;
KDES.FwdKern.zW = zW;
KDES.FwdKern.Kern = FS1Kern;

KDES.RvsKern.convscaleval = CSV;
KDES.RvsKern.W = W;
KDES.RvsKern.BW = BW;
KDES.RvsKern.beta = beta;
KDES.RvsKern.iKern = 1/res;
KDES.RvsKern.res = res;
KDES.RvsKern.zW = zW;
  
name = ['Kern_FSsqRw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3),'ss',num2str(SS,3),'r',num2str(res,3)];
KDES.name = regexprep(name,'\.', 'p');

Status2('done','',2);
Status2('done','',3);
