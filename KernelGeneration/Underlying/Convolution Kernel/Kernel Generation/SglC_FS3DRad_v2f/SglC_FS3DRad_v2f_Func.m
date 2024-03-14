%=========================================================
% 
%=========================================================

function [KDES,err] = SglC_FS3DRad_v2f_Func(KDES,INPUT)

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
hFig = figure(200); clf; 
subplot(1,2,1); hold on;
plot(u0,KernProf1,'k');
plot(u0,FTSphereProf1,'c');
plot(u0,KBProf1,'g');
box on;
legend('KernProf','FTSphereProf','KBProf');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

hFig.Units = 'inches';
hFig.Position = [6 5 14 5];
KDES.Figure(1).Name = 'KernImSpace';
KDES.Figure(1).Type = 'Graph';
KDES.Figure(1).hFig = hFig;
KDES.Figure(1).hAx = gca;

%--------------------------------------------
% Convolution Scale Value
%--------------------------------------------
Status2('busy','Calculate Convolution Scale Value of Single Kernel',1);
CSVres = 0.05;                                      % size of this value doesn't really matter -> same result just more accurate calculation
uT = (-(W/2):CSVres:(W/2)); 
tKern1 = zeros(length(uT),length(uT),length(uT));
N = length(uT);
L = length(u0)-1;
for a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((uT(a).^2 + uT(b).^2 + uT(c).^2).^(0.5))/(W/2);
            if rad <= 1
                tKern1(a,b,c) = lin_interp4(KernProf1,rad,L);
            end
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
CSV = sum(tKern1(:))*CSVres^3;

%--------------------------------------------
% Image Space Profiles
%--------------------------------------------
ZF = 600;
C = ZF/2 + 1;
Status2('busy','Determine Kern Image Shape',1);
ZFtKern1 = zerofill_isotropic_odd_doubles(ifftshift(tKern1),ZF);
TF1 = fftshift(fftn(ZFtKern1));
clear ZFtKern1
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
clear TF1
iCSVres = 1/((1/CSVres)+(1/W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));
subplot(1,2,2); hold on;
plot(f,real(TF1Prof),'k');
plot(f,imag(TF1Prof),'k:');
plot([-0.5 -0.5],[0 1],'k:');
plot([0.5 0.5],[0 1],'k:');
plot([-0.5*SS -0.5*SS],[0 1],'k:'); plot([0.5*SS 0.5*SS],[0 1],'k:');
xlim([-1.5 1.5]);
box on;
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

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
Kern = zeros(ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,ceil((zW/2)*(1/res))+1,'single');
N = length(u);
Status2('busy','Build 3D Kernel',2);  
parfor a = 1:N
    for b = 1:N
        for c = 1:N
            rad = ((u(a).^2 + u(b).^2 + u(c).^2).^(0.5))/(W/2);
            if rad <= 1
                Kern(a,b,c) = lin_interp4(KernProf1,rad,L);
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
Panel(3,:) = {'BW',BW,'Output'};
Panel(4,:) = {'Beta',beta,'Output'};
Panel(5,:) = {'Res',res,'Output'};
Panel(6,:) = {'DesForSS',KDES.DesforSS,'Output'};
Panel(7,:) = {'ConvScalVal',CSV,'Output'};
KDES.Panel = Panel;
KDES.PanelOutput = cell2struct(Panel,{'label','value','type'},2);

%--------------------------------------------
% Output
%--------------------------------------------
KDES.convscaleval = CSV;
KDES.DesforSS = SS;
KDES.posSS = posSS;
KDES.iKern = 1/res;
KDES.zW = zW;
KDES.Kern = Kern;
if SS == 1.25
    KDES.PossibleZeroFill = (50:10:1000);
elseif SS == 1.28
    KDES.PossibleZeroFill = (64:64:1000);
elseif SS == 1.6
    KDES.PossibleZeroFill = (64:16:1000);
elseif SS == 2.0
    KDES.PossibleZeroFill = (64:8:1000);
elseif SS == 2.5
    KDES.PossibleZeroFill = (50:10:1000);
elseif SS == 4.0
    KDES.PossibleZeroFill = (64:16:1000);
end
name = ['Kern_FSRw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3),'ss',num2str(SS,3),'r',num2str(res,3)];
KDES.name = regexprep(name,'\.', 'p');

Status2('done','',2);
Status2('done','',3);
