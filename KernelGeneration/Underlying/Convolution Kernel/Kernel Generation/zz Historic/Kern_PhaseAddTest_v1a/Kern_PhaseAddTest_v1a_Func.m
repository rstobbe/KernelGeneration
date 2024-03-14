%=========================================================
% 
%=========================================================

function [KDES,err] = Kern_PhaseAddTest_v1a_Func(KDES,INPUT)

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
zWtest = 2*(ceil((W*SS-2)/2)+1)/SS;
if zWtest ~= W
    error;  % look at other kernel functions, might need to deal with this
end

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
kernelsize = round(W*(1/res));
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
hFig = figure(200);
subplot(2,3,1); hold on;
plot(u0,KernProf1,'k');
plot(u0,FTSphereProf1,'c');
plot(u0,KBProf1,'g');
box on;
legend('KernProf','FTSphereProf','KBProf');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

hFig.Units = 'inches';
hFig.Position = [3 3 19 10];
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
ZF = 400;
C = ZF/2 + 1;
Status2('busy','Determine Kern Image Shape',1);
ZFtKern1 = zerofill_isotropic_odd_doubles(ifftshift(tKern1),ZF);
TF1 = fftshift(fftn(ZFtKern1));
clear ZFtKern1
TF1Prof = squeeze(TF1(C,C,:))/abs(TF1(C,C,C));
clear TF1
iCSVres = 1/((1/CSVres)+(1/W));
f = (-1/(2*iCSVres):1/(ZF*iCSVres):1/(2*iCSVres) - 1/(ZF*iCSVres));
subplot(2,3,2); hold on;
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
u = (-(W/2):res:(W/2)-res); 
Kern = zeros(ceil(W*(1/res)-1),ceil(W*(1/res)-1),ceil(W*(1/res)-1),'single');
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
sz = size(Kern);
Cen = sz(1)/2+1;
subplot(2,3,3); hold on;
plot(u,real(squeeze(Kern(:,Cen,Cen))),'r');
plot(u,imag(squeeze(Kern(:,Cen,Cen))),'b');
plot(u,abs(squeeze(Kern(:,Cen,Cen))),'k');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

%--------------------------------------------
% Generate Image
%--------------------------------------------
Im = fftshift(ifftn(ifftshift(Kern)));
MaxIm = max(abs(Im(:)));
subplot(2,3,4); hold on;
x = (-(1/res)/2:1/W:(1/res)/2-(1/W));
plot(x,real(squeeze(Im(Cen,Cen,:))),'r');
plot(x,imag(squeeze(Im(Cen,Cen,:))),'b');
plot(x,abs(squeeze(Im(Cen,Cen,:))),'k');
plot([-0.5 -0.5],[0 MaxIm],'k:');
plot([0.5 0.5],[0 MaxIm],'k:');
xlim([-(1/res)/2 (1/res)/2]);
ylabel('Kernel Image Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

%--------------------------------------------
% Add Phase
%--------------------------------------------
RelImDomain = (-1:1/W:1-1/W);
LenRelImDomain = length(RelImDomain);
%Phase = exp(1i*(pi/2)*RelImDomain).';
Phase = exp(1i*(4*pi)*RelImDomain).';
PhaseMat = repmat(Phase,1,LenRelImDomain,LenRelImDomain);
Bot = Cen-W;
Top = Cen+W-1;
Im(Bot:Top,Bot:Top,Bot:Top) = Im(Bot:Top,Bot:Top,Bot:Top) .* PhaseMat;
subplot(2,3,5); hold on;
plot(x,real(squeeze(Im(:,Cen,Cen))),'r');
plot(x,imag(squeeze(Im(:,Cen,Cen))),'b');
plot(x,abs(squeeze(Im(:,Cen,Cen))),'k');
plot([-0.5 -0.5],[0 MaxIm],'k:');
plot([0.5 0.5],[0 MaxIm],'k:');
xlim([-(1/res)/2 (1/res)/2]);
ylabel('Kernel Image Profile (Added Phase)','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Kern = fftshift(fftn(ifftshift(Im)));
subplot(2,3,6); hold on;
plot(u,real(squeeze(Kern(:,Cen,Cen))),'r');
plot(u,imag(squeeze(Kern(:,Cen,Cen))),'b');
plot(u,abs(squeeze(Kern(:,Cen,Cen))),'k');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('k-Steps (1/FoV)','fontsize',10,'fontweight','bold');

% %--------------------------------------------
% % Force Real
% %--------------------------------------------
% Kern = real(Kern);

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
KDES.iKern = 1/res;
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

name = ['Kern_PATw',num2str(W,2),'bw',num2str(BW,3),'b',num2str(beta,3),'ss',num2str(SS,3),'r',num2str(res,3)];
KDES.name = regexprep(name,'\.', 'p');

Status2('done','',2);
Status2('done','',3);
