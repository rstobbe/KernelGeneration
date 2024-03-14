%=========================================================
% 
%=========================================================

function [KDES,err] = SglC_KB3DCartFullComplexTest_v1g_Func(KDES,INPUT)

Status2('done','SglC_KB3DCart Kernel',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Get Input
%---------------------------------------------
W = KDES.W;
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

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
chW = ceil(((KDES.W*KDES.DesforSS)-2)/2);  
iKern = round(1e9*(1/(KDES.res*KDES.DesforSS)))/1e9;            % points between subsampling
KernCen = (chW+1)*iKern;
zW = (KernCen*2 + 1)*res;

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
kernelsize = ceil(zW*(1/res));
button = questdlg(['Kernel Dimensions: ',num2str(kernelsize),'.  Continue?']);
if strcmp(button,'No') || strcmp(button,'Cancel')
    err.flag = 4;
    err.msg = '';
    return
end
if kernelsize > 1000
    err.flag = 1;
    err.msg = 'excessive kernlsize';
    return
end

%--------------------------------------------
% Generate Kaiser profile
%--------------------------------------------
Status2('busy','Generate Kernel Profile',1);
u = (0:res:(W/2));                                                 
KaiserProf = besseli(0,beta*sqrt(1-(2*u/W).^2));
KaiserProfScale = KaiserProf(1);
KaiserProf = KaiserProf/KaiserProfScale;

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Calculate ConvScaleVal',1);
CSVres = res;
if rem((W/2),CSVres)
    error();
end
CSVu = (-(W/2):CSVres:(W/2));
CSVKaiserProf = besseli(0,beta*sqrt(1-(2*CSVu/W).^2));
CSVKaiserProf = CSVKaiserProf/KaiserProfScale;
tKern = zeros(length(CSVu),length(CSVu),length(CSVu));
N = length(CSVu);
Status2('busy','Workers in Parallel',2); 
for a = 1:N
    for b = 1:N
        for c = 1:N
            tKern(a,b,c) = CSVKaiserProf(a)*CSVKaiserProf(b)*CSVKaiserProf(c);
        end
    end 
end
Status2('done','',2);
CSV = sum(tKern(:))*CSVres^3;

%--------------------------------------------
% Generate Kaiser profile full
%--------------------------------------------
Status2('busy','Generate Kernel Profile Full',1);
u = (-(W/2):res:(W/2));                                                 
KaiserProfFull = besseli(0,beta*sqrt(1-(2*u/W).^2));
KaiserProfFull = KaiserProfFull/KaiserProfScale;
ZeroAdd = (kernelsize-length(KaiserProfFull))/2;
KaiserProfFull = [zeros(1,ZeroAdd),KaiserProfFull,zeros(1,ZeroAdd)];

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
Kern = zeros(ceil(zW*(1/res)),ceil(zW*(1/res)),ceil(zW*(1/res)),'single');
N = length(KaiserProfFull);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = KaiserProfFull(a)*KaiserProfFull(b)*KaiserProfFull(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);
Zeros = zeros(ceil(zW*(1/res)),ceil(zW*(1/res)),ceil(zW*(1/res)),'single');
Kern = complex(Kern,Zeros);
test = isreal(Kern)

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
KernCen = (N+1)/2;
Status2('busy','Visualize Kernel',1);
fh = figure(1100); 
subplot(2,2,1); hold on;
plot((1:length(Kern)),squeeze(Kern(KernCen,KernCen,:))); 
xlim([0 length(Kern)]); 
xlabel('Kernel Size','fontsize',10,'fontweight','bold');

subplot(2,2,2); hold on;
hW = (length(Kern)-1)/2;
plot((-hW:hW)*res,squeeze(Kern(KernCen,KernCen,:))); 
xlim([-zW/2 zW/2]); 
xlabel('Kernel Width','fontsize',10,'fontweight','bold');

subplot(2,2,3); hold on;
anlzzf = 64000;
KernProf = [KaiserProf zeros(1,anlzzf-2*length(KaiserProf)+1) flip(KaiserProf(2:length(KaiserProf)),2)];
FTProf = ifftshift(ifft(KernProf));
FTProf = FTProf/max(FTProf);
FOV = 1/(2*res);                                                        
t = (-FOV:(1/res)/anlzzf:FOV-(1/res)/anlzzf); 
plot(t,real(FTProf),'r');
plot(t,imag(FTProf),'b'); 
plot([0.5 0.5],[-0.2 1.2],'k:');
plot([-0.5 -0.5],[-0.2 1.2],'k:');
plot([SS/2 SS/2],[-0.2 1.2],'k:');
plot([-SS/2 -SS/2],[-0.2 1.2],'k:');
ylabel('Kernel Profile','fontsize',10,'fontweight','bold');
xlabel('Relative Image Space','fontsize',10,'fontweight','bold');
set(gca,'XLim',[-(SS/2+0.1) (SS/2+0.1)]);

fh.Name = 'Kernel Characteristics';
fh.NumberTitle = 'off';
fh.Position = [500 200 1000 800];

KDES.Figure(1).Name = 'Kernel Characteristics';
KDES.Figure(1).Type = 'Graph';
KDES.Figure(1).hFig = fh;
KDES.Figure(1).hAx = gca;

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
% Return
%--------------------------------------------
KDES.convscaleval = CSV;
KDES.DesforSS = SS;
KDES.iKern = 1/res;
KDES.zW = zW;
KDES.Kern = Kern;
KDES.PossibleZeroFill = (64:32:1000);
name = ['Kern_KBCFw',num2str(W,2),'b',num2str(beta,3),'ss',num2str(SS,3),'r',num2str(res,3),'_ComplexZeroTest'];
KDES.name = regexprep(name,'\.', 'p');


Status2('done','',2);
Status2('done','',3);
