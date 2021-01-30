%=========================================================
% 
%=========================================================

function [KDES,err] = SglC_KB3DCart_v1e_Func(KDES,INPUT)

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

%---------------------------------------------
% Find all SS values that are compatible with res
%---------------------------------------------
m = 1;
for n = 1:0.0005:6
    test = round(1e9/(res*n))/1e9;
    if not(rem(test,1))
        posSS(m) = n;
        m = m+1;
    end    
end

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
for n = 1:length(posSS)
    zWtest = 2*(ceil((W*posSS(n)-2)/2)+1)/posSS(n);
    zW = W + 2*res;
    while true
        if zW > zWtest
            break
        end
        zW = zW + 2*res;
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
parfor a = 1:N
    for b = 1:N
        for c = 1:N
            tKern(a,b,c) = CSVKaiserProf(a)*CSVKaiserProf(b)*CSVKaiserProf(c);
        end
    end 
end
Status2('done','',2);
CSV = sum(tKern(:))*CSVres^3;

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
if strcmp(KDES.prec,'singles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
elseif strcmp(KDES.prec,'doubles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'double');
end
N = length(KaiserProf);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = KaiserProf(a)*KaiserProf(b)*KaiserProf(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
Status2('busy','Visualize Kernel',1);
fh = figure(1100); 
subplot(2,2,1); hold on;
plot((1:length(Kern)),squeeze(Kern(1,1,:))); 
xlim([0 length(Kern)]); 
xlabel('Kernel Size','fontsize',10,'fontweight','bold');

subplot(2,2,2); hold on;
plot((0:length(Kern)-1)*res,squeeze(Kern(1,1,:))); 
xlim([0 zW/2]); 
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
KDES.posSS = posSS;
KDES.iKern = 1/res;
KDES.zW = zW;
KDES.Kern = Kern;
name = ['Kern_KBCw',num2str(W,2),'b',num2str(beta,3),'ss',num2str(SS,3)];
KDES.name = regexprep(name,'\.', 'p');

Status2('done','',2);
Status2('done','',3);
