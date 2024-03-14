%=========================================================
% 
%=========================================================

function [KDES,err] = Kern_FiltCube_v1g_Func(KDES,INPUT)

Status2('done','Kern_FiltCube',2);
Status2('done','',3);

err.flag = 0;
err.msg = '';

%---------------------------------------------
% Tests
%---------------------------------------------
if rem((KDES.CubeWid/2),KDES.Res)
    err.flag = 1;
    err.msg = 'CubeWid not a multiple of Res';
    return
end
if rem(round(1e9*(1/(KDES.Res*KDES.SubSamp)))/1e9,1)
    err.flag = 1;
    err.msg = '1/(Res*SubSamp) not an integer';
    return
end

%--------------------------------------------
% Determine How much to zero-fill for C-alg
%--------------------------------------------
KDES.TotalWid = KDES.CubeWid + KDES.FiltWid;
chW = ceil(((KDES.TotalWid*KDES.SubSamp)-2)/2);  
iKern = round(1e9*(1/(KDES.Res*KDES.SubSamp)))/1e9;            % points between subsampling
KernHalf = (chW+1)*iKern;
zW = (KernHalf*2 + 1)*KDES.Res;

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
kernelsize = ceil(zW*(1/KDES.Res));
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
% Generate Filter Profile
%--------------------------------------------
Status2('busy','Generate Filter Profile Full',1);
u = (-(KDES.FiltWid/2):KDES.Res:(KDES.FiltWid/2));                                                 
KaiserProfFull = besseli(0,KDES.FiltBeta*sqrt(1-(2*u/KDES.FiltWid).^2));
KaiserProfScale = max(KaiserProfFull(:));
KaiserProfFull = KaiserProfFull/KaiserProfScale;
ZeroAdd = (kernelsize-length(KaiserProfFull))/2;
KaiserProfFull = [zeros(1,ZeroAdd),KaiserProfFull,zeros(1,ZeroAdd)];

%--------------------------------------------
% Generate Filter
%--------------------------------------------
Status2('busy','Generate Filter',1);
Filt = zeros(ceil(zW*(1/KDES.Res)),ceil(zW*(1/KDES.Res)),ceil(zW*(1/KDES.Res)),'single');
N = length(KaiserProfFull);
KernCen = (N+1)/2;
for a = 1:N
    for b = 1:N
        for c = 1:N
            Filt(a,b,c) = KaiserProfFull(a)*KaiserProfFull(b)*KaiserProfFull(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%--------------------------------------------
% Generate Cube
%--------------------------------------------
Status2('busy','Generate Cube',1);
Cube = zeros(ceil(zW*(1/KDES.Res)),ceil(zW*(1/KDES.Res)),ceil(zW*(1/KDES.Res)),'single');
Wid = KDES.CubeWid/KDES.Res;
Bot = KernCen-Wid/2;
Top = Bot+Wid;
Cube(Bot:Top,Bot:Top,Bot:Top) = 1;

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
Kern = convn(Cube,Filt,'same');
Kern = Kern/max(Kern(:));
ConvScaleVal = sum(Kern(:))*KDES.Res^3;

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
Status2('busy','Visualize Kernel',1);
fh = figure(1100); 
subplot(2,1,1); hold on;
plot((1:length(Kern)),squeeze(Cube(KernCen,KernCen,:))); 
plot((1:length(Kern)),squeeze(Filt(KernCen,KernCen,:))); 
plot((1:length(Kern)),squeeze(Kern(KernCen,KernCen,:))); 
xlim([0 length(Kern)]); 
xlabel('Kernel Size','fontsize',10,'fontweight','bold');

subplot(2,1,2); hold on;
hW = (length(Kern)-1)/2;
plot((-hW:hW)*KDES.Res,squeeze(Kern(KernCen,KernCen,:))); 
xlim([-zW/2 zW/2]); 
xlabel('Kernel Width','fontsize',10,'fontweight','bold');

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
Panel(2,:) = {'CubeWidth',KDES.CubeWid,'Output'};
Panel(3,:) = {'FiltBeta',KDES.FiltBeta,'Output'};
Panel(4,:) = {'Res',KDES.Res,'Output'};
Panel(5,:) = {'SubSamp',KDES.SubSamp,'Output'};
Panel(6,:) = {'ConvScalVal',ConvScaleVal,'Output'};
KDES.Panel = Panel;
KDES.PanelOutput = cell2struct(Panel,{'label','value','type'},2);

%--------------------------------------------
% Return
%--------------------------------------------
KDES.ConvScaleVal = sum(Kern(:))*KDES.Res^3;
KDES.convscaleval = KDES.ConvScaleVal;
KDES.SubSamp = KDES.SubSamp;
KDES.DesforSS = KDES.SubSamp;
KDES.iKern = 1/KDES.Res;
KDES.zW = zW;
KDES.Kern = Kern;
KDES.res = KDES.Res;
KDES.W = KDES.TotalWid;
KDES.PossibleZeroFill = (64:32:1000);
name = ['Kern_FCw',num2str(KDES.CubeWid,2),'r',num2str(KDES.Res,3),'w',num2str(KDES.FiltWid,2),'b',num2str(KDES.FiltBeta,3),'ss',num2str(KDES.SubSamp,3)];
KDES.name = regexprep(name,'\.', 'p');


Status2('done','',2);
Status2('done','',3);
