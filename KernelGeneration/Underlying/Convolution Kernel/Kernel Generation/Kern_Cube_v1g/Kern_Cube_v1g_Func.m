%=========================================================
% 
%=========================================================

function [KDES,err] = Kern_Cube_v1g_Func(KDES,INPUT)

Status2('done','Kern_Cube',2);
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
KDES.TotalWid = KDES.CubeWid + 0.00001;
chW = ceil(((KDES.TotalWid*KDES.SubSamp)-2)/2);
if chW < 2
    KDES.TotalWid = 6/KDES.SubSamp;
    chW = ceil(((KDES.TotalWid*KDES.SubSamp)-2)/2);
    if chW ~= 2
        error;
    end
end
iKern = round(1e9*(1/(KDES.Res*KDES.SubSamp)))/1e9;            % points between subsampling
KernHalf = (chW+1)*iKern;
zW = (KernHalf*2 + 1)*KDES.Res;

%--------------------------------------------
% Test Kernel Size
%--------------------------------------------
kernelsize = (round(zW*(1/KDES.Res)*1000))/1000;
KernCen = (kernelsize+1)/2;
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
Kern = Cube;
%ConvScaleVal = sum(Kern(:))*KDES.Res^3;
ConvScaleVal = 1;

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
Status2('busy','Visualize Kernel',1);
fh = figure(1100); 
subplot(2,1,1); hold on;
plot((1:length(Kern)),squeeze(Kern(KernCen,KernCen,:))); 
xlim([1 length(Kern)]); 
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
Panel(3,:) = {'Res',KDES.Res,'Output'};
Panel(4,:) = {'SubSamp',KDES.SubSamp,'Output'};
Panel(5,:) = {'ConvScalVal',ConvScaleVal,'Output'};
KDES.Panel = Panel;
KDES.PanelOutput = cell2struct(Panel,{'label','value','type'},2);

%--------------------------------------------
% Return
%--------------------------------------------
KDES.ConvScaleVal = ConvScaleVal;
KDES.convscaleval = ConvScaleVal;
KDES.SubSamp = KDES.SubSamp;
KDES.DesforSS = KDES.SubSamp;
KDES.iKern = 1/KDES.Res;
KDES.zW = zW;
KDES.Kern = Kern;
KDES.res = KDES.Res;
KDES.W = KDES.TotalWid;
KDES.PossibleZeroFill = (64:32:1000);
name = ['Kern_FCw',num2str(KDES.CubeWid,3),'r',num2str(KDES.Res,3),'ss',num2str(KDES.SubSamp,3)];
KDES.name = regexprep(name,'\.', 'p');


Status2('done','',2);
Status2('done','',3);
