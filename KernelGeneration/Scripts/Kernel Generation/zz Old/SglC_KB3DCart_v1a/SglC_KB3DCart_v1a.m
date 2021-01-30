%==================================================
% (v1a) 
%
%==================================================

function [SCRPTipt,SCRPTGBL,err] = SglC_KB3DCart_v1a(SCRPTipt,SCRPTGBL)

err.flag = 0;
err.msg = '';

W = str2double(SCRPTGBL.CurrentScript.Width);
res = str2double(SCRPTGBL.CurrentScript.Res);
beta = str2double(SCRPTGBL.CurrentScript.Beta);
SS = str2double(SCRPTGBL.CurrentScript.DesForSS); 

KRNprms.prec = 'singles';

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
if zW/res > 1100
    err.flag = 1;
    err.msg = ['kernel size: ',num2str(zW/res)];
    return
end

%--------------------------------------------
% Generate Kaiser profile
%--------------------------------------------
Status2('busy','Generate Kernel Profile',1);
u = (0:res:(W/2));                                                 
M = beta * sqrt(1 - (2*u/W).^2);
Kern_1D = besseli(0,M);
Kern_1D = Kern_1D/Kern_1D(1);

%--------------------------------------------
% Calculate ConvScaleVal
%--------------------------------------------
Status2('busy','Calculate ConvScaleVal',1);
CSVres = 0.1;
%CSVres = 0.25;
%CSVres = 0.5; 
uT = (-(W/2):CSVres:(W/2)); 
tKern = zeros(length(uT),length(uT),length(uT));
N = length(uT);
for a = 1:N
    for b = 1:N
        for c = 1:N
            tKern(a,b,c) = Kern_1D(a)*Kern_1D(b)*Kern_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
CSV = sum(tKern(:))*CSVres^3;

%--------------------------------------------
% Generate Kernel
%--------------------------------------------
Status2('busy','Generate Kernel',1);
if strcmp(KRNprms.prec,'singles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'single');
elseif strcmp(KRNprms.prec,'doubles')
    Kern = zeros((zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,(zW/2)*round(1/res)+1,'double');
end
N = length(Kern_1D);
for a = 1:N
    for b = 1:N
        for c = 1:N
            Kern(a,b,c) = Kern_1D(a)*Kern_1D(b)*Kern_1D(c);
        end
    end
    Status2('busy',num2str(N-a),2);    
end
Status2('done','',2);

%--------------------------------------------
% Kernel Vis
%--------------------------------------------
Vis = 1;
if Vis == 1
    Status2('busy','Visualize Kernel',1);
    figure(1100); hold on; 
    plot((1:length(Kern)),squeeze(Kern(1,1,:))); 
    xlim([0 length(Kern)]); 
    xlabel('Kernel Size','fontsize',10,'fontweight','bold');
    
    figure(1101); hold on; 
    plot((0:length(Kern)-1)*res,squeeze(Kern(1,1,:))); 
    xlim([0 zW/2]); 
    xlabel('Kernel Width','fontsize',10,'fontweight','bold');
    
    figure(1102); hold on;
    anlzzf = 64000;
    KernProf = [Kern_1D zeros(1,anlzzf-2*length(Kern_1D)+1) flipdim(Kern_1D(2:length(Kern_1D)),2)];
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
end

%--------------------------------------------
% Display
%--------------------------------------------
SCRPTGBL.RWSUI.LocalOutput(1).label = 'ConvScaleVal';
SCRPTGBL.RWSUI.LocalOutput(1).value = num2str(CSV);

%--------------------------------------------
% Output
%--------------------------------------------
KRNprms.type = 'SglC_KB3DCart_v1a';
KRNprms.convscaleval = CSV;
KRNprms.DesforSS = SS;
KRNprms.iKern = 1/res;
KRNprms.W = W;
KRNprms.beta = beta;
KRNprms.res = res;
KRNprms.zW = zW;
KRNprms.Kern = Kern;
name = ['Kern_KBCw',num2str(W,2),'b',num2str(beta,3),'ss',num2str(SS,3)];
KRNprms.name = regexprep(name,'\.', 'p');

SCRPTGBL.RWSUI.SaveScriptOption = 'yes';
SCRPTGBL.RWSUI.SaveScriptPath = ['D:\1 Scripts\zs Shared\zy Convolution Kernels\',KRNprms.name];
SCRPTGBL.RWSUI.SaveVariables = {KRNprms};
SCRPTGBL.RWSUI.SaveVariableNames = {'KRNprms'};


